#include "hlbvh.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

// Runtime compilation keeps installation independent of the optional Xcode
// Metal command-line toolchain. All resources are ARC-owned; autoreleased
// command objects are drained at the end of each build.
namespace {
const char *kernels = R"MSL(
#include <metal_stdlib>
using namespace metal;
struct Bounds { float4 lo, hi; };
struct Morton { uint code, index; };
struct Treelet { uint first, count; };
struct Node { Bounds bounds; uint left, right, first, count, axis, padding[3]; };
struct Params { Bounds centers; uint n, blocks, shift, max_leaf; };
uint spread(uint v) {
    v=(v|(v<<16))&0x030000FF; v=(v|(v<<8))&0x0300F00F;
    v=(v|(v<<4))&0x030C30C3; return (v|(v<<2))&0x09249249;
}
Bounds empty_box() { return {float4(INFINITY), float4(-INFINITY)}; }
Bounds unite(Bounds a, Bounds b) { return {min(a.lo,b.lo), max(a.hi,b.hi)}; }
kernel void morton(device const Bounds* bounds [[buffer(0)]], device Morton* order [[buffer(1)]],
                   constant Params& p [[buffer(2)]], uint i [[thread_position_in_grid]]) {
    if(i>=p.n) return;
    uint code=0;
    for(uint d=0;d<3;++d) {
        float c=.5f*bounds[i].lo[d]+.5f*bounds[i].hi[d];
        float extent=.5f*p.centers.hi[d]-.5f*p.centers.lo[d];
        float offset=extent>0 ? (.5f*c-.5f*p.centers.lo[d])/extent : 0;
        code |= spread(uint(clamp(offset*1024.f,0.f,1023.f)))<<d;
    }
    order[i]={code,i};
}
// Stable radix sort with 256 records per cooperative threadgroup. Each SIMD
// group counts/ranks its digits; shared totals combine SIMD groups in input
// order. Reads and writes are per record rather than 64 serial records per lane.
kernel void histogram(device const Morton* order [[buffer(0)]], device uint* counts [[buffer(1)]],
                      constant Params& p [[buffer(2)]], uint block [[threadgroup_position_in_grid]],
                      uint lane [[thread_index_in_threadgroup]], uint sg [[simdgroup_index_in_threadgroup]],
                      uint sl [[thread_index_in_simdgroup]], uint width [[threads_per_simdgroup]]) {
    threadgroup uint hist[16*32];
    uint i=block*256+lane, groups=256/width;
    uint digit=i<p.n ? (order[i].code>>p.shift)&15 : 16;
    for(uint b=0;b<16;++b) {
        uint sum=simd_sum(uint(digit==b));
        if(sl==0) hist[b*32+sg]=sum;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(lane<16) {
        uint sum=0; for(uint g=0;g<groups;++g) sum+=hist[lane*32+g];
        counts[lane*p.blocks+block]=sum;
    }
}
// Scan each 256-block chunk in parallel. Only the much smaller chunk totals
// reach the final 16-lane scan; serial work is reduced by another factor of 256.
kernel void scan_chunks(device uint* counts [[buffer(0)]], device uint* chunks [[buffer(1)]],
                        constant Params& p [[buffer(2)]], uint group [[threadgroup_position_in_grid]],
                        uint lane [[thread_index_in_threadgroup]], uint sg [[simdgroup_index_in_threadgroup]],
                        uint sl [[thread_index_in_simdgroup]], uint width [[threads_per_simdgroup]]) {
    threadgroup uint sums[32];
    uint nchunks=(p.blocks+255)/256, b=group/nchunks, chunk=group%nchunks;
    uint i=chunk*256+lane;
    uint value=i<p.blocks ? counts[b*p.blocks+i] : 0;
    uint prefix=simd_prefix_exclusive_sum(value), sum=simd_sum(value);
    if(sl==0) sums[sg]=sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for(uint g=0;g<sg;++g) prefix+=sums[g];
    if(i<p.blocks) counts[b*p.blocks+i]=prefix;
    if(lane==255) chunks[b*nchunks+chunk]=prefix+value;
}
kernel void totals(device const uint* counts [[buffer(0)]], device uint* sums [[buffer(1)]],
                   constant Params& p [[buffer(2)]], uint b [[thread_position_in_grid]]) {
    if(b>=16) return;
    uint sum=0; for(uint i=0;i<p.blocks;++i) sum+=counts[b*p.blocks+i]; sums[b]=sum;
}
kernel void prefix(device uint* counts [[buffer(0)]], device const uint* sums [[buffer(1)]],
                   constant Params& p [[buffer(2)]], uint b [[thread_position_in_grid]]) {
    if(b>=16) return;
    uint base=0; for(uint i=0;i<b;++i) base+=sums[i];
    for(uint i=0;i<p.blocks;++i) { uint n=counts[b*p.blocks+i]; counts[b*p.blocks+i]=base; base+=n; }
}
kernel void scatter(device const Morton* input [[buffer(0)]], device Morton* output [[buffer(1)]],
                    device const uint* offsets [[buffer(2)]], device const uint* chunks [[buffer(3)]],
                    constant Params& p [[buffer(4)]], uint block [[threadgroup_position_in_grid]],
                    uint lane [[thread_index_in_threadgroup]], uint sg [[simdgroup_index_in_threadgroup]],
                    uint sl [[thread_index_in_simdgroup]], uint width [[threads_per_simdgroup]]) {
    threadgroup uint hist[16*32];
    uint i=block*256+lane;
    Morton m={0,0}; if(i<p.n) m=input[i];
    uint digit=i<p.n ? (m.code>>p.shift)&15 : 16, rank=0;
    for(uint b=0;b<16;++b) {
        uint local=simd_prefix_exclusive_sum(uint(digit==b));
        uint sum=simd_sum(uint(digit==b));
        if(digit==b) rank=local;
        if(sl==0) hist[b*32+sg]=sum;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if(i<p.n) {
        for(uint g=0;g<sg;++g) rank+=hist[digit*32+g];
        uint nchunks=(p.blocks+255)/256;
        output[offsets[digit*p.blocks+block]+chunks[digit*nchunks+block/256]+rank]=m;
    }
}
// The high 12 Morton bits define at most 4096 cells. Each lane finds its
// cell's sorted interval independently. Empty cells have count zero. Keeping
// this grouping on the GPU lets topology follow sorting without a CPU wait.
kernel void group_cells(device const Morton* order [[buffer(0)]], device Treelet* cells [[buffer(1)]],
                        constant Params& p [[buffer(2)]], uint cell [[thread_position_in_grid]]) {
    if(cell>=4096) return;
    uint lo=0, hi=p.n;
    while(lo<hi) {
        uint mid=lo+(hi-lo)/2;
        if((order[mid].code>>18)<cell) lo=mid+1; else hi=mid;
    }
    uint first=lo;
    hi=p.n;
    while(lo<hi) {
        uint mid=lo+(hi-lo)/2;
        if((order[mid].code>>18)<=cell) lo=mid+1; else hi=mid;
    }
    cells[cell]={first,lo-first};
}
// The renderer uses singleton leaves. Give every sorted primitive a lane and
// descend its Morton path independently. Only the leftmost primitive emits each
// internal node, so writes are disjoint without queues or atomic publication.
// A full binary subtree with L leaves occupies 2*L-1 preorder slots: this lets
// every lane derive child addresses, preserving the CPU builder's exact layout.
// Repeated Morton codes use the same balanced index split as the CPU.
kernel void parallel_topology(device const Bounds* bounds [[buffer(0)]], device const Morton* order [[buffer(1)]],
                              device const Treelet* cells [[buffer(2)]], device Node* nodes [[buffer(3)]],
                              constant Params& p [[buffer(4)]], uint i [[thread_position_in_grid]]) {
    if(i>=p.n) return;
    Treelet t=cells[order[i].code>>18];
    uint first=t.first, count=t.count, node=2*first;
    int bit=17;
    while(count>1) {
        while(bit>=0 && ((order[first].code^order[first+count-1].code)&(1u<<bit))==0) --bit;
        uint split=first+count/2;
        if(bit>=0) {
            uint lo=first, hi=first+count;
            while(lo<hi) { uint m=lo+(hi-lo)/2; if(order[m].code&(1u<<bit)) hi=m; else lo=m+1; }
            split=lo;
        }
        uint right=node+2*(split-first);
        if(i==first) {
            nodes[node].left=node+1; nodes[node].right=right;
            nodes[node].first=first; nodes[node].count=0;
            nodes[node].axis=bit>=0 ? uint(bit%3) : 0;
            nodes[node].padding[0]=count;
        }
        if(i<split) { ++node; count=split-first; }
        else { node=right; count=first+count-split; first=split; }
        --bit;
    }
    nodes[node].bounds=bounds[order[i].index];
    nodes[node].first=i; nodes[node].count=1; nodes[node].axis=0;
}
// Each SIMD group reduces one internal node's primitive interval. This repeats
// bound reads across ancestors, but exposes all nodes concurrently and avoids
// cross-threadgroup acquire/release assumptions or a per-level dispatch chain.
// Only min/max are used, so reduction order leaves the final bounds unchanged.
kernel void parallel_bounds(device const Bounds* bounds [[buffer(0)]], device const Morton* order [[buffer(1)]],
                            device const Treelet* cells [[buffer(2)]], device Node* nodes [[buffer(3)]],
                            constant Params& p [[buffer(4)]], uint group [[threadgroup_position_in_grid]],
                            uint sg [[simdgroup_index_in_threadgroup]], uint sl [[thread_index_in_simdgroup]],
                            uint width [[threads_per_simdgroup]]) {
    uint i=group*(128/width)+sg;
    if(i>=2*p.n) return;
    Treelet t=cells[order[i/2].code>>18];
    // Each independent treelet arena reserves one unused slot at its end.
    if(i==2*(t.first+t.count)-1 || nodes[i].count!=0) return;
    uint first=nodes[i].first, count=nodes[i].padding[0];
    float3 lo=float3(INFINITY), hi=float3(-INFINITY);
    for(uint j=sl;j<count;j+=width) {
        Bounds b=bounds[order[first+j].index];
        lo=min(lo,b.lo.xyz); hi=max(hi,b.hi.xyz);
    }
    lo=float3(simd_min(lo.x),simd_min(lo.y),simd_min(lo.z));
    hi=float3(simd_max(hi.x),simd_max(hi.y),simd_max(hi.z));
    if(sl==0) nodes[i].bounds={float4(lo,0),float4(hi,0)};
}

// General multi-primitive leaves retain the stack builder below.
// One lane per occupied Morton cell. Topology and postorder bound propagation
// both run here, not merely Morton-code generation. No inter-lane atomics or
// device-wide synchronization are needed between independent treelets.
kernel void treelets(device const Bounds* bounds [[buffer(0)]], device const Morton* order [[buffer(1)]],
                     device const Treelet* trees [[buffer(2)]], device Node* nodes [[buffer(3)]],
                     constant Params& p [[buffer(4)]], uint i [[thread_position_in_grid]]) {
    if(i>=p.blocks) return;
    struct Frame { uint node, first, count, split, stage; int bit; };
    // At most 18 Morton levels plus 30 balanced index levels under the host's
    // signed node-index limit, so this stack cannot overflow for valid inputs.
    Frame stack[64];
    Treelet t=trees[i];
    if(t.count==0) return;
    uint next=2*t.first;
    int top=0; stack[0]={next++,t.first,t.count,0,0,17};
    while(top>=0) {
        thread Frame& f=stack[top];
        if(f.stage==0) {
            if(f.count<=p.max_leaf) {
                Bounds box=empty_box();
                for(uint j=f.first;j<f.first+f.count;++j) box=unite(box,bounds[order[j].index]);
                nodes[f.node].bounds=box; nodes[f.node].first=f.first;
                nodes[f.node].count=f.count; nodes[f.node].axis=0; --top; continue;
            }
            while(f.bit>=0 && ((order[f.first].code^order[f.first+f.count-1].code)&(1u<<f.bit))==0) --f.bit;
            f.split=f.first+f.count/2;
            if(f.bit>=0) {
                uint lo=f.first, hi=f.first+f.count;
                while(lo<hi) { uint m=lo+(hi-lo)/2; if(order[m].code&(1u<<f.bit)) hi=m; else lo=m+1; }
                f.split=lo;
            }
            nodes[f.node].count=0; nodes[f.node].axis=f.bit>=0 ? uint(f.bit%3) : 0;
            nodes[f.node].left=next;
            Frame child={next++,f.first,f.split-f.first,0,0,f.bit-1};
            f.stage=1; stack[++top]=child;
        } else if(f.stage==1) {
            nodes[f.node].right=next;
            Frame child={next++,f.split,f.first+f.count-f.split,0,0,f.bit-1};
            f.stage=2; stack[++top]=child;
        } else {
            nodes[f.node].bounds=unite(nodes[nodes[f.node].left].bounds,nodes[nodes[f.node].right].bounds);
            --top;
        }
    }
}
)MSL";

struct MetalBuilder {
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
  id<MTLComputePipelineState> morton, histogram, scan_chunks, totals, prefix, scatter, group_cells, treelets, parallel_topology,
      parallel_bounds;
  MetalBuilder() {
    device = MTLCreateSystemDefaultDevice();
    if (!device)
      return;
    queue = [device newCommandQueue];
    NSError *error = nil;
    MTLCompileOptions *options = [MTLCompileOptions new];
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
    if (@available(macOS 15.0, *))
      options.mathMode = MTLMathModeSafe;
    else
#endif
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      options.fastMathEnabled = NO;
#pragma clang diagnostic pop
    }
    id<MTLLibrary> library = [device newLibraryWithSource:[NSString stringWithUTF8String:kernels]
                                                  options:options
                                                    error:&error];
    if (!library)
      throw std::runtime_error("Metal BVH shader compilation: " + std::string(error.localizedDescription.UTF8String));
    auto pipeline = [&](NSString *name) {
      NSError *failure = nil;
      id<MTLComputePipelineState> result =
          [device newComputePipelineStateWithFunction:[library newFunctionWithName:name] error:&failure];
      if (!result)
        throw std::runtime_error("Metal BVH pipeline: " + std::string(failure.localizedDescription.UTF8String));
      return result;
    };
    morton = pipeline(@"morton");
    histogram = pipeline(@"histogram");
    scan_chunks = pipeline(@"scan_chunks");
    totals = pipeline(@"totals");
    prefix = pipeline(@"prefix");
    scatter = pipeline(@"scatter");
    group_cells = pipeline(@"group_cells");
    treelets = pipeline(@"treelets");
    parallel_topology = pipeline(@"parallel_topology");
    parallel_bounds = pipeline(@"parallel_bounds");
    for (auto state : {histogram, scan_chunks, scatter, parallel_bounds}) {
      const NSUInteger width = state.threadExecutionWidth;
      if (width < 8 || width > 128 || 128 % width != 0)
        throw std::runtime_error("Metal BVH requires a supported SIMD group width; use hlbvh.");
    }
  }
};
MetalBuilder &builder() {
  static MetalBuilder state;
  return state;
}
struct Params {
  HLBVHBounds centers;
  uint32_t n, blocks, shift, max_leaf;
};
static_assert(sizeof(Params) == 48);
void finish(id<MTLCommandBuffer> command) {
  [command commit];
  [command waitUntilCompleted];
  if (command.status == MTLCommandBufferStatusError)
    throw std::runtime_error("Metal BVH command failed: " + std::string(command.error.localizedDescription.UTF8String));
}
} // namespace

bool MetalBVHAvailable() {
  @autoreleasepool {
    return MTLCreateSystemDefaultDevice() != nil;
  }
}

void BuildMetalTreelets(std::span<const HLBVHBounds> bounds, const HLBVHBounds &centroids, unsigned max_leaf,
                        HLBVHTree &result) {
  @autoreleasepool {
    const bool profile = std::getenv("RAYRENDER_METAL_BVH_PROFILE") != nullptr;
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    auto &gpu = builder();

    if (!gpu.device || !gpu.queue)
      throw std::runtime_error("No Metal device is available; use bvh_type = 'hlbvh'.");
    auto buffer = [&](size_t bytes) {
      id<MTLBuffer> value = [gpu.device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
      if (!value)
        throw std::runtime_error("Metal BVH buffer allocation failed; use bvh_type = 'hlbvh'.");
      return value;
    };
    Params p{centroids, uint32_t(bounds.size()), uint32_t((bounds.size() + 255) / 256), 0, max_leaf};
    id<MTLBuffer> boxes = buffer(bounds.size_bytes()), a = buffer(bounds.size() * sizeof(MortonPrimitive));
    id<MTLBuffer> b = buffer(a.length), counts = buffer(size_t(p.blocks) * 16 * sizeof(uint32_t)),
                  sums = buffer(16 * sizeof(uint32_t)),
                  chunks = buffer(size_t((p.blocks + 255) / 256) * 16 * sizeof(uint32_t));
    id<MTLBuffer> cells = buffer(4096 * sizeof(HLBVHTreelet));
    id<MTLBuffer> nodes = buffer((2 * bounds.size() + 4096) * sizeof(HLBVHNode));
    std::memcpy(boxes.contents, bounds.data(), bounds.size_bytes());
    id<MTLCommandBuffer> command = [gpu.queue commandBuffer];
    command.label = @"HLBVH sort, group, and build";
    auto dispatch = [&](id<MTLComputePipelineState> pipeline, NSArray<id<MTLBuffer>> *buffers, NSUInteger n,
                        bool cooperative = false) {
      if (!command)
        throw std::runtime_error("Metal BVH could not allocate a command buffer.");
      id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
      if (!encoder)
        throw std::runtime_error("Metal BVH could not allocate a command encoder.");
      [encoder setComputePipelineState:pipeline];
      for (NSUInteger i = 0; i < buffers.count; ++i)
        [encoder setBuffer:buffers[i] offset:0 atIndex:i];
      [encoder setBytes:&p length:sizeof(p) atIndex:buffers.count];
      NSUInteger group = cooperative ? 256
                         : pipeline == gpu.parallel_bounds
                             ? 128
                             : std::min(NSUInteger(128), pipeline.maxTotalThreadsPerThreadgroup);
      if (group > pipeline.maxTotalThreadsPerThreadgroup)
        throw std::runtime_error("Metal BVH kernel threadgroup exceeds the device pipeline limit; use hlbvh.");
      [encoder dispatchThreadgroups:MTLSizeMake((n + group - 1) / group, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(group, 1, 1)];
      [encoder endEncoding];
    };
    const auto ready = Clock::now();
    dispatch(gpu.morton, @[ boxes, a ], p.n);
    for (p.shift = 0; p.shift < 30; p.shift += 4) {
      const uint32_t blocks = p.blocks, nchunks = (blocks + 255) / 256;
      dispatch(gpu.histogram, @[ a, counts ], blocks * 256, true);
      dispatch(gpu.scan_chunks, @[ counts, chunks ], 16 * nchunks * 256, true);
      p.blocks = nchunks;
      dispatch(gpu.totals, @[ chunks, sums ], 16);
      dispatch(gpu.prefix, @[ chunks, sums ], 16);
      p.blocks = blocks;
      dispatch(gpu.scatter, @[ a, b, counts, chunks ], blocks * 256, true);
      std::swap(a, b);
    }
    // Separate encoders with tracked buffers preserve producer/consumer order.
    // Submit once: no host readback or allocation separates sort and topology.
    dispatch(gpu.group_cells, @[ a, cells ], 4096);
    p.blocks = 4096;
    if (max_leaf == 1) {
      dispatch(gpu.parallel_topology, @[ boxes, a, cells, nodes ], p.n);
      dispatch(gpu.parallel_bounds, @[ boxes, a, cells, nodes ],
               NSUInteger(2) * p.n * gpu.parallel_bounds.threadExecutionWidth);
    } else {
      dispatch(gpu.treelets, @[ boxes, a, cells, nodes ], p.blocks);
    }

    finish(command);
    const auto built = Clock::now();
    const double build_gpu = profile ? command.GPUEndTime - command.GPUStartTime : 0;
    std::memcpy(result.order.data(), a.contents, a.length);
    const auto *cell_table = static_cast<const HLBVHTreelet *>(cells.contents);
    for (unsigned i = 0; i < 4096; ++i)
      if (cell_table[i].count)
        result.treelets.push_back(cell_table[i]);
    // Metal shared storage is CPU-addressable after command completion. Retain
    // the ARC-owned buffer in the deleter instead of copying its entire arena.
    // Extra slots above 2*N are reserved for host upper-SAH stitching.
    result.nodes =
        std::shared_ptr<HLBVHNode[]>(static_cast<HLBVHNode *>(nodes.contents), [nodes](HLBVHNode *) { (void)nodes; });

    if (profile) {
      const auto copied = Clock::now();
      auto seconds = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
      std::fprintf(stderr,
                   "MetalBVH n=%u prepare=%.6f build_wall=%.6f build_gpu=%.6f readback_retain=%.6f\n",
                   p.n, seconds(start, ready), seconds(ready, built), build_gpu, seconds(built, copied));
    }
  }
}
