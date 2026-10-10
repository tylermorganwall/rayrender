R"MSL(

  struct TextureProbe {
    float4 uv, position, object_position, normal, object_normal, footprint;
  };
  kernel void probe_textures(device const TextureProbe *inputs [[buffer(0)]],
                             device float4 *outputs [[buffer(1)]],
                             device const Texture *textures [[buffer(2)]],
                             device const float4 *texels [[buffer(3)]],
                             uint i [[thread_position_in_grid]]) {
    TextureProbe p = inputs[i];
    TextureContext h = {p.uv.xy,      p.position.xyz,      p.object_position.xyz,
                        p.normal.xyz, p.object_normal.xyz, p.footprint};
    outputs[i] = float4(texture_value(uint(p.uv.z), h, textures, texels), 0);
  }
  struct BsdfProbe {
    Material material;
    float4 rho, transmission, view, direction, random, context;
  };
  struct BsdfResult {
    float4 evaluated, direction, weight;
  };
  kernel void probe_materials(device const BsdfProbe *inputs [[buffer(0)]],
                              device BsdfResult *outputs [[buffer(1)]],
                              device const float4 *data [[buffer(2)]],
                              uint i [[thread_position_in_grid]]) {
    BsdfProbe probe = inputs[i];
    BsdfValue value;
    BsdfSample sample;
    if (probe.material.type == 0) {
      auto mapped = mapped_diffuse(float3(0, 0, 1), probe.context.xyz, probe.view.xyz);
      value = mapped_eval(mapped, probe.direction.xyz, probe.rho.xyz, probe.material);
      float3 wi = mapped_sample(mapped, probe.random.xy, probe.random.z, probe.random.w);
      sample = {};
      sample.direction = dot(mapped.g, wi) > 0 ? wi : float3(0);
    } else if (probe.material.type == bsdf_openpbr) {
      PbrState b = prepare_pbr(probe.material, probe.rho.xyz, probe.view.xyz, float3(0, 0, 1),
                               float3(0, 0, 1), float3(1, 0, 0), probe.context.z, -1, data);
      value = pbr_eval(b, probe.direction.xyz);
      sample = pbr_sample(b, probe.random.xyz);
    } else {
      SurfaceBsdf b;
      b.material = probe.material;
      b.frame = {float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)};
      b.view = probe.view.xyz;
      b.rho = probe.rho.xyz;
      b.transmission = probe.transmission.xyz;
      b.distance = probe.context.x;
      b.impact = probe.context.y;
      value = surface_eval(b, probe.direction.xyz);
      sample = surface_sample(b, probe.random);
    }
    outputs[i] = {float4(value.value, value.pdf), float4(sample.direction, sample.pdf),
                  float4(sample.weight, sample.delta ? 1.f : 0.f)};
  }

)MSL"
R"MSL(

  struct GeometryProbe {
    float4 origin, direction;
    uint body, padding[3];
  };
  struct GeometryResult {
    float4 point, normal, uv;
    uint primitive, material, boundary, padding;
  };
  kernel void probe_geometry(device const Triangle *triangles [[buffer(0)]],
                             device const GeometryProbe *input [[buffer(1)]],
                             device GeometryResult *output [[buffer(2)]],
                             primitive_acceleration_structure triangle_scene [[buffer(16)]],
                             intersection_function_table<triangle_data> functions [[buffer(25)]],
                             device const Surface *surfaces [[buffer(26)]],
                             device const Quadric *quadrics [[buffer(28)]],
                             instance_acceleration_structure mixed_scene [[buffer(29)]],
                             intersection_function_table<triangle_data, instancing> mixed_functions
                             [[buffer(30)]],
                             uint i [[thread_position_in_grid]]) {
    GeometryScene scene = {triangle_scene, functions, mixed_scene, mixed_functions};
    GeometryProbe p = input[i];
    ray r;
    r.origin = p.origin.xyz;
    r.direction = p.direction.xyz;
    r.min_distance = p.origin.w;
    r.max_distance = p.direction.w;
    auto hit = trace_geometry(r, scene, quadrics, p.body);
    GeometryResult result = {};
    result.primitive = UINT_MAX;
    if (hit.type != intersection_type::none) {
      auto surface =
          geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
      result.point = float4(surface.position, hit.distance);
      result.normal = float4(surface.normal, 1);
      float3 bary = float3(1 - hit.coordinates.x - hit.coordinates.y, hit.coordinates.xy);
      result.uv.xy = surface_context(surface.face, bary, surface.position, surfaces).uv;
      result.primitive = hit.primitive_id;
      result.material = surface.face.material;
      result.boundary = surface.face.boundary;
    }
    output[i] = result;
  }

)MSL"
