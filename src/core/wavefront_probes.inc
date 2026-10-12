R"MSL(

  struct SubsurfaceProbe { float4 pole,axis,random,wi,rates; };
  struct SubsurfaceResult { float4 direction,guide; };
  kernel void probe_subsurface(device const SubsurfaceProbe *inputs [[buffer(0)]],
                               device SubsurfaceResult *outputs [[buffer(1)]],
                               uint i [[thread_position_in_grid]]) {
    SubsurfaceProbe s=inputs[i];
    float3 ordinary=phase_sample(float3(0,0,1),0,s.random.y,s.random.z);
    float3 direction=subsurface_direction(s.pole.xyz,s.axis.xyz,s.random.x,
                                          s.random.y,s.random.z,ordinary);
    outputs[i].direction=float4(direction,subsurface_direction_pdf(s.pole.xyz,s.axis.xyz,direction));
    float3 rates,ordinary_weights,guided_weights;
    subsurface_flight_mixture(s.pole.xyz,s.axis.xyz,s.wi.xyz,s.rates.xyz,
                              rates,ordinary_weights,guided_weights);
    outputs[i].guide=float4(subsurface_guide_pdf(s.pole.xyz,dot(s.axis.xyz,s.wi.xyz)),
                            dot(ordinary_weights+guided_weights,float3(1)));
  }

  kernel void probe_light_selection(device const float4 *inputs [[buffer(0)]],
                                    device float4 *outputs [[buffer(1)]],
                                    constant Params &p [[buffer(2)]],
                                    device const LightNode *nodes [[buffer(3)]],
                                    uint i [[thread_position_in_grid]]) {
    uint choice=select_light(inputs[i].w,inputs[i].xyz,p,nodes);
    float sum=0;
    for (uint j=0;j<p.lights;++j) sum+=light_selection_pdf(j,inputs[i].xyz,p,nodes);
    outputs[i]=float4(choice,light_selection_pdf(choice,inputs[i].xyz,p,nodes),sum,
                      light_selection_pdf(i%p.lights,inputs[i].xyz,p,nodes));
  }

  struct EnvironmentProbe {
    uint bits;
    float latitude;
    uint padding[2];
  };
  kernel void probe_environment(device const EnvironmentProbe *inputs [[buffer(0)]],
                                device float4 *outputs [[buffer(1)]],
                                constant Light &light [[buffer(2)]],
                                device const float4 *texels [[buffer(3)]],
                                uint i [[thread_position_in_grid]]) {
    float3 direction = sample_environment(light, inputs[i].bits, inputs[i].latitude, texels);
    outputs[i] = float4(direction, light_pdf(light, direction, texels));
  }
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
                        p.normal.xyz, p.object_normal.xyz, p.footprint, int(p.uv.w)};
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
                             device const GeometryStorage *storage [[buffer(28)]],
                             instance_acceleration_structure mixed_scene [[buffer(29)]],
                             intersection_function_table<triangle_data, instancing> mixed_functions
                             [[buffer(30)]],
                             uint i [[thread_position_in_grid]]) {
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const Quadric *quadrics =
        reinterpret_cast<device const Quadric *>(bytes + storage->quadric_offset);
    device const MeshInstance *instances =
        reinterpret_cast<device const MeshInstance *>(bytes + storage->instance_offset);
    GeometryScene scene = {triangle_scene, functions, mixed_scene, mixed_functions, instances,
      reinterpret_cast<device const float4 *>(bytes + storage->motion_offset),
      reinterpret_cast<device const Quadric *>(bytes + storage->motion_offset + 48)};
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
      auto surface = geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics,
                                      surfaces, hit.instance, scene.instances, r.origin, r.direction);
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
R"MSL(
  struct CameraProbe { float4 film_lens,time; };
  struct CameraResult { float4 origin,direction; };
  kernel void probe_camera(device const CameraProbe *input [[buffer(0)]],device CameraResult *output [[buffer(1)]],
                            device const float4 *data [[buffer(2)]],constant Camera &camera [[buffer(3)]],
                            uint i [[thread_position_in_grid]]) {
    auto p=input[i];float3 o,d;
    float2 film=p.film_lens.xy;if(camera.flip_x) film.x=1-film.x;
    float weight=camera_ray(camera,film,p.film_lens.zw,p.time.x*camera.shutter.z,data,o,d);
    output[i]={float4(o,weight),float4(d,0)};
  }
  struct HazeProbe { float4 values[2],start,delta; };
  kernel void probe_haze(device const HazeProbe *input [[buffer(0)]],device CameraResult *output [[buffer(1)]],
                          uint i [[thread_position_in_grid]]) {
    auto p=input[i];float values[8];for(uint j=0;j<8;++j) values[j]=p.values[j/4][j%4];
    float edges[5]={};float4 polynomial;
    uint count=haze_edges(values,p.start.xyz,p.delta.xyz,p.delta.w,edges,polynomial);
    uint selected=0;
    for(uint j=1;j<count;++j)
      if(haze_polynomial(polynomial,(edges[j-1]+edges[j])*.5f)<0) selected|=1u<<(j-1);
    output[i]={float4(edges[0],edges[1],edges[2],edges[3]),float4(edges[4],count,selected,0)};
  }
  struct MediumProbe { Medium medium; float4 origin, direction, uniforms[4]; };
  struct MediumResult { float4 distances[4], majorants[4]; };
  kernel void probe_medium(device const MediumProbe *input [[buffer(0)]],
                           device MediumResult *output [[buffer(1)]],
                           device const float *values [[buffer(2)]],
                           uint i [[thread_position_in_grid]]) {
    auto p = input[i];
    auto iterator = medium_iterator(p.medium, p.origin.xyz, p.direction.xyz, p.origin.w);
    for (uint j = 0; j < 16; ++j) {
      float majorant;
      output[i].distances[j/4][j%4] = medium_next(iterator, p.uniforms[j/4][j%4], values, majorant);
      output[i].majorants[j/4][j%4] = majorant;
    }
  }
)MSL"
