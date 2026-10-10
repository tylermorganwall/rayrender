R"MSL(
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
    if (probe.material.type == bsdf_openpbr) {
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
