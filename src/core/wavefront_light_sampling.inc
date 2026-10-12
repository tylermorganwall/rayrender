R"MSL(
struct LightNode {
  float4 center_radius2;
  float power;
  uint left, right, parent;
};
static_assert(sizeof(LightNode) == 32, "wavefront light BVH ABI mismatch");

float light_branch_probability(uint node, float3 position, device const LightNode *nodes) {
  LightNode a = nodes[nodes[node].left], b = nodes[nodes[node].right];
  float3 da = a.center_radius2.xyz - position, db = b.center_radius2.xyz - position;
  float a2 = max(max(dot(da, da), a.center_radius2.w), 1e-20f);
  float b2 = max(max(dot(db, db), b.center_radius2.w), 1e-20f);
  // Scale distances before division, avoiding overflow near point sources.
  float scale = min(a2, b2);
  float x = a.power * (scale / a2), y = b.power * (scale / b2);
  return x + y > 0 && isfinite(x + y) ? x / (x + y) : .5f;
}

float light_selection_pdf(uint light, float3 position, constant Params &p,
                          device const LightNode *nodes) {
  if (!(texture_features & 128) || light < p.environments)
    return 1.f / p.lights;
  float probability = 1;
  uint node = light - p.environments;
  while (nodes[node].parent != UINT_MAX) {
    uint parent = nodes[node].parent;
    float left = light_branch_probability(parent, position, nodes);
    probability *= nodes[parent].left == node ? left : 1 - left;
    node = parent;
  }
  return (.05f + .95f * (p.lights - p.environments) * probability) / p.lights;
}

uint select_light(float u, float3 position, constant Params &p, device const LightNode *nodes) {
  if (!(texture_features & 128))
    return min(uint(u * p.lights), p.lights - 1);
  // Infinite lights retain exactly their old 1/N probability, including the
  // atmosphere's separately evaluated source terms. Only the finite mass is
  // redistributed. Remap one uniform variate at each branch (no extra RNG).
  float scaled = u * p.lights;
  if (scaled < p.environments)
    return uint(scaled);
  uint count = p.lights - p.environments;
  u = (scaled - p.environments) / count;
  if (u < .05f)
    return p.environments + min(uint(u * (count / .05f)), count - 1);
  u = min((u - .05f) / .95f, 0x1.fffffep-1f);
  uint node = 2 * count - 2;
  while (nodes[node].right != UINT_MAX) {
    float left = light_branch_probability(node, position, nodes);
    bool first = u < left;
    u = first ? u / left : (u - left) / (1 - left);
    u = min(u, 0x1.fffffep-1f);
    node = first ? nodes[node].left : nodes[node].right;
  }
  return nodes[node].left;
}
)MSL"
