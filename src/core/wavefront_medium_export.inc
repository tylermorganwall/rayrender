// Called once per placement; native lattice arrays are shared between copies.
std::unordered_map<const SampledField *, WFField> volume_fields;
WFField blackbody_field;
std::unordered_map<const Medium *, WFMedium> sparse_media;
std::unordered_map<const Medium *, WFField> majorant_fields;
WFField ExportField(const SampledField &field, int channel = 0) {
  if (field.values.empty())
    return {};
  auto found = volume_fields.find(&field);
  WFField out;
  const size_t plane = size_t(field.dims[0]) * field.dims[1] * field.dims[2];
  if (found == volume_fields.end()) {
    if (scene.fields.size() + field.values.size() > UINT32_MAX)
      throw std::runtime_error("Metal medium field capacity exceeded");
    out = {uint32_t(scene.fields.size()), uint32_t(field.dims[0]), uint32_t(field.dims[1]), uint32_t(field.dims[2])};
    scene.fields.insert(scene.fields.end(), field.values.begin(), field.values.end());
    volume_fields.emplace(&field, out);
  } else
    out = found->second;
  out.offset += channel * plane;
  return out;
}

uint32_t ExportMedium(const Medium &medium, const Transform &placement) {
  scene.texture_features |= 16; // Specialize general-medium work out of ordinary/SSS shaders.
  WFMedium out;
  const auto inverse = placement.GetInverseMatrix();
  for (int k = 0; k < 3; ++k)
    out.inverse[k] = {float(inverse.m[k][0]), float(inverse.m[k][1]), float(inverse.m[k][2]), float(inverse.m[k][3])};
  out.emission = Pack(medium.Emission(point3f(0)));
  out.emission.w = medium.haze;
  out.controls = {float(medium.emission_scale), float(medium.temperature_scale), float(medium.temperature_offset), 0};
  if (const auto *grid = dynamic_cast<const GridMedium *>(&medium)) {
    out.lo = Pack(grid->majorants.lo);
    out.hi = Pack(grid->majorants.hi);
    out.hi.w = *std::max_element(grid->majorants.density.begin(), grid->majorants.density.end());
    out.density = ExportField(grid->density);
    out.temperature = ExportField(grid->temperatures);
    out.red = ExportField(grid->emissions, 0);
    out.green = ExportField(grid->emissions, 1);
    out.blue = ExportField(grid->emissions, 2);
    out.scale = ExportField(grid->emission_scales);
    out.controls.w = 1;
    auto cached = majorant_fields.find(&medium);
    if (cached == majorant_fields.end()) {
      uint32_t n = grid->majorants.resolution;
      out.majorant = {uint32_t(scene.fields.size()), n, n, n};
      scene.fields.insert(scene.fields.end(), grid->majorants.density.begin(), grid->majorants.density.end());
      majorant_fields.emplace(&medium, out.majorant);
    } else
      out.majorant = cached->second;
    if (!grid->temperatures.values.empty()) {
      for (float value : grid->temperatures.values)
        if ((value - medium.temperature_offset) * medium.temperature_scale > 20000)
          throw std::runtime_error("Metal grid temperature exceeds 20000 kelvin");
    }
  } else if (const auto *sparse = dynamic_cast<const NanoVDBMedium *>(&medium)) {
    auto found = sparse_media.find(&medium);
    if (found == sparse_media.end()) {
      sparse->ExportWavefront(out, scene.fields);
      sparse_media.emplace(&medium, out);
    } else {
      auto inverse = out;
      out = found->second;
      for (int k = 0; k < 3; ++k)
        out.inverse[k] = inverse.inverse[k];
    }
  } else if (!medium.IsHomogeneous())
    throw std::runtime_error("unrecognized medium representation");
  if (out.temperature.nx) {
    if (!blackbody_field.nx) {
      blackbody_field = {uint32_t(scene.fields.size()), 1991, 1, 1};
      for (int i = 0; i < 1991; ++i) {
        auto rgb = BlackbodyRGB(100 + 10 * i);
        for (int c = 0; c < 3; ++c)
          scene.fields.push_back(rgb[c]);
      }
    }
    out.blackbody = blackbody_field;
  }
  // Density samples are unscaled; extinction already contains density_scale.
  // Keep the cutoff in those same unscaled units, including zero-density media.
  out.lo.w = medium.haze_density_threshold > 0 && medium.density_scale > 0
                 ? float(std::max(double(std::numeric_limits<float>::min()),
                                  medium.haze_density_threshold / medium.density_scale)) : 0;
  uint32_t id = scene.media.size();
  scene.media.push_back(out);
  scene.native_media.push_back(&medium);
  return id;
}
