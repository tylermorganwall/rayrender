#ifndef RAYRENDER_PREVIEW_TEXTURES_H
#define RAYRENDER_PREVIEW_TEXTURES_H

#include "../materials/texturecache.h"
#include "preview_validation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

// Editor texture recipes retain their source so imported/embedded pixels stay
// alive, and own any replacement files independently of the renderer's cache.
namespace PreviewTextures {
// Texture decoders report backend errors; attach the owning field so the
// inspector can keep its draft red without discarding a working material.
template <typename Load>
auto CheckedFile(const char* field, Load load) -> decltype(load()) {
  try {
    return load();
  } catch (const Rcpp::internal::InterruptedException&) {
    throw;
  } catch (const std::exception& error) {
    throw PreviewFieldError({field}, error.what());
  }
}
inline std::string FilePath(const std::string& value) {
  if (value.empty()) {
    return value;
  }
  const std::string path = R_ExpandFileName(value.c_str());
  if (!fs::is_regular_file(fs::path(path))) {
    throw std::runtime_error("Texture file does not exist: " + value);
  }
  return path;
}

class ColorTexture final : public texture {
public:
  enum Mode { Solid, Checker, Noise, Gradient, WorldGradient, Image, Imported };
  int mode = Solid;
  point3f color = point3f(1), secondary = point3f(0);
  double period = 3, noise_scale = 1, phase = 0, noise_intensity = 10;
  bool transpose = false, hsv = true;
  point3f start = point3f(0), end = point3f(0, 1, 0);
  std::array<double, 3> repeat{1, 1, 0};
  std::string path, loaded_path;
  std::shared_ptr<texture> source, image_source, evaluated;
  std::shared_ptr<TextureCache> file_owner;

  explicit ColorTexture(const std::shared_ptr<texture>& input) : source(input) {
    if (auto p = dynamic_cast<ColorTexture*>(input.get())) {
      *this = *p;
      return;
    }
    if (auto p = dynamic_cast<constant_texture*>(input.get())) {
      color = p->color;
    } else if (auto p = dynamic_cast<checker_texture*>(input.get())) {
      auto even = dynamic_cast<constant_texture*>(p->even.get());
      auto odd = dynamic_cast<constant_texture*>(p->odd.get());
      if (even && odd) {
        mode = Checker;
        color = odd->color;
        secondary = even->color;
        period = p->period;
      } else {
        mode = Imported;
      }
    } else if (auto p = dynamic_cast<noise_texture*>(input.get())) {
      mode = Noise;
      color = p->color;
      secondary = p->color2;
      noise_scale = p->scale;
      phase = p->phase * 180 / M_PI;
      noise_intensity = p->intensity;
    } else if (auto p = dynamic_cast<gradient_texture*>(input.get())) {
      mode = Gradient;
      hsv = p->hsv;
      color = hsv ? HSVtoRGB(p->gamma_color1) : p->gamma_color1;
      secondary = hsv ? HSVtoRGB(p->gamma_color2) : p->gamma_color2;
      transpose = p->aligned_v;
    } else if (auto p = dynamic_cast<world_gradient_texture*>(input.get())) {
      mode = WorldGradient;
      hsv = p->hsv;
      color = hsv ? HSVtoRGB(p->gamma_color1) : p->gamma_color1;
      secondary = hsv ? HSVtoRGB(p->gamma_color2) : p->gamma_color2;
      start = p->point1;
      end = p->point1 + p->dir;
    } else if (auto p = dynamic_cast<image_texture_float*>(input.get())) {
      mode = Image;
      repeat = {p->repeatu, p->repeatv, 0};
      image_source = input;
      path = loaded_path = p->preview_path;
    } else if (auto p = dynamic_cast<image_texture_char*>(input.get())) {
      mode = Image;
      repeat = {p->repeatu, p->repeatv, 0};
      image_source = input;
      path = loaded_path = p->preview_path;
    } else {
      mode = Imported;
    }
    evaluated = source;
  }

  // Build only after every field has been validated and copied into the recipe.
  // This permits changing both gradient endpoints in one atomic material edit.
  void Build() {
    switch (mode) {
    case Solid:
      evaluated = std::make_shared<constant_texture>(color);
      break;
    case Checker:
      if (std::abs(period) < 1e-8) {
        throw std::runtime_error("Checker period must be nonzero.");
      }
      evaluated = std::make_shared<checker_texture>(
          std::make_shared<constant_texture>(secondary),
          std::make_shared<constant_texture>(color),
          period);
      break;
    case Noise:
      evaluated = std::make_shared<noise_texture>(
          noise_scale, color, secondary, phase * M_PI / 180, noise_intensity);
      break;
    case Gradient:
      evaluated = std::make_shared<gradient_texture>(color, secondary, transpose, hsv);
      break;
    case WorldGradient:
      if ((end - start).squared_length() < 1e-12) {
        throw std::runtime_error("Gradient start and end must be different points.");
      }
      evaluated =
          std::make_shared<world_gradient_texture>(start, end, color, secondary, hsv);
      break;
    case Image: {
      if (!path.empty() && path != loaded_path) {
        auto owner = std::make_shared<TextureCache>();
        int width, height, channels;
        auto pixels = CheckedFile("Color texture file", [&] {
          return owner->LookupFloat(FilePath(path), width, height, channels, 4);
        });
        image_source = std::make_shared<image_texture_float>(pixels, width, height, 4);
        file_owner = owner;
        loaded_path = path;
      }
      if (auto p = dynamic_cast<image_texture_float*>(image_source.get())) {
        auto copy = std::make_shared<image_texture_float>(*p);
        copy->repeatu = repeat[0];
        copy->repeatv = repeat[1];
        evaluated = copy;
      } else if (auto p = dynamic_cast<image_texture_char*>(image_source.get())) {
        auto copy = std::make_shared<image_texture_char>(*p);
        copy->repeatu = repeat[0];
        copy->repeatv = repeat[1];
        evaluated = copy;
      } else {
        throw PreviewFieldError({"Color texture file"},
                                "Choose an image texture file.");
      }
      break;
    }
    case Imported:
      evaluated = source;
      break;
    default:
      throw std::runtime_error("Unknown color texture mode.");
    }
  }

  point3f value(Float u, Float v, const point3f& p) const override {
    const auto result = evaluated->value(u, v, p);
    return mode == Image || mode == Imported ? result * color : result;
  }
};

// Keep source bytes immutable. Range/flip edits remap samples rather than
// rewriting a shared texture cache, so sibling instances retain their roughness.
class RoughnessTexture final : public roughness_texture {
public:
  std::shared_ptr<roughness_texture> source;
  std::shared_ptr<TextureCache> file_owner;
  std::string path, loaded_path;

  explicit RoughnessTexture(const std::shared_ptr<roughness_texture>& input)
      : roughness_texture(nullptr, 0, 0, 0), source(input) {
    if (auto p = dynamic_cast<RoughnessTexture*>(input.get())) {
      *this = *p;
    } else if (input) {
      data = input->data;
      nx = input->nx;
      ny = input->ny;
      channels = input->channels;
      minimum = input->minimum;
      maximum = input->maximum;
      flip = input->flip;
      path = loaded_path = input->preview_path;
    }
  }

  void Build(bool enabled) {
    if (minimum > maximum) {
      throw PreviewFieldError({"Roughness map range"},
                              "Roughness map minimum must not exceed its maximum.");
    }
    if (enabled && !path.empty() && path != loaded_path) {
      auto owner = std::make_shared<TextureCache>();
      int width, height, count;
      auto pixels = CheckedFile("Roughness map file", [&] {
        return owner->LookupChar(FilePath(path), width, height, count, 3);
      });
      data = pixels;
      nx = width;
      ny = height;
      channels = 3;
      file_owner = owner;
      preview_path = loaded_path = path;
    }
    if (enabled && !data) {
      throw PreviewFieldError({"Roughness map file"}, "Choose a roughness map file.");
    }
  }
};
}
#endif
