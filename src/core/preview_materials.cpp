#include "preview_scene.h"
#include "../materials/material.h"
#include "../materials/constant.h"
#include "../hitables/instance.h"
#include "../hitables/box.h"
#include "../hitables/trimesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/plymesh.h"
#include "../hitables/mesh3d.h"
#include "../volumes/boundary.h"
#include <set>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
// Preserve image/procedural texture detail while exposing an editable RGB multiplier.
class PreviewTintTexture final : public texture {
public:
  std::shared_ptr<texture> source;
  point3f tint;
  PreviewTintTexture(std::shared_ptr<texture> source, point3f tint)
      : source(std::move(source)), tint(tint) {
  }

  point3f value(Float u, Float v, const point3f& p) const override {
    return source->value(u, v, p) * tint;
  }
};
std::array<double, 3> Values(const point3f& p) {
  return {p[0], p[1], p[2]};
}

point3f Color(const std::array<double, 3>& p) {
  return point3f(p[0], p[1], p[2]);
}

// Walk through placement and acceleration wrappers to find the materials beneath
// a selected root. Deduplicate shared pointers while retaining first-seen slot order;
// this traversal discovers material slots without changing object selection depth.
void Gather(hitable* h, std::vector<material*>& out, std::set<material*>& seen) {
  if (!h) {
    return;
  }

  auto add = [&](material* m) {
    if (m && seen.insert(m).second) {
      out.push_back(m);
    }
  };
  if (auto p = dynamic_cast<AnimatedHitable*>(h)) {
    return Gather(p->primitive.get(), out, seen);
  }

  if (auto p = dynamic_cast<MediumBoundary*>(h)) {
    return Gather(p->geometry.get(), out, seen);
  }

  if (auto p = dynamic_cast<instance*>(h)) {
    return Gather(p->original_scene, out, seen);
  }

  if (auto p = dynamic_cast<constant_medium*>(h)) {
    add(p->phase_function.get());
    return;
  }
  // File/mesh material order is stable even when its triangle BVH is rebuilt.
  if (auto p = dynamic_cast<trimesh*>(h)) {
    for (auto& m : p->mesh->mesh_materials) {
      add(m.get());
    }
    return;
  }

  if (auto p = dynamic_cast<raymesh*>(h)) {
    for (auto& m : p->mesh->mesh_materials) {
      add(m.get());
    }
    return;
  }

  if (auto p = dynamic_cast<plymesh*>(h)) {
    for (auto& m : p->mesh->mesh_materials) {
      add(m.get());
    }
    return;
  }

  if (auto p = dynamic_cast<mesh3d*>(h)) {
    for (auto& m : p->mesh->mesh_materials) {
      add(m.get());
    }
    return;
  }

  if (auto p = dynamic_cast<BVHAggregate*>(h)) {
    for (auto& child : p->Primitives()) {
      Gather(child.get(), out, seen);
    }
    return;
  }

  if (auto p = dynamic_cast<hitable_list*>(h)) {
    for (auto& child : p->objects) {
      Gather(child.get(), out, seen);
    }
    return;
  }

  add(h->mat_ptr.get());
}
}

#include "preview_surface_maps.h"

// Material-specific bindings stay on the render thread. Only field values and
// widget metadata enter the GUI; all setters operate on an unpublished rebuild.
struct PreviewMaterialAccess {
  std::vector<PreviewMaterialBinding> fields, texture_fields;
  std::vector<std::function<void()>> finish;
  std::string section = "Surface";

  void vector(const std::string& name, std::array<double, 3> value, unsigned count,
              double lo, double hi,
              std::function<void(const std::array<double, 3>&)> set) {
    PreviewField field;
    field.name = name;
    field.section = section;
    field.values = value;
    field.count = count;
    field.minimum = lo;
    field.maximum = hi;
    for (unsigned i = 0; i < count; ++i) {
      field.minimum = std::min(field.minimum, value[i]);
      field.maximum = std::max(field.maximum, value[i]);
    }
    fields.push_back({field, std::move(set), {}});
  }

  void scalar(const std::string& name, double value, double lo, double hi,
              std::function<void(double)> set) {
    vector(name, {value, 0, 0}, 1, lo, hi, [set](const auto& values) {
      set(values[0]);
    });
  }

  void color(const std::string& name, point3f value, std::function<void(point3f)> set) {
    vector(name, Values(value), 3, 0, 1, [set](const auto& values) {
      set(Color(values));
    });
    auto& field = fields.back().field;
    field.color = field.minimum >= 0 && field.maximum <= 1;
    if (!field.color) {
      field.maximum = std::max(field.maximum, 100000.0);
    }
  }

  void boolean(const std::string& name, bool value, std::function<void(bool)> set) {
    scalar(name, value, 0, 1, [set](double value) {
      set(value != 0);
    });
    fields.back().field.boolean = true;
  }

  void choice(const std::string& name, int value, std::vector<std::string> options,
              std::function<void(int)> set) {
    scalar(name, value, 0, options.size() - 1, [set](double value) {
      set(static_cast<int>(value));
    });
    fields.back().field.choices = std::move(options);
  }

  void file(const std::string& name, const std::string& value,
            std::function<void(const std::string&)> set) {
    PreviewField field;
    field.name = name;
    field.section = section;
    field.text_input = true;
    field.text = value;
    fields.push_back({field, {}, std::move(set)});
  }

  void visible(const std::string& control, std::vector<int> choices) {
    fields.back().field.condition = control;
    fields.back().field.visible_choices = std::move(choices);
  }

  // A stable recipe schema permits switching texture modes without changing
  // saved field identities or dropping edits when a material is rebuilt.
  void texture_colors(std::shared_ptr<texture>& target) {
    auto draft = std::make_shared<PreviewTextures::ColorTexture>(target);
    color("Color", draft->color, [draft](point3f value) {
      draft->color = value;
    });
    const size_t begin = fields.size();
    section = "Color texture";
    choice("Texture mode",
           draft->mode,
           {"Solid",
            "Checker",
            "Noise",
            "UV gradient",
            "World gradient",
            "Image",
            "Imported"},
           [draft](int value) {
             draft->mode = value;
           });
    color("Secondary color", draft->secondary, [draft](point3f value) {
      draft->secondary = value;
    });
    visible("Texture mode", {1, 2, 3, 4});
    scalar("Checker period", draft->period, .0001, 1000, [draft](double value) {
      draft->period = value;
    });
    visible("Texture mode", {1});
    scalar("Noise scale", draft->noise_scale, 0, 1000, [draft](double value) {
      draft->noise_scale = value;
    });
    visible("Texture mode", {2});
    scalar("Noise phase (degrees)", draft->phase, -360, 360, [draft](double value) {
      draft->phase = value;
    });
    visible("Texture mode", {2});
    scalar("Noise intensity", draft->noise_intensity, 0, 1000, [draft](double value) {
      draft->noise_intensity = value;
    });
    visible("Texture mode", {2});
    boolean("Gradient transpose", draft->transpose, [draft](bool value) {
      draft->transpose = value;
    });
    visible("Texture mode", {3});
    choice(
        "Gradient color space", draft->hsv ? 0 : 1, {"HSV", "RGB"}, [draft](int value) {
          draft->hsv = value == 0;
        });
    visible("Texture mode", {3, 4});
    vector("Gradient start XYZ",
           Values(draft->start),
           3,
           -100000,
           100000,
           [draft](const auto& value) {
             draft->start = Color(value);
           });
    visible("Texture mode", {4});
    vector("Gradient end XYZ",
           Values(draft->end),
           3,
           -100000,
           100000,
           [draft](const auto& value) {
             draft->end = Color(value);
           });
    visible("Texture mode", {4});
    file("Color texture file", draft->path, [draft](const std::string& value) {
      draft->path = value;
    });
    visible("Texture mode", {5});
    vector(
        "Image repeat U/V", draft->repeat, 2, .0001, 1000, [draft](const auto& value) {
          draft->repeat = value;
        });
    visible("Texture mode", {5});
    texture_fields.insert(texture_fields.end(), fields.begin() + begin, fields.end());
    fields.erase(fields.begin() + begin, fields.end());
    section = "Surface";
    finish.push_back([draft, &target] {
      draft->Build();
      target = draft;
    });
  }

  // Reconstruct the distribution once after all controls have been applied.
  // This refreshes cached coefficients and keeps its map storage instance-local.
  void roughness(MicrofacetDistribution*& target, bool gloss = false) {
    struct Distribution {
      int kind = 0;
      point2f input;
      bool visible = true, mapped = false;
      std::shared_ptr<roughness_texture> source;
    };
    auto draft = std::make_shared<Distribution>();
    if (auto p = dynamic_cast<TrowbridgeReitzDistribution*>(target)) {
      draft->input = p->roughness_input;
      draft->mapped = p->has_roughness;
      draft->source = p->roughness;
      draft->visible = p->sampleVisibleArea;
    } else if (auto p = dynamic_cast<BeckmannDistribution*>(target)) {
      draft->kind = 1;
      draft->input = p->roughness_input;
      draft->mapped = p->has_roughness;
      draft->source = p->roughness;
      draft->visible = p->sampleVisibleArea;
    } else {
      return;
    }
    auto map = std::make_shared<PreviewTextures::RoughnessTexture>(draft->source);
    choice("Microfacet distribution",
           draft->kind,
           {"Trowbridge-Reitz", "Beckmann"},
           [draft](int value) {
             draft->kind = value;
           });
    std::array<double, 3> value{};
    for (int i = 0; i < 2; ++i) {
      value[i] =
          gloss ? 1 - 2 * std::sqrt(draft->input[i]) : std::sqrt(draft->input[i]);
    }
    vector(gloss ? "Gloss X/Y" : "Roughness X/Y",
           value,
           2,
           0,
           1,
           [draft, gloss](const auto& value) {
             for (int i = 0; i < 2; ++i) {
               const double roughness = gloss ? (1 - value[i]) / 2 : value[i];
               draft->input.e[i] = roughness * roughness;
             }
           });
    const size_t begin = fields.size();
    section = "Roughness map";
    boolean("Use roughness map", draft->mapped, [draft](bool value) {
      draft->mapped = value;
    });
    file("Roughness map file", map->path, [map](const std::string& value) {
      map->path = value;
    });
    visible("Use roughness map", {1});
    vector("Roughness map range",
           {map->minimum, map->maximum, 0},
           2,
           0,
           1,
           [map](const auto& value) {
             map->minimum = value[0];
             map->maximum = value[1];
           });
    visible("Use roughness map", {1});
    boolean("Flip roughness map", map->flip, [map](bool value) {
      map->flip = value;
    });
    visible("Use roughness map", {1});
    texture_fields.insert(texture_fields.end(), fields.begin() + begin, fields.end());
    fields.erase(fields.begin() + begin, fields.end());
    section = "Surface";
    finish.push_back([draft, map, &target] {
      map->Build(draft->mapped);
      std::unique_ptr<MicrofacetDistribution> next;
      if (draft->kind == 0) {
        next = std::make_unique<TrowbridgeReitzDistribution>(
            draft->input[0], draft->input[1], map, draft->mapped, draft->visible);
      } else {
        next = std::make_unique<BeckmannDistribution>(
            draft->input[0], draft->input[1], map, draft->mapped, draft->visible);
      }
      delete target;
      target = next.release();
    });
  }

  void surface_maps(material* mat, hitable* root) {
    auto targets = std::make_shared<PreviewTextures::SurfaceTargets>();
    targets->Gather(root, mat);
    if (!targets->densities.empty()) {
      section = "Volume";
      scalar("Fog density",
             *targets->densities.front(),
             .000001,
             1000,
             [targets](double value) {
               for (auto density : targets->densities) {
                 *density = value;
               }
             });
    }
    if (targets->alpha.empty() && targets->bump.empty()) {
      return;
    }
    auto draft = std::make_shared<PreviewSurfaceMapSettings>();
    if (mat->preview_maps) {
      *draft = *mat->preview_maps;
    } else {
      if (!targets->alpha.empty()) {
        draft->alpha = *targets->alpha.front();
        draft->alpha_enabled = bool(draft->alpha);
        draft->alpha_path = draft->alpha ? draft->alpha->preview_path : "";
      }
      if (!targets->bump.empty()) {
        draft->bump = *targets->bump.front();
        draft->bump_enabled = bool(draft->bump);
        if (draft->bump) {
          draft->bump_path = draft->bump->preview_path;
          draft->bump_intensity = draft->bump->intensity;
          draft->bump_repeat = {draft->bump->repeatu, draft->bump->repeatv, 0};
        }
      }
    }
    if (!targets->alpha.empty()) {
      section = "Alpha map";
      boolean("Use alpha map", draft->alpha_enabled, [draft](bool value) {
        draft->alpha_enabled = value;
      });
      file("Alpha map file", draft->alpha_path, [draft](const std::string& value) {
        draft->alpha_path = value;
      });
      visible("Use alpha map", {1});
    }
    if (!targets->bump.empty()) {
      section = "Bump map";
      boolean("Use bump map", draft->bump_enabled, [draft](bool value) {
        draft->bump_enabled = value;
      });
      file("Bump map file", draft->bump_path, [draft](const std::string& value) {
        draft->bump_path = value;
      });
      visible("Use bump map", {1});
      scalar("Bump intensity", draft->bump_intensity, -100, 100, [draft](double value) {
        draft->bump_intensity = value;
      });
      visible("Use bump map", {1});
      vector("Bump repeat U/V",
             draft->bump_repeat,
             2,
             .0001,
             1000,
             [draft](const auto& value) {
               draft->bump_repeat = value;
             });
      visible("Use bump map", {1});
    }
    finish.push_back([draft, targets, mat] {
      draft->Build();
      for (auto target : targets->alpha) {
        *target = draft->alpha_enabled ? draft->alpha : nullptr;
      }
      for (auto target : targets->bump) {
        *target = draft->bump_enabled ? draft->bump : nullptr;
      }
      mat->preview_maps = draft;
      // Refresh inner mesh/instance acceleration metadata before the enclosing
      // scene builds its own shadow classification from these objects.
      for (auto aggregate : targets->aggregates) {
        if (aggregate) {
          aggregate->RefreshShadowType();
        }
      }
    });
  }

  void optical(const std::string& name, point3f value,
               std::function<void(point3f)> set) {
    vector(name, Values(value), 3, 0, 100, [set](const auto& value) {
      set(Color(value));
    });
  }

  void describe(material* mat, hitable* root = nullptr) {
    if (auto p = dynamic_cast<lambertian*>(mat)) {
      texture_colors(p->albedo);
      scalar("Sigma (degrees)", p->preview_sigma, 0, 180, [p](double value) {
        p->preview_sigma = value;
      });
      finish.push_back([p] {
        p->preview_rough_model =
            p->preview_sigma > 0
                ? std::make_shared<orennayar>(p->albedo, p->preview_sigma * M_PI / 180)
                : nullptr;
      });
    } else if (auto p = dynamic_cast<metal*>(mat)) {
      texture_colors(p->albedo);
      scalar("Fuzz", p->fuzz, 0, 1, [p](double value) {
        p->fuzz = value;
      });
      optical("Eta RGB", p->eta, [p](point3f value) {
        p->eta = value;
      });
      optical("Kappa RGB", p->k, [p](point3f value) {
        p->k = value;
      });
    } else if (auto p = dynamic_cast<dielectric*>(mat)) {
      color("Tint", p->albedo, [p](point3f value) {
        p->albedo = value;
      });
      scalar("Index of refraction", p->ref_idx, .001, 5, [p](double value) {
        p->ref_idx = value;
      });
      for (int i = 0; i < 3; ++i) {
        scalar(std::string("Absorption ") + "RGB"[i],
               p->attenuation[i],
               0,
               100,
               [p, i](double value) {
                 p->attenuation[i] = value;
               });
      }
      scalar("Priority", p->priority, 0, 100000, [p](double value) {
        p->priority = static_cast<size_t>(value);
      });
      fields.back().field.integer = true;
      fields.back().field.speed = 1;
    } else if (auto p = dynamic_cast<orennayar*>(mat)) {
      texture_colors(p->albedo);
      const double sigma =
          std::sqrt(std::max(0.0, .66 * (1 - p->A) / std::max(1e-8, 2.0 * p->A - 1)));
      scalar("Sigma (degrees)", sigma * 180 / M_PI, 0, 180, [p](double value) {
        orennayar next(p->albedo, value * M_PI / 180);
        p->A = next.A;
        p->B = next.B;
      });
    } else if (auto p = dynamic_cast<diffuse_light*>(mat)) {
      texture_colors(p->emit);
      scalar("Light intensity", p->intensity, 0, 100000, [p](double value) {
        p->intensity = value;
      });
      boolean("Invisible", p->invisible, [p](bool value) {
        p->invisible = value;
      });
    } else if (auto p = dynamic_cast<spot_light*>(mat)) {
      texture_colors(p->emit);
      scalar("Light intensity", p->intensity, 0, 100000, [p](double value) {
        p->intensity = value;
      });
      boolean("Invisible", p->invisible, [p](bool value) {
        p->invisible = value;
      });
      vector("Spotlight direction XYZ",
             {p->spot_direction[0], p->spot_direction[1], p->spot_direction[2]},
             3,
             -100000,
             100000,
             [p](const auto& value) {
               vec3f direction(value[0], value[1], value[2]);
               if (direction.squared_length() < 1e-12) {
                 throw std::runtime_error("Spotlight direction must be nonzero.");
               }
               p->spot_direction = unit_vector(direction);
             });
      scalar("Spotlight width (degrees)",
             std::acos(std::clamp(double(p->cosTotalWidth), -1.0, 1.0)) * 180 / M_PI,
             0,
             180,
             [p](double value) {
               p->cosTotalWidth = std::cos(value * M_PI / 180);
             });
      scalar("Falloff start (degrees)",
             std::acos(std::clamp(double(p->cosFalloffStart), -1.0, 1.0)) * 180 / M_PI,
             0,
             180,
             [p](double value) {
               p->cosFalloffStart = std::cos(value * M_PI / 180);
             });
      finish.push_back([p] {
        if (p->cosFalloffStart < p->cosTotalWidth) {
          throw std::runtime_error("Spotlight falloff must start within its width.");
        }
      });
    } else if (auto p = dynamic_cast<MicrofacetReflection*>(mat)) {
      texture_colors(p->albedo);
      roughness(p->distribution);
      optical("Eta RGB", p->eta, [p](point3f value) {
        p->eta = value;
      });
      optical("Kappa RGB", p->k, [p](point3f value) {
        p->k = value;
      });
    } else if (auto p = dynamic_cast<MicrofacetTransmission*>(mat)) {
      texture_colors(p->albedo);
      roughness(p->distribution);
      scalar("Index of refraction", p->eta, .001, 5, [p](double value) {
        p->eta = value;
      });
    } else if (auto p = dynamic_cast<glossy*>(mat)) {
      texture_colors(p->albedo);
      roughness(p->distribution, true);
      vector("Specular reflectance", Values(p->Rs), 3, 0, 1, [p](const auto& value) {
        p->Rs = Color(value);
      });
    } else if (auto p = dynamic_cast<isotropic*>(mat)) {
      texture_colors(p->albedo);
    } else if (auto p = dynamic_cast<hair*>(mat)) {
      // Match hair()'s pigment and reflectance conversions; absorption is the
      // canonical mode for imported materials whose original inputs are unknown.
      auto reflectance_scale = [](double beta) {
        return 5.969 - .215 * beta + 2.532 * std::pow(beta, 2) -
               10.73 * std::pow(beta, 3) + 5.574 * std::pow(beta, 4) +
               .245 * std::pow(beta, 5);
      };
      auto color_value = p->preview_color;
      if (p->preview_color_mode == 0) {
        for (int i = 0; i < 3; ++i) {
          color_value[i] = std::exp(-std::sqrt(std::max(Float(0), p->sigma_a[i])) *
                                    reflectance_scale(p->beta_n));
        }
      }
      choice("Hair color mode",
             p->preview_color_mode,
             {"Absorption", "Color", "Pigment"},
             [p](int value) {
               p->preview_color_mode = value;
             });
      optical("Absorption RGB", p->sigma_a, [p](point3f value) {
        p->sigma_a = value;
      });
      visible("Hair color mode", {0});
      color("Hair color", color_value, [p](point3f value) {
        p->preview_color = value;
      });
      visible("Hair color mode", {1});
      scalar("Pigment", p->preview_pigment, 0, 20, [p](double value) {
        p->preview_pigment = value;
      });
      visible("Hair color mode", {2});
      scalar("Red pigment", p->preview_red_pigment, 0, 20, [p](double value) {
        p->preview_red_pigment = value;
      });
      visible("Hair color mode", {2});
      scalar("Index of refraction", p->eta, .001, 5, [p](double value) {
        p->eta = value;
      });
      scalar("Longitudinal roughness", p->beta_m, .01, 1, [p](double value) {
        p->beta_m = value;
      });
      scalar("Azimuthal roughness", p->beta_n, .01, 1, [p](double value) {
        p->beta_n = value;
      });
      scalar("Scale angle (degrees)", p->alpha, -90, 90, [p](double value) {
        p->alpha = value;
      });
      finish.push_back([p, reflectance_scale] {
        if (p->preview_color_mode == 1) {
          for (int i = 0; i < 3; ++i) {
            const double ratio = std::log(std::max(Float(1e-6), p->preview_color[i])) /
                                 reflectance_scale(p->beta_n);
            p->sigma_a[i] = ratio * ratio;
          }
        } else if (p->preview_color_mode == 2) {
          p->sigma_a = p->preview_pigment * point3f(.419, .697, 1.37) +
                       p->preview_red_pigment * point3f(.187, .4, 1.05);
        }
        hair next(p->sigma_a, p->eta, p->beta_m, p->beta_n, p->alpha);
        next.preview_color_mode = p->preview_color_mode;
        next.preview_color = p->preview_color;
        next.preview_pigment = p->preview_pigment;
        next.preview_red_pigment = p->preview_red_pigment;
        *p = next;
      });
    }
    fields.insert(fields.end(), texture_fields.begin(), texture_fields.end());
    surface_maps(mat, root);
  }
};

std::vector<material*> PreviewMaterials(hitable* root) {
  std::vector<material*> result;
  std::set<material*> seen;
  Gather(root, result, seen);
  return result;
}

std::vector<PreviewMaterialBinding> PreviewMaterialFields(material* mat,
                                                          hitable* root) {
  PreviewMaterialAccess access;
  access.describe(mat, root);
  return std::move(access.fields);
}

// Reject stale layouts and invalid values before running any setters. Setters and
// finalizers touch only the unpublished candidate; failures discard that rebuild.
void PreviewApplyMaterial(material* mat, const PreviewMaterialEdit& edit,
                          hitable* root) {
  if (mat->GetName() != edit.type) {
    throw std::runtime_error("The selected material changed; select the object again.");
  }
  PreviewMaterialAccess access;
  access.describe(mat, root);
  auto& bindings = access.fields;
  if (bindings.size() != edit.fields.size()) {
    throw std::runtime_error("Material options changed; select the object again.");
  }
  for (size_t i = 0; i < bindings.size(); ++i) {
    const auto& field = edit.fields[i];
    const auto& expected = bindings[i].field;
    if (field.name != expected.name || field.count != expected.count ||
        field.text_input != expected.text_input) {
      throw std::runtime_error("Material option mismatch.");
    }
    if (expected.text_input) {
      if (field.text.size() >= 4096 || field.text.find('\0') != std::string::npos) {
        throw std::runtime_error("Texture paths must be shorter than 4096 bytes.");
      }
      continue;
    }
    for (unsigned j = 0; j < field.count; ++j) {
      const double value = field.values[j];
      if (!std::isfinite(value) || value < expected.minimum ||
          value > expected.maximum ||
          ((expected.integer || expected.boolean || !expected.choices.empty()) &&
           value != std::floor(value))) {
        throw std::runtime_error(
            "Material values must be finite and within their displayed range.");
      }
    }
  }
  for (size_t i = 0; i < bindings.size(); ++i) {
    if (bindings[i].field.text_input) {
      bindings[i].set_text(edit.fields[i].text);
    } else {
      bindings[i].set(edit.fields[i].values);
    }
  }
  for (auto& finalize : access.finish) {
    finalize();
  }
}
