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

// Narrow editor access keeps material types and their shading implementation fixed.
struct PreviewMaterialAccess {
  std::vector<PreviewMaterialBinding> fields;
  void color(const std::string& name, point3f value, std::function<void(point3f)> set) {
    PreviewField f;
    f.name = name;
    f.color = true;
    f.count = 3;
    f.values = Values(value);
    for (double v : f.values) {
      if (v < 0 || v > 1) {
        f.color = false;
      }
    }
    // Keep existing values outside [0,1] editable as numbers; the bounded RGB
    // picker would otherwise reject valid high-intensity material data.
    if (!f.color) {
      f.minimum = std::min(0.0, *std::min_element(f.values.begin(), f.values.end()));
      f.maximum =
          std::max(100000.0, *std::max_element(f.values.begin(), f.values.end()));
    }
    fields.push_back({f, [set](const auto& v) {
                        set(Color(v));
                      }});
  }

  void scalar(const std::string& name, double value, double lo, double hi,
              std::function<void(double)> set) {
    PreviewField f;
    f.name = name;
    f.values[0] = value;
    // Include the current value even when it lies outside the suggested UI range.
    f.minimum = std::min(lo, value);
    f.maximum = std::max(hi, value);
    fields.push_back({f, [set](const auto& v) {
                        set(v[0]);
                      }});
  }

  // Replace texture objects when editing colors so other owners of a shared
  // texture retain their values. Checker patterns and gradient settings survive.
  void texture_colors(std::shared_ptr<texture>& target) {
    if (auto p = dynamic_cast<constant_texture*>(target.get())) {
      color("Color", p->color, [&target](point3f c) {
        target = std::make_shared<constant_texture>(c);
      });
    } else if (auto p = dynamic_cast<checker_texture*>(target.get())) {
      auto even = dynamic_cast<constant_texture*>(p->even.get());
      auto odd = dynamic_cast<constant_texture*>(p->odd.get());
      if (even && odd) {
        color("Checker color", even->color, [&target](point3f c) {
          auto copy = std::make_shared<checker_texture>(
              *static_cast<checker_texture*>(target.get()));
          copy->even = std::make_shared<constant_texture>(c);
          target = copy;
        });
        color("Base color", odd->color, [&target](point3f c) {
          auto copy = std::make_shared<checker_texture>(
              *static_cast<checker_texture*>(target.get()));
          copy->odd = std::make_shared<constant_texture>(c);
          target = copy;
        });
      }
    } else if (auto p = dynamic_cast<gradient_texture*>(target.get())) {
      for (int i = 0; i < 2; ++i) {
        auto c = i ? p->gamma_color2 : p->gamma_color1;
        color(i ? "Gradient end" : "Gradient start",
              p->hsv ? HSVtoRGB(c) : c,
              [&target, i](point3f value) {
                auto copy = std::make_shared<gradient_texture>(
                    *static_cast<gradient_texture*>(target.get()));
                (i ? copy->gamma_color2 : copy->gamma_color1) =
                    copy->hsv ? RGBtoHSV(value) : value;
                target = copy;
              });
      }
    } else {
      // Reuse the original source on later edits rather than stacking multipliers.
      auto tinted = dynamic_cast<PreviewTintTexture*>(target.get());
      color("Texture tint", tinted ? tinted->tint : point3f(1), [&target](point3f c) {
        auto previous = dynamic_cast<PreviewTintTexture*>(target.get());
        target = std::make_shared<PreviewTintTexture>(
            previous ? previous->source : target, c);
      });
    }
  }

  // Map the stored distribution inputs back to R's roughness/gloss controls.
  // Texture-driven roughness stays fixed because two scalar controls cannot
  // represent its spatial variation. Reconstructing the distribution updates
  // its cached coefficients while preserving the chosen distribution class.
  void roughness(MicrofacetDistribution*& target, bool gloss = false) {
    point2f value;
    bool supported = false;
    if (auto d = dynamic_cast<BeckmannDistribution*>(target)) {
      value = d->roughness_input;
      supported = !d->has_roughness;
    }
    if (auto d = dynamic_cast<TrowbridgeReitzDistribution*>(target)) {
      value = d->roughness_input;
      supported = !d->has_roughness;
    }
    if (!supported) {
      return;
    }
    PreviewField f;
    f.name = gloss ? "Gloss X/Y" : "Roughness X/Y";
    f.count = 2;
    f.minimum = 0;
    f.maximum = 1;
    for (int i = 0; i < 2; ++i) {
      f.values[i] = gloss ? 1 - 2 * std::sqrt(value[i]) : std::sqrt(value[i]);
    }
    fields.push_back(
        {f, [&target, gloss](const auto& v) {
           const double x = gloss ? (1 - v[0]) / 2 : v[0],
                        y = gloss ? (1 - v[1]) / 2 : v[1];
           MicrofacetDistribution* next = nullptr;
           if (auto d = dynamic_cast<BeckmannDistribution*>(target)) {
             next = new BeckmannDistribution(
                 x * x, y * y, d->roughness, d->has_roughness, d->sampleVisibleArea);
           }
           if (auto d = dynamic_cast<TrowbridgeReitzDistribution*>(target)) {
             next = new TrowbridgeReitzDistribution(
                 x * x, y * y, d->roughness, d->has_roughness, d->sampleVisibleArea);
           }
           if (next) {
             delete target;
             target = next;
           }
         }});
  }

  // Expose parameters supported by the existing material class. Setters update
  // that class's native representation; the UI cannot switch material types.
  void describe(material* m) {
    if (auto p = dynamic_cast<lambertian*>(m)) {
      texture_colors(p->albedo);
    } else if (auto p = dynamic_cast<metal*>(m)) {
      texture_colors(p->albedo);
      scalar("Fuzz", p->fuzz, 0, 1, [p](double v) {
        p->fuzz = v;
      });
    } else if (auto p = dynamic_cast<dielectric*>(m)) {
      color("Tint", p->albedo, [p](point3f c) {
        p->albedo = c;
      });
      scalar("Index of refraction", p->ref_idx, 1, 5, [p](double v) {
        p->ref_idx = v;
      });
      for (int i = 0; i < 3; ++i) {
        scalar(std::string("Absorption ") + "RGB"[i],
               p->attenuation[i],
               0,
               100,
               [p, i](double v) {
                 p->attenuation[i] = v;
               });
      }
    } else if (auto p = dynamic_cast<orennayar*>(m)) {
      texture_colors(p->albedo);
      // Oren-Nayar stores derived coefficients; recover sigma for editing and
      // let its constructor recompute both coefficients when sigma changes.
      const double sigma =
          std::sqrt(std::max(0.0, .66 * (1 - p->A) / std::max(1e-8, 2.0 * p->A - 1)));
      scalar("Roughness (radians)", sigma, 0, 1.57, [p](double v) {
        orennayar next(p->albedo, v);
        p->A = next.A;
        p->B = next.B;
      });
    } else if (auto p = dynamic_cast<diffuse_light*>(m)) {
      texture_colors(p->emit);
      scalar("Light intensity", p->intensity, 0, 100000, [p](double v) {
        p->intensity = v;
      });
    } else if (auto p = dynamic_cast<spot_light*>(m)) {
      texture_colors(p->emit);
      scalar("Light intensity", p->intensity, 0, 100000, [p](double v) {
        p->intensity = v;
      });
    } else if (auto p = dynamic_cast<MicrofacetReflection*>(m)) {
      texture_colors(p->albedo);
      roughness(p->distribution);
    } else if (auto p = dynamic_cast<MicrofacetTransmission*>(m)) {
      texture_colors(p->albedo);
      roughness(p->distribution);
      scalar("Index of refraction", p->eta, 1, 5, [p](double v) {
        p->eta = v;
      });
    } else if (auto p = dynamic_cast<glossy*>(m)) {
      texture_colors(p->albedo);
      roughness(p->distribution, true);
      color("Specular reflectance", p->Rs, [p](point3f c) {
        p->Rs = c;
      });
    } else if (auto p = dynamic_cast<isotropic*>(m)) {
      texture_colors(p->albedo);
    } else if (auto p = dynamic_cast<hair*>(m)) {
      // Present absorption as transmission exp(-sigma_a). Clamp before taking
      // the inverse logarithm so selecting black never creates infinite absorption.
      color("Transmission color",
            point3f(std::exp(-p->sigma_a[0]),
                    std::exp(-p->sigma_a[1]),
                    std::exp(-p->sigma_a[2])),
            [p](point3f c) {
              p->sigma_a = point3f(-std::log(std::max(c[0], Float(1e-6))),
                                   -std::log(std::max(c[1], Float(1e-6))),
                                   -std::log(std::max(c[2], Float(1e-6))));
            });
      // Hair roughness affects cached scattering terms, so rebuild them together.
      scalar("Longitudinal roughness", p->beta_m, .01, 1, [p](double v) {
        *p = hair(p->sigma_a, p->eta, v, p->beta_n, p->alpha);
      });
      scalar("Azimuthal roughness", p->beta_n, .01, 1, [p](double v) {
        *p = hair(p->sigma_a, p->eta, p->beta_m, v, p->alpha);
      });
    }
  }
};

std::vector<material*> PreviewMaterials(hitable* root) {
  std::vector<material*> result;
  std::set<material*> seen;
  Gather(root, result, seen);
  return result;
}

std::vector<PreviewMaterialBinding> PreviewMaterialFields(material* mat) {
  PreviewMaterialAccess access;
  access.describe(mat);
  return std::move(access.fields);
}

// Rebind fields to the candidate material and reject stale layouts or invalid
// values. A later setter may fail after earlier setters ran, so callers apply
// this only to unpublished geometry that can be discarded on failure.
void PreviewApplyMaterial(material* mat, const PreviewMaterialEdit& edit) {
  if (mat->GetName() != edit.type) {
    throw std::runtime_error("The selected material changed; select the object again.");
  }

  auto bindings = PreviewMaterialFields(mat);
  if (bindings.size() != edit.fields.size()) {
    throw std::runtime_error("Material options changed; select the object again.");
  }

  for (size_t i = 0; i < bindings.size(); ++i) {
    auto& f = edit.fields[i];
    auto& expected = bindings[i].field;
    if (f.name != expected.name || f.count != expected.count) {
      throw std::runtime_error("Material option mismatch.");
    }
    for (unsigned j = 0; j < f.count; ++j) {
      if (!std::isfinite(f.values[j]) || f.values[j] < expected.minimum ||
          f.values[j] > expected.maximum) {
        throw std::runtime_error(
            "Material values must be finite and within their displayed range.");
      }
    }
    bindings[i].set(f.values);
  }
}
