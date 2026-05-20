#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION 
#endif

#include "math/float.h"
#include "math/vectypes.h"
#include "math/mathinline.h"
#include "math/transform.h"
#include "math/transformcache.h"
#include "core/camera.h"
#include "math/float.h"
#include "core/buildscene.h"
#include "rng.h"
#include "hitables/infinite_area_light.h"
#include "core/adaptivesampler.h"
#include "math/sampler.h"
#include "core/color.h"
#include "core/integrator.h"
#include "utils/debug.h"

#include "materials/texturecache.h"
#include "hitables/box.h"
#include "hitables/sphere.h"
#include "core/PreviewDisplay.h"
#include "utils/raylog.h"
#include <cfenv>

// #define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

#ifdef HAS_OIDN
#undef None
#include <OpenImageDenoise/oidn.hpp>
#endif


#include "RProgress.h"
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "math/RayMatrix.h"
#include "RcppThread.h"

using namespace std;

static std::vector<PreviewTextOverlay> parse_preview_text_overlays(List render_info) {
  std::vector<PreviewTextOverlay> overlays;
  if(!render_info.containsElementNamed("screen_text_preview")) {
    return overlays;
  }
  List screen_text_preview = as<List>(render_info["screen_text_preview"]);
  if(!screen_text_preview.containsElementNamed("active") ||
     !as<bool>(screen_text_preview["active"])) {
    return overlays;
  }
  if(!screen_text_preview.containsElementNamed("overlays")) {
    return overlays;
  }
  List overlay_list = as<List>(screen_text_preview["overlays"]);
  overlays.reserve(overlay_list.size());
  for(int i = 0; i < overlay_list.size(); i++) {
    List overlay_r = as<List>(overlay_list[i]);
    NumericVector image = as<NumericVector>(overlay_r["image"]);
    IntegerVector dims = image.attr("dim");
    if(dims.size() < 3 || dims[2] < 4) {
      continue;
    }
    unsigned int image_height = static_cast<unsigned int>(dims[0]);
    unsigned int image_width = static_cast<unsigned int>(dims[1]);
    PreviewTextOverlay overlay;
    overlay.anchor = point3f(as<Float>(overlay_r["x"]),
                             as<Float>(overlay_r["y"]),
                             as<Float>(overlay_r["z"]));
    overlay.x_offset = static_cast<int>(std::round(as<Float>(overlay_r["x_offset"])));
    overlay.y_offset = static_cast<int>(std::round(as<Float>(overlay_r["y_offset"])));
    overlay.hjust = as<Float>(overlay_r["hjust"]);
    overlay.vjust = as<Float>(overlay_r["vjust"]);
    overlay.clip = as<bool>(overlay_r["clip"]);
    overlay.occlusion = as<bool>(overlay_r["occlusion"]);
    overlay.partial_occlusion = overlay_r.containsElementNamed("partial_occlusion") ?
      as<bool>(overlay_r["partial_occlusion"]) : false;
    overlay.occlusion_tolerance = as<Float>(overlay_r["occlusion_tolerance"]);
    overlay.width = image_width;
    overlay.height = image_height;
    overlay.rgba.resize(static_cast<size_t>(image_width) *
                        static_cast<size_t>(image_height) * 4);
    for(unsigned int y = 0; y < image_height; y++) {
      for(unsigned int x = 0; x < image_width; x++) {
        for(unsigned int channel = 0; channel < 4; channel++) {
          size_t r_idx = y +
            static_cast<size_t>(image_height) * x +
            static_cast<size_t>(image_height) *
              static_cast<size_t>(image_width) * channel;
          size_t c_idx = 4 * (x + static_cast<size_t>(image_width) * y) + channel;
          Float value = clamp(static_cast<Float>(image[r_idx]), 0.f, 1.f);
          overlay.rgba[c_idx] = static_cast<unsigned char>(std::round(255.f * value));
        }
      }
    }
    overlays.push_back(std::move(overlay));
  }
  return overlays;
}

static std::vector<PreviewLineOverlay> parse_preview_line_overlays(List render_info) {
  std::vector<PreviewLineOverlay> overlays;
  if(!render_info.containsElementNamed("screen_line_preview")) {
    return overlays;
  }
  List screen_line_preview = as<List>(render_info["screen_line_preview"]);
  if(!screen_line_preview.containsElementNamed("active") ||
     !as<bool>(screen_line_preview["active"])) {
    return overlays;
  }
  if(!screen_line_preview.containsElementNamed("lines")) {
    return overlays;
  }
  List line_list = as<List>(screen_line_preview["lines"]);
  overlays.reserve(line_list.size());
  for(int i = 0; i < line_list.size(); i++) {
    List line_r = as<List>(line_list[i]);
    PreviewLineOverlay overlay;
    overlay.start = point3f(as<Float>(line_r["x"]),
                            as<Float>(line_r["y"]),
                            as<Float>(line_r["z"]));
    overlay.end = point3f(as<Float>(line_r["xend"]),
                          as<Float>(line_r["yend"]),
                          as<Float>(line_r["zend"]));
    overlay.x_offset = static_cast<int>(std::round(as<Float>(line_r["x_offset"])));
    overlay.y_offset = static_cast<int>(std::round(as<Float>(line_r["y_offset"])));
    overlay.xend_offset = static_cast<int>(std::round(as<Float>(line_r["xend_offset"])));
    overlay.yend_offset = static_cast<int>(std::round(as<Float>(line_r["yend_offset"])));
    overlay.width = as<Float>(line_r["width"]);
    overlay.red = as<Float>(line_r["red"]);
    overlay.green = as<Float>(line_r["green"]);
    overlay.blue = as<Float>(line_r["blue"]);
    overlay.alpha = as<Float>(line_r["alpha"]);
    std::string lineend = as<std::string>(line_r["lineend"]);
    overlay.lineend = lineend == "butt" ? 1 : lineend == "square" ? 2 : 0;
    overlay.clip = as<bool>(line_r["clip"]);
    overlay.occlusion = as<bool>(line_r["occlusion"]);
    overlay.partial_occlusion = line_r.containsElementNamed("partial_occlusion") ?
      as<bool>(line_r["partial_occlusion"]) : false;
    overlay.occlusion_tolerance = as<Float>(line_r["occlusion_tolerance"]);
    overlays.push_back(overlay);
  }
  return overlays;
}

static bool is_text_anchor_occluded(const point3f& anchor,
                                    Float occlusion_tolerance,
                                    RayCamera* cam,
                                    hitable* world,
                                    random_gen& rng) {
  if(!cam || !world) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = anchor - origin;
  Float distance_to_anchor = relative.length();
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    distance_to_anchor * occlusion_tolerance : occlusion_tolerance;
  if(distance_to_anchor <= endpoint_tolerance) {
    return false;
  }
  vec3f direction = relative / distance_to_anchor;
  Float t_max = distance_to_anchor - endpoint_tolerance;

  if(cam->get_fov() == 0) {
    vec3f right = cam->get_u();
    vec3f up = cam->get_v();
    vec3f forward = cam->get_w();
    Float x_camera = dot(relative, right);
    Float y_camera = dot(relative, up);
    Float z_camera = dot(relative, forward);
    endpoint_tolerance = occlusion_tolerance < 1.f ?
      z_camera * occlusion_tolerance : occlusion_tolerance;
    if(z_camera <= endpoint_tolerance) {
      return false;
    }
    origin = origin + x_camera * right + y_camera * up;
    direction = forward;
    t_max = z_camera - endpoint_tolerance;
  }

  hit_record hrec;
  Ray visibility_ray(origin, direction, 0.5f);
  return world->hit(visibility_ray, 0.001f, t_max, hrec, rng) &&
    hrec.shape->GetName() != "EnvironmentLight";
}

static LogicalVector compute_screen_text_visibility(List render_info,
                                                    RayCamera* cam,
                                                    hitable* world,
                                                    random_gen& rng) {
  if(!render_info.containsElementNamed("screen_text_occlusion")) {
    return LogicalVector();
  }
  List screen_text_occlusion = as<List>(render_info["screen_text_occlusion"]);
  if(!screen_text_occlusion.containsElementNamed("active")) {
    return LogicalVector();
  }
  if(!screen_text_occlusion.containsElementNamed("occlusion")) {
    return LogicalVector();
  }
  LogicalVector occlusion = as<LogicalVector>(screen_text_occlusion["occlusion"]);
  int n = occlusion.size();
  LogicalVector visible(n, true);
  if(!as<bool>(screen_text_occlusion["active"])) {
    return visible;
  }
  NumericVector x = as<NumericVector>(screen_text_occlusion["x"]);
  NumericVector y = as<NumericVector>(screen_text_occlusion["y"]);
  NumericVector z = as<NumericVector>(screen_text_occlusion["z"]);
  NumericVector occlusion_tolerance =
    as<NumericVector>(screen_text_occlusion["occlusion_tolerance"]);

  for(int i = 0; i < n; i++) {
    if(!occlusion[i]) {
      continue;
    }
    visible[i] = !is_text_anchor_occluded(point3f(x[i], y[i], z[i]),
                                          occlusion_tolerance[i],
                                          cam,
                                          world,
                                          rng);
  }
  return visible;
}

static LogicalVector compute_screen_line_visibility(List render_info,
                                                    RayCamera* cam,
                                                    hitable* world,
                                                    random_gen& rng) {
  if(!render_info.containsElementNamed("screen_line_occlusion")) {
    return LogicalVector();
  }
  List screen_line_occlusion = as<List>(render_info["screen_line_occlusion"]);
  if(!screen_line_occlusion.containsElementNamed("active")) {
    return LogicalVector();
  }
  if(!screen_line_occlusion.containsElementNamed("occlusion")) {
    return LogicalVector();
  }
  LogicalVector occlusion = as<LogicalVector>(screen_line_occlusion["occlusion"]);
  int n = occlusion.size();
  LogicalVector visible(n, true);
  if(!as<bool>(screen_line_occlusion["active"])) {
    return visible;
  }
  NumericVector x = as<NumericVector>(screen_line_occlusion["x"]);
  NumericVector y = as<NumericVector>(screen_line_occlusion["y"]);
  NumericVector z = as<NumericVector>(screen_line_occlusion["z"]);
  NumericVector occlusion_tolerance =
    as<NumericVector>(screen_line_occlusion["occlusion_tolerance"]);

  for(int i = 0; i < n; i++) {
    if(!occlusion[i]) {
      continue;
    }
    visible[i] = !is_text_anchor_occluded(point3f(x[i], y[i], z[i]),
                                          occlusion_tolerance[i],
                                          cam,
                                          world,
                                          rng);
  }
  return visible;
}

static bool should_return_screen_text_overlay(List render_info) {
  if(!render_info.containsElementNamed("screen_text_native_overlay")) {
    return false;
  }
  return as<bool>(render_info["screen_text_native_overlay"]);
}

static bool should_return_screen_line_overlay(List render_info) {
  if(!render_info.containsElementNamed("screen_line_native_overlay")) {
    return false;
  }
  return as<bool>(render_info["screen_line_native_overlay"]);
}

static List get_screen_camera_info(RayCamera* cam) {
  if(!cam) {
    return List::create();
  }
  point3f origin = cam->get_origin();
  vec3f right = cam->get_u();
  vec3f up = cam->get_v();
  vec3f forward = cam->get_w();
  point2f ortho = cam->get_ortho();
  return List::create(
    _["origin"] = NumericVector::create(origin.xyz.x, origin.xyz.y, origin.xyz.z),
    _["u"] = NumericVector::create(right.xyz.x, right.xyz.y, right.xyz.z),
    _["v"] = NumericVector::create(up.xyz.x, up.xyz.y, up.xyz.z),
    _["w"] = NumericVector::create(forward.xyz.x, forward.xyz.y, forward.xyz.z),
    _["fov"] = cam->get_fov(),
    _["ortho_dimensions"] = NumericVector::create(ortho.xy.x, ortho.xy.y)
  );
}

static bool project_text_anchor_to_screen(const PreviewTextOverlay& overlay,
                                          RayCamera* cam,
                                          unsigned int display_width,
                                          unsigned int display_height,
                                          Float& screen_x,
                                          Float& screen_y) {
  if(!cam || display_width == 0 || display_height == 0) {
    return false;
  }
  Float fov = cam->get_fov();
  if(fov < 0) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = overlay.anchor - origin;
  Float rel_len2 = relative.squared_length();
  if(rel_len2 <= 0) {
    return false;
  }
  vec3f right = cam->get_u();
  vec3f up = cam->get_v();
  vec3f forward = cam->get_w();
  Float x_camera = dot(relative, right);
  Float y_camera = dot(relative, up);
  Float z_camera = dot(relative, forward);
  Float s = 0.5f;
  Float t = 0.5f;
  bool in_front = true;
  const Float pi_val = static_cast<Float>(3.14159265358979323846);

  if(fov == 0) {
    point2f ortho = cam->get_ortho();
    s = 0.5f + x_camera / ortho.xy.x;
    t = 0.5f + y_camera / ortho.xy.y;
    in_front = z_camera >= 0;
  } else if(fov == 360) {
    vec3f direction = relative / std::sqrt(rel_len2);
    Float local_x = dot(direction, right);
    Float local_y = dot(direction, up);
    Float local_z = dot(direction, forward);
    Float theta = std::acos(clamp(local_y, -1.f, 1.f));
    Float phi = std::atan2(local_x, local_z);
    s = std::fmod((phi - pi_val) / (2.f * pi_val) + 1.f, 1.f);
    t = 1.f - theta / pi_val;
  } else {
    if(z_camera <= 0) {
      return false;
    }
    Float aspect = static_cast<Float>(display_width) /
      static_cast<Float>(display_height);
    Float half_height = std::tan(fov * pi_val / 360.f);
    Float half_width = aspect * half_height;
    s = 0.5f + x_camera / (2.f * z_camera * half_width);
    t = 0.5f + y_camera / (2.f * z_camera * half_height);
  }
  if(!in_front) {
    return false;
  }
  if(overlay.clip && (s < 0 || s > 1 || t < 0 || t > 1)) {
    return false;
  }
  screen_x = (1.f - s) * static_cast<Float>(display_width - 1);
  screen_y = (1.f - t) * static_cast<Float>(display_height - 1);
  return std::isfinite(screen_x) && std::isfinite(screen_y);
}

static bool is_text_pixel_occluded(const PreviewTextOverlay& overlay,
                                   Float screen_x,
                                   Float screen_y,
                                   RayCamera* cam,
                                   hitable* world,
                                   unsigned int display_width,
                                   unsigned int display_height,
                                   random_gen& rng) {
  if(!overlay.partial_occlusion || !cam || !world || cam->get_fov() < 0) {
    return false;
  }
  vec3f forward = cam->get_w();
  Float label_depth = dot(overlay.anchor - cam->get_origin(), forward);
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    label_depth * occlusion_tolerance : occlusion_tolerance;
  if(label_depth <= endpoint_tolerance) {
    return false;
  }

  Float denom_x = static_cast<Float>(std::max(1u, display_width - 1));
  Float denom_y = static_cast<Float>(std::max(1u, display_height - 1));
  Float s = 1.f - screen_x / denom_x;
  Float t = 1.f - screen_y / denom_y;
  Ray visibility_ray = cam->get_ray(s, t, point3f(0.f, 0.f, 0.f), 0.5f);

  hit_record hrec;
  if(!world->hit(visibility_ray, 0.001f, MaxT, hrec, rng) ||
     hrec.shape->GetName() == "EnvironmentLight") {
    return false;
  }
  Float scene_depth = dot(hrec.p - cam->get_origin(), forward);
  return scene_depth < label_depth - endpoint_tolerance;
}

static NumericVector composite_screen_text_overlay(
    const std::vector<PreviewTextOverlay>& overlays,
    unsigned int display_width,
    unsigned int display_height,
    RayCamera* cam,
    hitable* world,
    random_gen& rng) {
  NumericVector output(static_cast<R_xlen_t>(display_width) *
                       static_cast<R_xlen_t>(display_height) * 4);
  output.attr("dim") = IntegerVector::create(
    static_cast<int>(display_height),
    static_cast<int>(display_width),
    4
  );
  if(!cam || !world || display_width == 0 || display_height == 0) {
    return output;
  }
  for(const PreviewTextOverlay& overlay : overlays) {
    if(!overlay.partial_occlusion) {
      continue;
    }
    Float screen_x;
    Float screen_y;
    if(!project_text_anchor_to_screen(overlay,
                                      cam,
                                      display_width,
                                      display_height,
                                      screen_x,
                                      screen_y)) {
      continue;
    }
    int overlay_width = static_cast<int>(overlay.width);
    int overlay_height = static_cast<int>(overlay.height);
    int left = static_cast<int>(std::round(screen_x + overlay.x_offset -
                                           overlay.hjust * overlay_width));
    int top = static_cast<int>(std::round(screen_y + overlay.y_offset -
                                          overlay.vjust * overlay_height));
    int right_bound = std::min<int>(left + overlay_width, display_width);
    int bottom_bound = std::min<int>(top + overlay_height, display_height);
    int x_start = std::max(0, left);
    int y_start = std::max(0, top);
    if(x_start >= right_bound || y_start >= bottom_bound) {
      continue;
    }
    for(int y = y_start; y < bottom_bound; y++) {
      unsigned int source_y = static_cast<unsigned int>(y - top);
      for(int x = x_start; x < right_bound; x++) {
        unsigned int source_x = static_cast<unsigned int>(x - left);
        size_t src_idx = 4 * (source_x + overlay.width * source_y);
        Float source_alpha = overlay.rgba[src_idx + 3] / 255.f;
        if(source_alpha <= 0) {
          continue;
        }
        if(is_text_pixel_occluded(overlay,
                                  static_cast<Float>(x),
                                  static_cast<Float>(y),
                                  cam,
                                  world,
                                  display_width,
                                  display_height,
                                  rng)) {
          continue;
        }
        R_xlen_t dst_idx = y +
          static_cast<R_xlen_t>(display_height) * x;
        R_xlen_t channel_offset =
          static_cast<R_xlen_t>(display_height) * display_width;
        Float dest_alpha = output[dst_idx + 3 * channel_offset];
        Float output_alpha = source_alpha + dest_alpha * (1.f - source_alpha);
        if(output_alpha <= 0) {
          continue;
        }
        for(int channel = 0; channel < 3; channel++) {
          Float source_value = overlay.rgba[src_idx + channel] / 255.f;
          Float dest_value = output[dst_idx + channel * channel_offset];
          output[dst_idx + channel * channel_offset] =
            (source_value * source_alpha +
             dest_value * dest_alpha * (1.f - source_alpha)) / output_alpha;
        }
        output[dst_idx + 3 * channel_offset] = output_alpha;
      }
    }
  }
  return output;
}

static bool project_line_point_to_screen(const point3f& point,
                                         bool clip,
                                         RayCamera* cam,
                                         unsigned int display_width,
                                         unsigned int display_height,
                                         Float& screen_x,
                                         Float& screen_y,
                                         Float& depth) {
  if(!cam || display_width == 0 || display_height == 0) {
    return false;
  }
  Float fov = cam->get_fov();
  if(fov < 0) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = point - origin;
  Float rel_len2 = relative.squared_length();
  if(rel_len2 <= 0) {
    return false;
  }
  vec3f right = cam->get_u();
  vec3f up = cam->get_v();
  vec3f forward = cam->get_w();
  Float x_camera = dot(relative, right);
  Float y_camera = dot(relative, up);
  Float z_camera = dot(relative, forward);
  Float s = 0.5f;
  Float t = 0.5f;
  bool in_front = true;
  const Float pi_val = static_cast<Float>(3.14159265358979323846);

  if(fov == 0) {
    point2f ortho = cam->get_ortho();
    s = 0.5f + x_camera / ortho.xy.x;
    t = 0.5f + y_camera / ortho.xy.y;
    in_front = z_camera >= 0;
  } else if(fov == 360) {
    vec3f direction = relative / std::sqrt(rel_len2);
    Float local_x = dot(direction, right);
    Float local_y = dot(direction, up);
    Float local_z = dot(direction, forward);
    Float theta = std::acos(clamp(local_y, -1.f, 1.f));
    Float phi = std::atan2(local_x, local_z);
    s = std::fmod((phi - pi_val) / (2.f * pi_val) + 1.f, 1.f);
    t = 1.f - theta / pi_val;
  } else {
    if(z_camera <= 0) {
      return false;
    }
    Float aspect = static_cast<Float>(display_width) /
      static_cast<Float>(display_height);
    Float half_height = std::tan(fov * pi_val / 360.f);
    Float half_width = aspect * half_height;
    s = 0.5f + x_camera / (2.f * z_camera * half_width);
    t = 0.5f + y_camera / (2.f * z_camera * half_height);
  }
  if(!in_front) {
    return false;
  }
  if(clip && (s < 0 || s > 1 || t < 0 || t > 1)) {
    return false;
  }
  screen_x = (1.f - s) * static_cast<Float>(display_width - 1);
  screen_y = (1.f - t) * static_cast<Float>(display_height - 1);
  depth = z_camera;
  return std::isfinite(screen_x) && std::isfinite(screen_y);
}

static Float line_coverage_and_depth(Float x0,
                                     Float y0,
                                     Float depth0,
                                     Float x1,
                                     Float y1,
                                     Float depth1,
                                     Float px,
                                     Float py,
                                     Float width,
                                     int lineend,
                                     Float& line_depth) {
  Float radius = width / 2.f;
  Float dx = x1 - x0;
  Float dy = y1 - y0;
  Float len2 = dx * dx + dy * dy;
  if(len2 <= static_cast<Float>(1e-12)) {
    Float dist = std::sqrt((px - x0) * (px - x0) + (py - y0) * (py - y0));
    line_depth = (depth0 + depth1) * 0.5f;
    return clamp(radius + 0.5f - dist, 0.f, 1.f);
  }

  Float t_depth = ((px - x0) * dx + (py - y0) * dy) / len2;
  line_depth = depth0 + clamp(t_depth, 0.f, 1.f) * (depth1 - depth0);
  Float cx0 = x0;
  Float cy0 = y0;
  Float cx1 = x1;
  Float cy1 = y1;
  if(lineend == 2) {
    Float len = std::sqrt(len2);
    Float ux = dx / len;
    Float uy = dy / len;
    cx0 -= ux * radius;
    cy0 -= uy * radius;
    cx1 += ux * radius;
    cy1 += uy * radius;
    dx = cx1 - cx0;
    dy = cy1 - cy0;
    len2 = dx * dx + dy * dy;
  }

  Float t = ((px - cx0) * dx + (py - cy0) * dy) / len2;
  if(lineend == 1 && (t < 0.f || t > 1.f)) {
    return 0.f;
  }
  t = clamp(t, 0.f, 1.f);
  Float closest_x = cx0 + t * dx;
  Float closest_y = cy0 + t * dy;
  Float dist = std::sqrt((px - closest_x) * (px - closest_x) +
                         (py - closest_y) * (py - closest_y));
  return clamp(radius + 0.5f - dist, 0.f, 1.f);
}

static bool is_line_pixel_occluded(const PreviewLineOverlay& overlay,
                                   Float screen_x,
                                   Float screen_y,
                                   Float line_depth,
                                   RayCamera* cam,
                                   hitable* world,
                                   unsigned int display_width,
                                   unsigned int display_height,
                                   random_gen& rng) {
  if(!overlay.partial_occlusion || !cam || !world || cam->get_fov() < 0) {
    return false;
  }
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    line_depth * occlusion_tolerance : occlusion_tolerance;
  if(line_depth <= endpoint_tolerance) {
    return false;
  }
  Float denom_x = static_cast<Float>(std::max(1u, display_width - 1));
  Float denom_y = static_cast<Float>(std::max(1u, display_height - 1));
  Float s = 1.f - screen_x / denom_x;
  Float t = 1.f - screen_y / denom_y;
  Ray visibility_ray = cam->get_ray(s, t, point3f(0.f, 0.f, 0.f), 0.5f);

  hit_record hrec;
  if(!world->hit(visibility_ray, 0.001f, MaxT, hrec, rng) ||
     hrec.shape->GetName() == "EnvironmentLight") {
    return false;
  }
  Float scene_depth = dot(hrec.p - cam->get_origin(), cam->get_w());
  return scene_depth < line_depth - endpoint_tolerance;
}

static NumericVector composite_screen_line_overlay(
    const std::vector<PreviewLineOverlay>& overlays,
    unsigned int display_width,
    unsigned int display_height,
    RayCamera* cam,
    hitable* world,
    random_gen& rng) {
  NumericVector output(static_cast<R_xlen_t>(display_width) *
                       static_cast<R_xlen_t>(display_height) * 4);
  output.attr("dim") = IntegerVector::create(
    static_cast<int>(display_height),
    static_cast<int>(display_width),
    4
  );
  if(!cam || !world || display_width == 0 || display_height == 0) {
    return output;
  }
  for(const PreviewLineOverlay& overlay : overlays) {
    if(!overlay.partial_occlusion) {
      continue;
    }
    Float x0;
    Float y0;
    Float depth0;
    Float x1;
    Float y1;
    Float depth1;
    if(!project_line_point_to_screen(overlay.start,
                                     overlay.clip,
                                     cam,
                                     display_width,
                                     display_height,
                                     x0,
                                     y0,
                                     depth0) ||
       !project_line_point_to_screen(overlay.end,
                                     overlay.clip,
                                     cam,
                                     display_width,
                                     display_height,
                                     x1,
                                     y1,
                                     depth1)) {
      continue;
    }
    x0 += overlay.x_offset;
    y0 += overlay.y_offset;
    x1 += overlay.xend_offset;
    y1 += overlay.yend_offset;
    Float pad = overlay.width / 2.f + 1.f;
    int left = std::max(0, static_cast<int>(std::floor(std::min(x0, x1) - pad)));
    int right = std::min<int>(display_width - 1,
                              static_cast<int>(std::ceil(std::max(x0, x1) + pad)));
    int top = std::max(0, static_cast<int>(std::floor(std::min(y0, y1) - pad)));
    int bottom = std::min<int>(display_height - 1,
                               static_cast<int>(std::ceil(std::max(y0, y1) + pad)));
    if(left > right || top > bottom) {
      continue;
    }
    R_xlen_t channel_offset =
      static_cast<R_xlen_t>(display_height) * display_width;
    for(int y = top; y <= bottom; y++) {
      for(int x = left; x <= right; x++) {
        Float line_depth;
        Float source_alpha = overlay.alpha * line_coverage_and_depth(
          x0, y0, depth0, x1, y1, depth1,
          static_cast<Float>(x), static_cast<Float>(y),
          overlay.width, overlay.lineend, line_depth
        );
        if(source_alpha <= 0) {
          continue;
        }
        if(is_line_pixel_occluded(overlay,
                                  static_cast<Float>(x),
                                  static_cast<Float>(y),
                                  line_depth,
                                  cam,
                                  world,
                                  display_width,
                                  display_height,
                                  rng)) {
          continue;
        }
        R_xlen_t dst_idx = y +
          static_cast<R_xlen_t>(display_height) * x;
        Float dest_alpha = output[dst_idx + 3 * channel_offset];
        Float output_alpha = source_alpha + dest_alpha * (1.f - source_alpha);
        if(output_alpha <= 0) {
          continue;
        }
        Float source_values[3] = {overlay.red, overlay.green, overlay.blue};
        for(int channel = 0; channel < 3; channel++) {
          Float dest_value = output[dst_idx + channel * channel_offset];
          output[dst_idx + channel * channel_offset] =
            (source_values[channel] * source_alpha +
             dest_value * dest_alpha * (1.f - source_alpha)) / output_alpha;
        }
        output[dst_idx + 3 * channel_offset] = output_alpha;
      }
    }
  }
  return output;
}


// [[Rcpp::export]]
List render_scene_rcpp(List scene, List camera_info, List scene_info, List render_info) {
  RESET_RAYLOG();
  START_TIMER("Overall Time");
  feclearexcept(FE_ALL_EXCEPT);

  //Unpack scene info
  IntegerVector shape = as<IntegerVector>(scene_info["shape"]);
  
  //Unpack render info
  bool ambient_light = as<bool>(render_info["ambient_light"]);
  NumericVector bghigh  = as<NumericVector>(render_info["bghigh"]);
  NumericVector bglow = as<NumericVector>(render_info["bglow"]);
  Float clampval = as<Float>(render_info["clampval"]);
  bool progress_bar = as<bool>(render_info["progress_bar"]);
  int numbercores = as<int>(render_info["numbercores"]);
  bool hasbackground = as<bool>(render_info["hasbackground"]);
  std::string background = as<std::string>(render_info["background"]);
  Float rotate_env = as<Float>(render_info["rotate_env"]);
  Float intensity_env = as<Float>(render_info["intensity_env"]);
  bool verbose = as<bool>(render_info["verbose"]);
  int debug_channel = as<int>(render_info["debug_channel"]);
  Float min_variance = as<Float>(render_info["min_variance"]);
  int min_adaptive_size = as<int>(render_info["min_adaptive_size"]);
  IntegratorType integrator_type = static_cast<IntegratorType>(as<int>(render_info["integrator_type"]));
  bool print_debug_info = as<bool>(render_info["print_debug_info"]);
#ifdef HAS_OIDN
  bool denoise = as<bool>(render_info["denoise"]);
#endif

  Environment pkg = Environment::namespace_env("rayrender");
  Function print_time = pkg["print_time"];


  //Unpack Camera Info
  int nx = as<int>(camera_info["nx"]);
  int ny = as<int>(camera_info["ny"]);
  int ns = as<int>(camera_info["ns"]);
  Float fov = as<Float>(camera_info["fov"]);
  NumericVector lookfromvec = as<NumericVector>(camera_info["lookfrom"]);
  NumericVector lookatvec = as<NumericVector>(camera_info["lookat"]);
  Float aperture = as<Float>(camera_info["aperture"]);
  NumericVector camera_up = as<NumericVector>(camera_info["camera_up"]);
  Float shutteropen = as<Float>(camera_info["shutteropen"]);
  Float shutterclose = as<Float>(camera_info["shutterclose"]);
  Float focus_distance = as<Float>(camera_info["focal_distance"]);
  NumericVector ortho_dimensions = as<NumericVector>(camera_info["ortho_dimensions"]);
  std::size_t max_depth = as<std::size_t>(camera_info["max_depth"]);
  std::size_t roulette_active = as<std::size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  NumericVector light_direction = as<NumericVector>(camera_info["light_direction"]);
  int stratified_x = static_cast<int>(stratified_dim(0));
  int stratified_y = static_cast<int>(stratified_dim(1));
  vec3f preview_light_direction(light_direction(0), light_direction(1), light_direction(2));
  Float preview_exponent = light_direction.size() > 3 ? static_cast<Float>(light_direction(3)) : 0;
  NumericMatrix realCameraInfo = as<NumericMatrix>(camera_info["real_camera_info"]);
  Float film_size = as<Float>(camera_info["film_size"]);
  Float camera_scale = as<Float>(camera_info["camera_scale"]);
  Float sample_dist = as<Float>(camera_info["sample_dist"]);
  bool keep_colors = as<bool>(camera_info["keep_colors"]);
  bool preview     = as<bool>(camera_info["preview"]);
  bool interactive = as<bool>(camera_info["interactive"]);
  bool deferred_render = as<bool>(camera_info["deferred_render"]);
  bool auto_exposure = as<bool>(camera_info["auto_exposure"]);
  Float iso = as<Float>(camera_info["iso"]);
  int bvh_type = as<int>(camera_info["bvh"]);

  
  //Initialize transformation cache
  TransformCache transformCache;
  TransformCache transformCacheBg;
  
  //Initialize texture cache
  TextureCache texCache;

  //Initialize output matrices
  RayMatrix rgb_output(nx,ny, 3);
  RayMatrix draw_rgb_output(nx,ny, 3);
  RayMatrix alpha_output(nx,ny, 1);
  RayMatrix normalOutput(nx,ny, 3);
  RayMatrix albedoOutput(nx,ny, 3);

#ifdef HAS_OIDN
  // Create an Open Image Denoise device
  oidn::DeviceRef device = oidn::newDevice(); // CPU or GPU if available
  // oidn::DeviceRef device = oidn::newDevice(oidn::DeviceType::CPU);
  device.commit();
  // Create buffers for input/output images accessible by both host (CPU) and device (CPU/GPU)
  oidn::BufferRef colorBuf  = device.newBuffer(rgb_output.begin(), nx * ny * 3 * sizeof(Float));
  oidn::BufferRef albedoBuf = device.newBuffer(albedoOutput.begin(), nx * ny * 3 * sizeof(Float));
  oidn::BufferRef normalBuf = device.newBuffer(normalOutput.begin(), nx * ny * 3 * sizeof(Float));
  oidn::BufferRef colorBuf2 = device.newBuffer(draw_rgb_output.begin(), nx * ny * 3 * sizeof(Float));

  // Create a filter for denoising a beauty (color) image using optional auxiliary images too
  // This can be an expensive operation, so try no to create a new filter for every image!
  oidn::FilterRef filter = device.newFilter("RT"); // generic ray tracing filter
  filter.setImage("color",  colorBuf,  oidn::Format::Float3, nx, ny); // beauty
  filter.setImage("albedo", albedoBuf, oidn::Format::Float3, nx, ny); // auxiliary
  filter.setImage("normal", normalBuf, oidn::Format::Float3, nx, ny); // auxiliary
  filter.setImage("output", colorBuf2,  oidn::Format::Float3, nx, ny); // denoised beauty
  filter.set("hdr", true); // beauty image is HDR
  filter.commit();
#endif
  
  point3f lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  point3f lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  point3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  point3f backgroundlow(bglow[0],bglow[1],bglow[2]);
  Float dist_to_focus = focus_distance;
  
  std::vector<bool> has_image;
  std::vector<bool> has_alpha;
  std::vector<bool> has_bump;
  std::vector<bool> has_roughness;
  
  std::unique_ptr<RayCamera> cam;


  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  if(fov < 0) {
    Transform CamTransform = LookAt(lookfrom,
                                    lookat,
                                    vec3f(camera_up(0),camera_up(1),camera_up(2))).GetInverseMatrix();
    Transform* CameraTransform = transformCache.Lookup(CamTransform);
    
    AnimatedTransform CamTr(CameraTransform,0,CameraTransform,0);
    
    std::vector<Float> lensData;
    for(int i = 0; i < realCameraInfo.rows(); i++) {
      for(int j = 0; j < realCameraInfo.cols(); j++) {
        lensData.push_back(realCameraInfo.at(i,j));
      }
    }
    
    if(fov < 0 && lensData.size() == 0) {
      throw std::runtime_error("No lens data passed in lens descriptor file.");
    }
    
    cam = std::unique_ptr<RayCamera>(new RealisticCamera(CamTr,shutteropen, shutterclose,
                         aperture, nx,ny, focus_distance, false, lensData,
                         film_size, camera_scale, iso, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                         CamTransform, lookat));
  } else if(fov == 0) {
    cam = std::unique_ptr<RayCamera>(new ortho_camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                      ortho_dimensions(0), ortho_dimensions(1),
                      shutteropen, shutterclose, iso));
  } else if (fov == 360) {
    cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                            shutteropen, shutterclose, iso));
  } else {
    cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)), fov, Float(nx)/Float(ny),
                                     aperture, dist_to_focus,
                                     shutteropen, shutterclose, iso));
  }
  print_time(verbose, "Generated Camera" );


  int nx1, ny1, nn1;

  std::vector<Float* > textures;
  std::vector<unsigned char * > alpha_textures;
  std::vector<unsigned char * > bump_textures;
  std::vector<unsigned char * > roughness_textures;
  
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;

  
  hitable_list imp_sample_objects;
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  std::vector<std::shared_ptr<alpha_texture> > alpha;
  std::vector<std::shared_ptr<bump_texture> > bump;
  std::vector<std::shared_ptr<roughness_texture> > roughness;
  std::vector<int> texture_idx;

  std::shared_ptr<hitable> worldbvh = build_scene(scene, 
                                                   shape, 
                                                   shutteropen,
                                                   shutterclose,
                                                   textures, 
                                                   alpha_textures,
                                                   bump_textures,
                                                   roughness_textures, 
                                                   shared_materials, 
                                                   alpha, bump, roughness,
                                                   bvh_type,
                                                   transformCache, 
                                                   texCache,
                                                   imp_sample_objects,
                                                   instanced_objects,
                                                   instance_importance_sampled,
                                                   texture_idx,
                                                   verbose, 
                                                   rng);
  print_time(verbose, "Built Scene BVH" );
  if(print_debug_info) {
    worldbvh->hitable_info_bounds(shutteropen,shutterclose);
  }
  //Calculate world bounds and ensure camera is inside infinite area light
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = convert_to_vec3(bounding_box_world.Centroid());
  world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : 
     1.1*(lookfrom - world_center).length();
  world_radius *= interactive ? 100 : 1;
  
  if(fov == 0) {
    Float ortho_diag = sqrt(pow(ortho_dimensions(0),2) + pow(ortho_dimensions(1),2));
    world_radius += ortho_diag;
  }

  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;
  Float *background_texture_data = nullptr;

  //Background rotation
  Matrix4x4 Identity;
  Transform BackgroundAngle(Identity);
  if(rotate_env != 0) {
    BackgroundAngle = Translate(world_center) * RotateY(rotate_env);
  } else {
    BackgroundAngle = Translate(world_center);
  }

  Transform* BackgroundTransform = transformCacheBg.Lookup(BackgroundAngle);
  Transform* BackgroundTransformInv = transformCacheBg.Lookup(BackgroundAngle.GetInverseMatrix());
  if(hasbackground) {
    background_texture_data = texCache.LookupFloat(background, nx1, ny1, nn1, 3);
    // nn1 = 3;
    // texture_bytes += nx1 * ny1 * nn1;
    
    if(background_texture_data) {
      bool has_env_light = false;
      const std::size_t env_size = static_cast<std::size_t>(nx1) *
        static_cast<std::size_t>(ny1) *
        static_cast<std::size_t>(nn1);
      for(std::size_t i = 0; i < env_size; i++) {
        if(background_texture_data[i] > 0) {
          has_env_light = true;
          break;
        }
      }
      if(has_env_light) {
        background_texture = std::make_shared<image_texture_float>(background_texture_data, nx1, ny1, nn1,
                                                                   1, 1, intensity_env);
        background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
        background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, convert_to_point3(world_center),
                                                                background_texture, background_material,
                                                                BackgroundTransform,
                                                                BackgroundTransformInv, false);
      } else {
        hasbackground = false;
        ambient_light = true;
        backgroundhigh = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
        backgroundlow = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
        background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
        background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
        background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                                                background_texture, background_material,
                                                                BackgroundTransform,BackgroundTransformInv,false);
      }
    } else {
      Rcpp::Rcout << "Failed to load background image at " << background << "\n";
      hasbackground = false;
      ambient_light = true;
      backgroundhigh = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                                              background_texture, background_material,
                                                              BackgroundTransform,BackgroundTransformInv,false);
    }
  } else if(ambient_light) {
    //Check if both high and low are black, and set to FLT_MIN
    if(backgroundhigh.length() == 0 && backgroundlow.length() == 0) {
      backgroundhigh = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
    }
    background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                              background_texture, background_material,
                                              BackgroundTransform, BackgroundTransformInv, false);

  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(point3f(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                              background_texture, background_material,
                                              BackgroundTransform,
                                              BackgroundTransformInv, false);
  }
  print_time(verbose, "Loaded background" );
  hitable_list world;
  world.add(worldbvh);

  bool impl_only_bg = false;
  world.add(background_sphere);
  if((imp_sample_objects.size() == 0 || hasbackground || ambient_light || interactive) && debug_channel != 18) {
    impl_only_bg = true;
  }
  LogicalVector screen_text_visible = compute_screen_text_visibility(
    render_info,
    cam.get(),
    worldbvh.get(),
    rng
  );
  LogicalVector screen_line_visible = compute_screen_line_visibility(
    render_info,
    cam.get(),
    worldbvh.get(),
    rng
  );
  std::vector<PreviewTextOverlay> text_overlays =
    parse_preview_text_overlays(render_info);
  std::vector<PreviewLineOverlay> line_overlays =
    parse_preview_line_overlays(render_info);
  preview = preview && debug_channel == 0;
#ifdef HAS_OIDN
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         deferred_render, (lookat-lookfrom).length(), cam.get(),
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         filter, denoise, auto_exposure);
#else
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         deferred_render, (lookat-lookfrom).length(), cam.get(),
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         auto_exposure);
#endif
  Display.SetTextOverlays(text_overlays);
  Display.SetLineOverlays(line_overlays);
  
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  QUERY_MEMORY_USAGE();
  PRINT_CURRENT_MEMORY("Before raytracing");
  
  // Rcpp::Rcout << "Total world size: " << world.GetSize() + texture_bytes << " (Textures: " << texture_bytes << ") \n";
  if(debug_channel != 0) {
    debug_scene(numbercores, nx, ny, ns, debug_channel,
                min_variance, min_adaptive_size,
                rgb_output, normalOutput, albedoOutput,
                progress_bar, sample_method, stratified_x, stratified_y,
                verbose, cam.get(), fov,
                world, imp_sample_objects, 
                clampval, max_depth, roulette_active,
                preview_light_direction, preview_exponent, rng, sample_dist, keep_colors,
                backgroundhigh);
  } else {
    pathtracer(numbercores, nx, ny, ns, debug_channel,
               min_variance, min_adaptive_size,
               rgb_output, normalOutput, albedoOutput,
               alpha_output,
               draw_rgb_output,
               progress_bar, sample_method, stratified_x, stratified_y,
               verbose, cam.get(),  fov,
               world, imp_sample_objects,
               clampval, max_depth, roulette_active, Display, integrator_type);
  }
  PRINT_CURRENT_MEMORY("After raytracing");
#ifdef HAS_OIDN
  if(denoise) {
    filter.execute();
    const char* errorMessage;
    if (device.getError(errorMessage) != oidn::Error::None) {
      Rcpp::Rcout << "Error: " << errorMessage << std::endl;
    }
  }
#endif
  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  RayMatrix final_output = rgb_output;
#ifdef HAS_OIDN
  if(denoise) {
    final_output = draw_rgb_output;
  }
#endif
  #ifdef RAY_COLOR_DEBUG
  final_output.print();
  #endif
  List final_image = List::create(_["r"] = final_output.ConvertRcpp(0), 
                                  _["g"] = final_output.ConvertRcpp(1), 
                                  _["b"] = final_output.ConvertRcpp(2),

                                  _["nx"] = normalOutput.ConvertRcpp(0), 
                                  _["ny"] = normalOutput.ConvertRcpp(1), 
                                  _["nz"] = normalOutput.ConvertRcpp(2),

                                  _["cx"] = albedoOutput.ConvertRcpp(0), 
                                  _["cy"] = albedoOutput.ConvertRcpp(1), 
                                  _["cz"] = albedoOutput.ConvertRcpp(2),

                                  _["a"] = alpha_output.ConvertRcpp());
  if(Display.Keyframes.size() > 0) {
    List keyframes(Display.Keyframes.size());
    for(unsigned int i = 0; i < Display.Keyframes.size(); i++ ) {
      keyframes(i) = Display.Keyframes[i];
    }
    final_image.attr("keyframes") = keyframes;
  }
  final_image.attr("preview_exposure") = Display.preview_exposure_adjustment;
  final_image.attr("screen_camera_info") = get_screen_camera_info(cam.get());
  if(screen_text_visible.size() > 0) {
    final_image.attr("screen_text_visible") = screen_text_visible;
  }
  if(screen_line_visible.size() > 0) {
    final_image.attr("screen_line_visible") = screen_line_visible;
  }
  if(should_return_screen_text_overlay(render_info)) {
    final_image.attr("screen_text_overlay") = composite_screen_text_overlay(
      text_overlays,
      static_cast<unsigned int>(nx),
      static_cast<unsigned int>(ny),
      cam.get(),
      worldbvh.get(),
      rng
    );
  }
  if(should_return_screen_line_overlay(render_info)) {
    final_image.attr("screen_line_overlay") = composite_screen_line_overlay(
      line_overlays,
      static_cast<unsigned int>(nx),
      static_cast<unsigned int>(ny),
      cam.get(),
      worldbvh.get(),
      rng
    );
  }
  STOP_TIMER("Overall Time");

  PRINT_LOG_REPORT(numbercores);

  PRINT_CURRENT_MEMORY("After cleanup");
  
  return(final_image);
}
