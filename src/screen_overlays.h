#ifndef RAYRENDER_SCREEN_OVERLAYS_H
#define RAYRENDER_SCREEN_OVERLAYS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "core/camera.h"
#include "core/PreviewDisplay.h"
#include "hitables/hitable.h"
#include "math/rng.h"

// Shared by still and animation rendering so screen placement, pixel widths,
// and occlusion use the same camera and world-space conventions.
namespace rayrender_screen {
using namespace Rcpp;

inline std::vector<PreviewTextOverlay> parse_preview_text_overlays(List render_info) {
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

inline std::vector<PreviewLineOverlay> parse_preview_line_overlays(List render_info) {
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

inline bool is_text_anchor_occluded(const point3f& anchor,
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

inline LogicalVector compute_screen_text_visibility(List render_info,
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

inline LogicalVector compute_screen_line_visibility(List render_info,
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

inline bool should_return_screen_text_overlay(List render_info) {
  if(!render_info.containsElementNamed("screen_text_native_overlay")) {
    return false;
  }
  return as<bool>(render_info["screen_text_native_overlay"]);
}

inline bool should_return_screen_line_overlay(List render_info) {
  if(!render_info.containsElementNamed("screen_line_native_overlay")) {
    return false;
  }
  return as<bool>(render_info["screen_line_native_overlay"]);
}

inline List get_screen_camera_info(RayCamera* cam) {
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

inline bool project_text_anchor_to_screen(const PreviewTextOverlay& overlay,
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

inline bool is_text_pixel_occluded(const PreviewTextOverlay& overlay,
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

inline NumericVector composite_screen_text_overlay(
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
    if(overlay.occlusion &&
       is_text_anchor_occluded(overlay.anchor,
                               overlay.occlusion_tolerance,
                               cam,
                               world,
                               rng)) {
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

inline bool project_line_point_to_screen(const point3f& point,
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

inline Float line_coverage_and_depth(Float x0,
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

inline bool is_line_pixel_occluded(const PreviewLineOverlay& overlay,
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

inline NumericVector composite_screen_line_overlay(
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
    if(overlay.occlusion) {
      point3f anchor = overlay.start + (overlay.end - overlay.start) * 0.5f;
      if(is_text_anchor_occluded(anchor,
                                 overlay.occlusion_tolerance,
                                 cam,
                                 world,
                                 rng)) {
        continue;
      }
    }
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


} // namespace rayrender_screen

#endif
