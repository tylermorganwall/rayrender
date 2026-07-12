#ifndef RAYRENDER_OIDN_AUX_H
#define RAYRENDER_OIDN_AUX_H

#include <cstddef>
#include <functional>

#include "../core/camera.h"
#include "../hitables/hitablelist.h"
#include "../math/RayMatrix.h"

struct OidnAuxRenderOptions {
  std::size_t samples = 1;
  std::size_t max_depth = 8;
  std::size_t max_dielectric_splits = 1;
  int sample_method = 0;
  int stratified_x = 1;
  int stratified_y = 1;
};

void render_oidn_aux_features(std::size_t numbercores,
                              std::size_t nx,
                              std::size_t ny,
                              RayCamera* cam,
                              Float fov,
                              hitable* world,
                              const OidnAuxRenderOptions& options,
                              RayMatrix& normalOutput,
                              RayMatrix& albedoOutput,
                              std::function<bool()> poll_cancel = std::function<bool()>());

#endif
