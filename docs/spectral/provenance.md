# Spectral Data and Source Provenance

This file tracks copied or generated spectral assets and source-derived implementation material. PR 0 adds the inventory and notice policy only; it does not copy pbrt source, pbrt data tables, RGB-to-spectrum tables, or measured spectral datasets into rayrender.

## Current Inventory

| Asset or source | Location in rayrender | Origin | License | Status | Required notice |
|---|---|---|---|---|---|
| pbrt-v4 reference source | `mmp/pbrt-v4` | `https://github.com/mmp/pbrt-v4.git` at `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Reference only; not copied into package | Preserve copyright and Apache-2.0 notices before copying or adapting source |
| pbrt-v4 sampled spectrum and wavelength utilities | `src/base/sampled_spectrum.h`; `tools/spectral-tests/pr3-spectrum-tests.cpp` | Adapted from `src/pbrt/util/sampling.h` and `src/pbrt/util/spectrum.h` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Source-derived implementation material introduced in PR 3 | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 spectrum representations and setup helpers | `src/base/spectrum.h`; `tools/spectral-tests/pr4-spectrum-assets-tests.cpp` | Adapted from `src/pbrt/util/spectrum.h` and `src/pbrt/util/spectrum.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Source-derived implementation material introduced in PR 4 | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 named spectral data | `inst/extdata/spectral/named-spectra-v1.tsv`; `inst/extdata/spectral/named-spectra-v1-metadata.tsv`; `tools/spectral-tests/generate-pr4-spectral-assets.R` | Generated from `src/pbrt/util/spectrum.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`; glass source rows originate from pbrt's embedded refractiveindex.info CC0 data | Apache-2.0 and CC0-1.0 | Packaged named spectra introduced in PR 4; file checksums recorded in `docs/spectral/assets-manifest.csv`; per-spectrum checksums recorded in metadata TSV | `inst/COPYRIGHTS` records pbrt and refractiveindex.info notices |
| pbrt-v4 RGB color space and reconstruction code | `src/base/color_types.h`; `src/base/rgb_spectrum.h`; `tools/spectral-tests/pr5-rgb-spectrum-tests.cpp` | Adapted from `src/pbrt/util/color.h`, `src/pbrt/util/color.cpp`, `src/pbrt/util/colorspace.h`, and `src/pbrt/util/colorspace.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Source-derived implementation material introduced in PR 5 | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 sRGB RGB-to-spectrum coefficient table | `inst/extdata/spectral/rgb-to-spectrum-srgb-v1.bin`; `tools/spectral-tests/generate-pr5-rgb-table.R` | Generated from pinned pbrt build output `build/rgbspectrum_srgb.cpp`, which is produced by pbrt's RGB-to-spectrum table generation algorithm at commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Packaged binary table introduced in PR 5; SHA-256 checksum recorded in `docs/spectral/assets-manifest.csv`; binary header stores version, dimensions, color-space id, scalar precision, byte-order marker, and payload Adler-32 | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 Film and PixelSensor contracts | `src/render/spectral_film.h`; `src/render/spectral_film.cpp`; `tools/spectral-tests/pr8-film-tests.cpp` | Adapted from `src/pbrt/film.h` and `src/pbrt/film.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Source-derived implementation material introduced in PR 8; no pbrt ColorChecker swatch data or measured sensor datasets are copied | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 camera contracts | `src/render/spectral_camera.h`; `src/render/spectral_camera.cpp`; `tools/spectral-tests/pr9-camera-tests.cpp` | Adapted from `src/pbrt/base/camera.h`, `src/pbrt/base/sampler.h`, `src/pbrt/cameras.h`, and `src/pbrt/cameras.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`, with rayrender legacy camera frame conventions preserved | Apache-2.0 | Source-derived implementation material introduced in PR 9; no pbrt lens data, realistic-camera tables, or measured camera datasets are copied | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 interaction, shape, and primitive contracts | `src/render/spectral_scene.h`; `src/render/spectral_scene.cpp`; `src/render/spectral_legacy_shape.cpp`; `tools/spectral-tests/pr10-scene-tests.cpp` | Adapted from `src/pbrt/interaction.h`, `src/pbrt/ray.h`, `src/pbrt/shapes.h`, and `src/pbrt/cpu/primitive.cpp` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`, with rayrender legacy hit-record fields and CSG extension constraints preserved | Apache-2.0 | Source-derived implementation material introduced in PR 10; no pbrt geometry assets or third-party shape data are copied | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 light and light-sampler contracts | `src/render/spectral_light.h`; `src/render/spectral_light.cpp`; `tools/spectral-tests/pr11-light-tests.cpp` | Adapted from `src/pbrt/base/light.h`, `src/pbrt/lights.h`, `src/pbrt/lights.cpp`, and `src/pbrt/lightsamplers.h` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`, with rayrender schema-v2 adapter hooks and legacy environment-sphere exclusion preserved | Apache-2.0 | Source-derived implementation material introduced in PR 11; no image environment assets or measured light data are copied | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 BxDF/BSDF and scattering contracts | `src/render/spectral_bsdf.h`; `src/render/spectral_bsdf.cpp`; `tools/spectral-tests/pr12-bsdf-tests.cpp` | Adapted from `src/pbrt/base/bxdf.h`, `src/pbrt/bxdfs.h`, `src/pbrt/bsdf.h`, and `src/pbrt/util/scattering.h` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`, with only the standalone diffuse/null foundation, Fresnel, Trowbridge-Reitz, and BSDF frame wrapper ported in PR 12 | Apache-2.0 | Source-derived implementation material introduced in PR 12; no measured optical constants or material assets are copied | `inst/COPYRIGHTS` records pbrt notice |
| pbrt-v4 Material closure contracts | `src/materials/spectral_material.h`; `src/materials/spectral_material.cpp`; `tools/spectral-tests/pr13-material-tests.cpp` | Adapted from `src/pbrt/base/material.h` and `src/pbrt/materials.h` at pinned pbrt commit `8c19f304558fd7681e2fef2c395a689d0106fb05`, with rayrender texture/context and `ScratchBuffer` APIs preserved | Apache-2.0 | Source-derived implementation material introduced in PR 13; no material asset data is copied | `inst/COPYRIGHTS` records pbrt notice |
| Spectral implementation plan | `docs/spectral/rayrender_spectral_rendering_codex_plan.md` | Project-local plan document | Project documentation | Present before PR 0 implementation | None beyond repository license policy |

## Package Size Impact

PR 5 adds `inst/extdata/spectral/rgb-to-spectrum-srgb-v1.bin`, a 9.0 MiB installed-package asset. The table keeps pbrt's `64^3` coefficient resolution and stores float32 scale nodes plus `[3][64][64][64][3]` coefficients in a compact versioned binary container. Additional RGB color spaces are deferred so PR 5 adds only the sRGB table.

## Expected Future Asset Records

Every copied or generated spectral asset must have a row in `docs/spectral/assets-manifest.csv` before it is used by implementation code. The row must record:

- repository path;
- semantic type;
- source citation or generator;
- source license;
- SHA-256 checksum;
- whether it must be included in the source package;
- any required notice or `inst/COPYRIGHTS` update.

## Notice Policy

Before copying pbrt source, pbrt-generated tables, measured spectra, or third-party sensor/lens data:

1. confirm the source license and compatibility with rayrender's GPL-3 package license;
2. preserve SPDX, copyright, and Apache-2.0 notice text where required;
3. record generator script versions and input checksums for generated binary assets;
4. update `inst/COPYRIGHTS` in the same PR that introduces the asset;
5. verify packaged assets with `tools/spectral-tests/check-spectral-assets.R`.

This document is an engineering provenance record, not legal advice.
