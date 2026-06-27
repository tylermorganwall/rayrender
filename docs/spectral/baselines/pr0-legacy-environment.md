# PR 0 Legacy Baseline Environment

- Captured at: `2026-06-26T19:48:39-0400`
- rayrender package version: `0.41.3`
- rayrender commit: `0ab9b6a2de4b075f91daefee6d38487782cabd9e`
- pbrt-v4 commit: `8c19f304558fd7681e2fef2c395a689d0106fb05`
- R: `R version 4.6.0 (2026-04-24)`
- R platform: `aarch64-apple-darwin23`
- System: `Darwin 24.6.0 arm64`
- CXX: `clang++ -arch arm64 -g -std=gnu++17`
- CXXFLAGS: `-falign-functions=64 -Wall -g -O2`
- CPPFLAGS: `-I/opt/R/arm64/include`

## Captured Scenes

- `cornell_diffuse_still`: render_scene; cornell diffuse sphere; random; nee; denoise=FALSE; parallel=FALSE
- `material_mix_still`: render_scene; mixed diffuse/metal/cube; random; nee; denoise=FALSE; parallel=FALSE
- `dielectric_still`: render_scene; dielectric plus diffuse sphere; random; nee; denoise=FALSE; parallel=FALSE
- `orbit_animation_frames`: render_animation; two orbit frames; random; nee; denoise=FALSE; parallel=FALSE
