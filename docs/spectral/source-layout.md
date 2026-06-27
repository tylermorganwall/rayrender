# Renderer Source Layout

PR 1 introduces a shared render-session layer without changing RGB transport. PR 2 through PR 5 add
isolated base types, packaged spectral assets, and pbrt-style sRGB spectrum reconstruction. PR 6 adds
the R schema-v2 descriptor and compiler-input validation shell. PR 7 adds the first C++ spectral
texture layer. PR 8 adds the spectral PixelSensor and Film accumulation layer. PR 9 adds the
wavelength-aware spectral camera interface. PR 10 adds spectral shape, primitive, aggregate, and
scene scaffolding. PR 11 adds explicit spectral lights, infinite lights, light registries, and
light samplers.

## Current Layout

- `src/render/`: renderer setup orchestration shared by still and animation entry points.
- `src/render/render_session.*`: legacy RGB setup bridge for scene compilation, camera construction, environment light setup, output buffers, and explicit frame inputs.
- `src/render/spectral_film.*`: PR 8 PixelSensor, Film, FilmTile, output-color conversion, visible-surface hooks, and deterministic test accumulation.
- `src/render/spectral_camera.*`: PR 9 spectral-only camera sample, ray, handle, perspective/orthographic/thin-lens generation, ray differentials, initial-medium recording, dielectric-region initialization stubs, and Film wavelength-sampling bridge.
- `src/render/spectral_scene.*`: PR 10 spectral Interaction, SurfaceInteraction, Shape callbacks, native sphere and SDF CSG adapters, primitive bindings, transformed primitives, aggregates, Scene ownership, ray spawning, capability validation, and legacy hitable wrapping.
- `src/render/spectral_light.*`: PR 11 spectral Light flags/types, point/spot/distant/diffuse-area lights, uniform/image infinite lights, visibility endpoints, scalar light power estimates, UniformLightSampler/PowerLightSampler, and schema/compiler adapter hooks for area lights, free lights, environments, and legacy emissive materials.
- `src/base/`: explicit RGB/XYZ/color-space types, sRGB transfer functions, RGB-to-spectrum table loading, RGB spectrum wrappers, `SampledSpectrum`, `SampledWavelengths`, spectrum representations and named-spectrum registry loading, tagged dispatch handles, `ScratchBuffer`, BxDF flags, and optional sample-result conventions for future pbrt-style transport.
- `src/materials/spectral_texture.*`: PR 7 spectral texture descriptors, evaluators, semantic image cache keys, image loading adapters, and hit-record derivative adapters.
- `R/schema_v2_descriptors.R`: PR 6 schema-v2 spectral descriptors, scene decorators, validation, legacy schema-v1 adaptation, serialization helpers, and compiler-input shell.
- `R/render_scene.R`, `R/add_object.R`, `R/ray_scene.R`: legacy-compatible entry points that detect schema-v2 scenes and expose the spectral validation shell without changing `rgb_legacy` rendering.
- `inst/extdata/spectral/`: generated named spectral assets, per-spectrum metadata, and the pbrt-layout sRGB RGB-to-spectrum coefficient table used by the spectral base layer.
- `src/render_scene_rcpp.cpp`: Rcpp still-render entry point; it may own R-facing result assembly and preview overlays.
- `src/render_animation_rcpp.cpp`: Rcpp animation entry point; it may own frame iteration and R `post_process_frame` calls.
- `src/core/`: existing RGB integrators, cameras, scene building, preview display, and low-level renderer code.

## Rule

Do not add new top-level renderer setup duplication to the Rcpp entry points. New shared setup belongs under `src/render/`, and entry points should pass explicit inputs into `RenderSession` or a narrower helper there.

The PR 1 render-session layer is intentionally RGB-only. PR 2 through PR 11 spectral code is
scaffolding and must remain unused by legacy transport except for isolated conversion, packet,
asset-loading, RGB reconstruction, descriptor validation, texture evaluation, Film/sensor conversion,
spectral camera generation, spectral scene/shape adaptation, explicit spectral light construction,
and allocation tests.
