# Renderer Source Layout

PR 1 introduces a shared render-session layer without changing RGB transport. PR 2 through PR 5 add
isolated base types, packaged spectral assets, and pbrt-style sRGB spectrum reconstruction. PR 6 adds
the R schema-v2 descriptor and compiler-input validation shell. PR 7 adds the first C++ spectral
texture layer. PR 8 adds the spectral PixelSensor and Film accumulation layer. PR 9 adds the
wavelength-aware spectral camera interface. PR 10 adds spectral shape, primitive, aggregate, and
scene scaffolding. PR 11 adds explicit spectral lights, infinite lights, light registries, and
light samplers. PR 12 adds the isolated spectral BxDF/BSDF scattering foundation. PR 13 adds
the isolated pbrt-style spectral Material closure layer for diffuse and interface/null materials.
PR 14 adds the isolated spectral RandomWalk vertical slice for diffuse surfaces and explicit
area/infinite emitters. PR 15 adds the isolated spectral PathIntegrator vertical slice with
direct-light sampling, BSDF-sampled emitter-hit MIS, and Russian roulette scaffolding for vacuum
diffuse scenes. PR 16 adds spectral conductor BxDF/material support, measured eta/k and
reflectance-compatibility descriptor paths, and conductor coverage in the isolated path tests.
PR 17 adds constant-IOR dielectric region state, smooth dielectric BSDF/material support, and
null/skipped dielectric boundary traversal for the isolated RandomWalk and Path integrators.
PR 18 adds proof-aware spectral eta resolution, rough dielectric and thin-dielectric closure
support, and pbrt-style wavelength termination diagnostics for dispersive dielectric paths.

## Current Layout

- `src/render/`: renderer setup orchestration shared by still and animation entry points.
- `src/render/render_session.*`: legacy RGB setup bridge for scene compilation, camera construction, environment light setup, output buffers, and explicit frame inputs.
- `src/render/spectral_film.*`: PR 8 PixelSensor, Film, FilmTile, output-color conversion, visible-surface hooks, and deterministic test accumulation.
- `src/render/spectral_camera.*`: PR 9 spectral-only camera sample, ray, handle, perspective/orthographic/thin-lens generation, ray differentials, initial-medium recording, dielectric-region initialization hooks, and Film wavelength-sampling bridge.
- `src/render/spectral_scene.*`: PR 10 spectral Interaction, SurfaceInteraction, Shape callbacks, native sphere and SDF CSG adapters, primitive bindings, transformed primitives, aggregates, Scene ownership, ray spawning, capability validation, and legacy hitable wrapping; PR 17 adds multi-boundary dielectric region attachments, binding lookup for camera-origin containment, and coherent CSG region-boundary capability support.
- `src/render/spectral_light.*`: PR 11 spectral Light flags/types, point/spot/distant/diffuse-area lights, uniform/image infinite lights, visibility endpoints, scalar light power estimates, UniformLightSampler/PowerLightSampler, and schema/compiler adapter hooks for area lights, free lights, environments, and legacy emissive materials.
- `src/render/spectral_bsdf.*`: PR 12 spectral `BSDFSample`, local scattering-coordinate helpers, Fresnel dielectric/conductor functions, `TrowbridgeReitzDistribution`, `DiffuseBxDF`, `NullBxDF`, non-owning `BxDF` dispatch, and the pbrt-style `BSDF` frame wrapper; PR 16 `ConductorBxDF`, conductor `rho()` estimates, smooth specular and rough microfacet conductor sampling/PDF behavior; PR 17 smooth constant-ratio `DielectricBxDF` reflection/transmission/TIR sampling; PR 18 full rough `DielectricBxDF` f/Sample_f/PDF behavior and `ThinDielectricBxDF`.
- `src/render/spectral_dielectric.*`: PR 17 `DielectricRegionTable`, `DielectricPathState`, side-aware region membership analysis, checked transition tokens, exterior region handling, priority-skipped and index-matched null traversal classification, and analytic/CSG point-containment initialization; PR 18 sampled, Cauchy, Sellmeier, and named-glass eta handles plus proof-based effective-ratio constancy detection.
- `src/render/spectral_integrator.*`: PR 14 deterministic sample streams, `RandomWalkIntegrator`, per-worker scratch state, Film render loop, vacuum-only/delta-light validation, and radiance diagnostics; PR 15 `PathIntegrator`, direct-light sampling through `UniformLightSampler`, BSDF-sampled area/infinite emitter-hit MIS, max-depth ordering, Russian roulette accounting, and Path Film render loop; PR 17 copied visibility region traversal, dielectric diagnostics, null/skipped boundary traversal, and commit-on-transmission state updates; PR 18 wavelength-termination counters for dispersive dielectric and thin-dielectric material evaluation.
- `src/base/`: explicit RGB/XYZ/color-space types, sRGB transfer functions, RGB-to-spectrum table loading, RGB spectrum wrappers, `SampledSpectrum`, `SampledWavelengths`, spectrum representations and named-spectrum registry loading, PR 18 Cauchy/Sellmeier IOR spectra, tagged dispatch handles, `ScratchBuffer`, BxDF flags, and optional sample-result conventions for future pbrt-style transport.
- `src/materials/spectral_texture.*`: PR 7 spectral texture descriptors, evaluators, semantic image cache keys, image loading adapters, and hit-record derivative adapters.
- `src/materials/spectral_material.*`: PR 13 spectral `MaterialEvalContext`, tagged Material dispatch, diffuse and interface/null material parameter objects, alpha and bump hooks, scratch-allocated `BSDF` closure construction, and a spectral material handle table; PR 16 `ConductorMaterial` with eta/k spectra, sampled reflectance compatibility, anisotropic roughness textures, remapping, and packet-preserving closure construction; PR 17 constant-ratio `DielectricMaterial` fed by `ResolvedDielectricInterface`; PR 18 rough/dispersive `DielectricMaterial` wavelength termination and `ThinDielectricMaterial`.
- `R/schema_v2_descriptors.R`: PR 6 schema-v2 spectral descriptors, scene decorators, validation, legacy schema-v1 adaptation, serialization helpers, and compiler-input shell; PR 17 rejects nonzero dielectric absorption until media support and maps legacy dielectric refraction/priority into optical-region metadata.
- `R/render_scene.R`, `R/add_object.R`, `R/ray_scene.R`: legacy-compatible entry points that detect schema-v2 scenes and expose the spectral validation shell without changing `rgb_legacy` rendering.
- `inst/extdata/spectral/`: generated named spectral assets, per-spectrum metadata, and the pbrt-layout sRGB RGB-to-spectrum coefficient table used by the spectral base layer.
- `src/render_scene_rcpp.cpp`: Rcpp still-render entry point; it may own R-facing result assembly and preview overlays.
- `src/render_animation_rcpp.cpp`: Rcpp animation entry point; it may own frame iteration and R `post_process_frame` calls.
- `src/core/`: existing RGB integrators, cameras, scene building, preview display, and low-level renderer code.

## Rule

Do not add new top-level renderer setup duplication to the Rcpp entry points. New shared setup belongs under `src/render/`, and entry points should pass explicit inputs into `RenderSession` or a narrower helper there.

The PR 1 render-session layer is intentionally RGB-only. PR 2 through PR 18 spectral code is
scaffolding and must remain unused by legacy transport except for isolated conversion, packet,
asset-loading, RGB reconstruction, descriptor validation, texture evaluation, Film/sensor conversion,
spectral camera generation, spectral scene/shape adaptation, explicit spectral light construction,
BSDF/BxDF scattering tests, material closure construction, alpha/bump tests, allocation tests,
RandomWalk transport tests, PathIntegrator estimator tests, conductor material tests, and Film
render-loop tests, plus constant-IOR dielectric region-state, smooth dielectric path tests, rough
dielectric scattering tests, and dispersive/thin-dielectric wavelength-termination tests.
