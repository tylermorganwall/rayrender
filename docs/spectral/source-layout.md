# Renderer Source Layout

PR 1 introduces a shared render-session layer without changing RGB transport.

## Current Layout

- `src/render/`: renderer setup orchestration shared by still and animation entry points.
- `src/render/render_session.*`: legacy RGB setup bridge for scene compilation, camera construction, environment light setup, output buffers, and explicit frame inputs.
- `src/render_scene_rcpp.cpp`: Rcpp still-render entry point; it may own R-facing result assembly and preview overlays.
- `src/render_animation_rcpp.cpp`: Rcpp animation entry point; it may own frame iteration and R `post_process_frame` calls.
- `src/core/`: existing RGB integrators, cameras, scene building, preview display, and low-level renderer code.

## Rule

Do not add new top-level renderer setup duplication to the Rcpp entry points. New shared setup belongs under `src/render/`, and entry points should pass explicit inputs into `RenderSession` or a narrower helper there.

The PR 1 layer is intentionally RGB-only. Spectral types, wavelength state, and pbrt spectral assets are introduced in later PRs.
