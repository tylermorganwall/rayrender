# ADR 0001: pbrt-Style Spectral Renderer Architecture

- Status: Accepted for implementation plan
- Date: 2026-06-26
- Scope: PRs 0-24 of the spectral renderer migration

## Context

rayrender 0.41.3 is an RGB path tracer with materials, emission, hitables, media-like constant volumes, and dielectric priority behavior coupled through legacy interfaces. The spectral migration needs pbrt-v4-style transport while preserving existing RGB behavior until the new path reaches parity.

## Decision

The new renderer will be a separate `render_mode = "spectral"` path modeled on the pinned pbrt-v4 commit in `docs/spectral/versions.md`.

Core decisions:

1. Spectral transport uses pbrt's four-sample wavelength packets, `SampledSpectrum`, role-aware RGB reconstruction, Film, and PixelSensor model.
2. The current RGB renderer remains available as `render_mode = "rgb_legacy"` throughout the migration.
3. Materials become immutable parameter objects whose per-hit evaluated closures are scratch-allocated BxDFs wrapped by a BSDF.
4. Emission is represented by explicit Light objects, including Primitive-bound area lights and scene-level infinite lights.
5. Textures split into scalar and spectral texture interfaces with explicit color encoding and RGB reconstruction role.
6. Scene compilation moves toward schema-v2 descriptors with named, validated fields and centralized schema-v1 adaptation.
7. Camera, Film, Material, Light, Medium, and integrator APIs receive wavelength packets explicitly. Geometric Ray state does not own wavelengths.
8. Participating media use pbrt-style Medium, PhaseFunction, and VolPathIntegrator contracts, not the legacy stochastic surface model.
9. rayrender's lower-number-wins overlapping dielectric priority behavior is retained as an explicit extension, implemented as path-local region membership state rather than mutable dielectric pointers on rays.
10. CSG remains supported as a rayrender-specific Shape extension with explicit capability flags and validation limits.

## Consequences

Positive consequences:

- pbrt conventions for BSDFs, PDFs, MIS, wavelength sampling, Film integration, and media remain testable in isolation.
- Existing scenes can continue to render through `rgb_legacy` while spectral behavior is developed.
- Region and medium state become value-semantic path state, which avoids aliasing between camera paths and visibility paths.
- Public R APIs can expose physical units and roles without forcing existing constructors to change immediately.

Costs and constraints:

- The migration requires new renderer infrastructure instead of incrementally converting all `point3f` radiometry to spectra.
- Legacy material and light conveniences need centralized adapters and compatibility diagnostics in spectral mode.
- pbrt source/data/table licensing and provenance must be reviewed before copying or generating derived assets.
- Feature parity is staged across many PRs, so the spectral renderer remains experimental until the final gates pass.

## Non-Decisions

- The spectral renderer does not become the default in PR 0.
- No pbrt source, RGB-to-spectrum table, measured spectrum, or generated binary asset is copied in PR 0.
- Fluorescence, polarization, diffraction, wavelength redistribution, GPU wavefront rendering, and bidirectional methods remain out of scope for the initial completion target.

