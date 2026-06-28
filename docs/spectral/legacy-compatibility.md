# Legacy Compatibility Guide

This guide records PR24 compatibility policy for the spectral renderer rollout.

## Current Default

`rgb_legacy` remains the default render mode. Existing scenes should continue to
render through the legacy RGB renderer unless a user explicitly requests:

```r
render_mode = "spectral"
```

In PR24, spectral mode returns compiler input only when `return_result = TRUE`.
It does not silently fall back to RGB rendering.

## Legacy Conversion Policy

Legacy constructors are adapted once through `legacy_scene_to_schema_v2()`.
Ambiguous mappings are reported in a scene conversion report and warning text:

- `metal()` maps to a compatibility conductor from legacy RGB reflectance and
  fuzz. Prefer measured `conductor(eta = ..., k = ...)`.
- `glossy()` maps to a spectral coated diffuse material using legacy RGB base
  color, roughness, and normal-incidence reflectance.
- `light()` as a material is a legacy area-light shorthand. Prefer
  `with_light(area_light(...))`.
- `dielectric()` maps to a dielectric interface plus optical-region metadata.
  Legacy absorption remains deferred until participating media are fully wired.

The converter must not treat emitted RGB as albedo or albedo RGB as emission.

## Expected Appearance Changes

Spectral mode is not expected to match RGB legacy images exactly. Differences
can come from:

- RGB-to-spectrum reconstruction roles: albedo, illuminant, and unbounded RGB
  are different physical quantities.
- Spectral Fresnel and measured optical constants.
- CIE or RGB sensor integration and output color-space conversion.
- Corrected BSDF cosine/PDF conventions.
- Explicit light sampling and infinite-light treatment.
- Absorption and dispersion once participating media and spectral glass are
  fully routed through rendering.

Do not add arbitrary correction factors to force spectral output toward legacy
RGB output.

## Warning Schedule

- Current release: spectral mode is opt-in and compiler-input only.
- Next stable release: keep `rgb_legacy` as the default while collecting
  spectral acceptance data.
- Default-switch candidate: requires accepted PR24 gates, release notes, and
  maintainer approval.
- Deprecation: open a separate future plan after at least one stable spectral
  release.

`render_mode = "rgb_legacy"` remains the compatibility escape hatch for at
least one stable release after any future default switch.
