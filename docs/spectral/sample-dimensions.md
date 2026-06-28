# Spectral Sample Dimensions

This file records stable sampler-dimension assignments for the spectral renderer.
Any intentional reorder must update deterministic references and this document in
the same change.

| Dimension | Use |
|---|---|
| 0 | Pixel filter sample, 2D |
| 1 | Wavelength sample, 1D |
| 2 | Lens sample, 2D |
| 3 | Time sample, 1D |
| 4 | Light selection, 1D |
| 5 | Light position/direction, 2D |
| 6 | BSDF component selection, 1D |
| 7 | BSDF direction sample, 2D |
| 8 | Alpha test, 1D |
| 9 | Russian roulette, 1D |
| 10 | Medium distance/collision sample, 1D |
| 11 | Phase-function direction sample, 2D |
| 12 | Layered-BxDF interface selection, 1D |
| 13 | Layered-BxDF direction sample, 2D |

Rules:

- Adaptive sampling must not reuse wavelength packet components as color
  channels.
- Visibility rays use the same alpha policy with their own local dimension
  stream.
- Medium and phase-function dimensions remain reserved even while `VolPath` is
  deferred.
- Additional feature dimensions must be appended unless a reference update
  intentionally changes the sequence.
