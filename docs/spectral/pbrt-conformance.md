# pbrt Conformance Matrix

This table is the running record of which pbrt-v4 contracts have been ported, which rayrender files implement them, and which deviations are intentional. The pinned pbrt reference is recorded in `docs/spectral/versions.md`.

| Subsystem | pbrt file/function | rayrender implementation | Status | Tests | Deviation |
|---|---|---|---|---|---|
| Project pin | Git commit `8c19f304558fd7681e2fef2c395a689d0106fb05` | `docs/spectral/versions.md` | Pinned in PR 0 | Manual pin check | None |
| Architecture decisions | pbrt-v4 renderer architecture | `docs/spectral/adr/0001-spectral-renderer.md` | Recorded in PR 0 | Documentation review | rayrender keeps legacy RGB renderer during migration |
| Render session and source layout | pbrt-v4 renderer setup concepts | `src/render/render_session.*`, `docs/spectral/source-layout.md` | Extracted in PR 1 | Local install, deterministic source tests, PR 0 legacy baseline compare | Legacy RGB setup only; spectral state is deferred |
| Base utility types | pbrt utility color types, typed dispatch handles, scratch allocation, BxDF flags | `src/base/*`, `tools/spectral-tests/pr2-base-tests.cpp` | Foundation in PR 2 | `Rscript tools/spectral-tests/run-pr2-base-tests.R`, local install, PR 0 legacy baseline compare | Isolated from legacy transport; `RGBBSDFSample` is a temporary convention until `SampledSpectrum` lands |
| Wavelength packet | `src/pbrt/util/sampling.h::SampleVisibleWavelengths`, `VisibleWavelengthsPDF`; `src/pbrt/util/spectrum.h::SampledWavelengths` | `src/base/sampled_spectrum.h`, `tools/spectral-tests/pr3-spectrum-tests.cpp` | Foundation in PR 3 | `Rscript tools/spectral-tests/run-pr3-spectrum-tests.R`, PR 0 legacy baseline compare | Isolated from legacy transport; Film sampling integration deferred |
| Sampled spectrum | `src/pbrt/util/spectrum.h::SampledSpectrum` arithmetic and diagnostics | `src/base/sampled_spectrum.h`, `src/base/base_test.cpp` | Foundation in PR 3 | `Rscript tools/spectral-tests/run-pr3-spectrum-tests.R`, local install | XYZ/RGB sensor conversion deferred until Film and named spectra |
| Spectrum hierarchy | `Spectrum`, RGB spectrum classes, spectral data | Not implemented | Not started | Pending PRs 4-5 | None recorded |
| RGB color space and reconstruction | `RGBColorSpace`, `RGBToSpectrumTable` | Not implemented | Not started | Pending PR 5 | None recorded |
| Texture split | `FloatTexture`, `SpectrumTexture`, texture evaluators | Not implemented | Not started | Pending PR 7 | None recorded |
| Film and sensor | `Film`, `PixelSensor` | Not implemented | Not started | Pending PR 8 | None recorded |
| Camera interface | `Camera::GenerateRay()` with mutable wavelengths | Not implemented | Not started | Pending PR 9 | None recorded |
| Shape and primitive split | `Shape`, `Primitive`, `SurfaceInteraction` | Not implemented | Not started | Pending PR 10 | CSG will be a rayrender-specific Shape extension |
| Lights | `Light`, `LightSampler`, infinite lights | Not implemented | Not started | Pending PR 11 | None recorded |
| BxDF and BSDF | `BxDF`, `BSDF`, `BSDFSample`, Fresnel, microfacet | Not implemented | Not started | Pending PR 12 | None recorded |
| Materials | `Material::GetBSDF()` closure model | Not implemented | Not started | Pending PR 13 | Legacy materials remain isolated to `rgb_legacy` |
| Integrators | `RandomWalkIntegrator`, `PathIntegrator`, `VolPathIntegrator` | Not implemented | Not started | Pending PRs 14-15 and 20 | None recorded |
| Nested dielectric regions | rayrender extension over pbrt dielectric boundaries | Not implemented | Not started | Pending PRs 17-18 | Lower numeric priority wins; implemented as explicit extension |
| Media | `Medium`, `PhaseFunction`, `VolPathIntegrator` | Not implemented | Not started | Pending PR 20 | None recorded |
| Realistic camera dispersion | pbrt realistic camera and measured sensors | Not implemented | Not started | Pending PR 21 | rayrender lens-file differences must be documented |
| Spectral assets | pbrt data/tables and generated assets | `docs/spectral/provenance.md`, `docs/spectral/assets-manifest.csv` | Inventory scaffolded in PR 0 | `tools/spectral-tests/check-spectral-assets.R` | No spectral assets copied yet |
| Legacy baselines | Current rayrender RGB renderer | `docs/spectral/baselines/` | Captured in PR 0 | `tools/spectral-tests/capture-pr0-baselines.R` | None; no renderer code changed |

## Status Values

- `Pinned`: reference or decision is recorded and immutable until a later ADR changes it.
- `Not started`: implementation has not begun.
- `In progress`: implementation exists but the PR gate has not passed.
- `Conformant`: implementation and tests match the pinned pbrt behavior.
- `Extension`: intentionally differs from pbrt and is documented by ADR or plan section.
- `Blocked`: implementation cannot proceed without a documented decision.
