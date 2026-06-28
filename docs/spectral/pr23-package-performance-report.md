# PR23 Package and Performance Report

PR23 adds DCI-P3, Rec.2020, and ACES2065-1 color-space definitions, name
normalization, and table filename/loader policy without changing spectral
estimator behavior. Only the sRGB RGB-to-spectrum table remains packaged.

## Asset Footprint

Measured on the PR23 working tree:

- `inst/extdata/spectral`: 9.3 MiB total.
- `named-spectra-v1.tsv`: 265,629 bytes.
- `named-spectra-v1-metadata.tsv`: 23,610 bytes.
- `rgb-to-spectrum-srgb-v1.bin`: 9,437,528 bytes.
- source tarball built from PR23 tree: 9,880,710 bytes.
- installed package from `R CMD check`: 17.4 MiB total, including 9.5 MiB
  `extdata` and 6.6 MiB `libs`.
- local pbrt generated source tables:
  - `rgbspectrum_srgb.cpp`: 37,528,351 bytes;
  - `rgbspectrum_dci_p3.cpp`: 37,407,551 bytes;
  - `rgbspectrum_rec2020.cpp`: 37,235,278 bytes;
  - `rgbspectrum_aces.cpp`: 36,930,965 bytes.

The compact packaged sRGB table is 9.0 MiB. Adding DCI-P3, Rec.2020, and
ACES2065-1 at the same binary layout would add roughly 27 MiB of package data.
PR23 therefore treats those tables as optional assets until package-size policy
explicitly permits bundling them. The loader has deterministic filenames and
clear missing-table diagnostics for those future assets:

- `rgb-to-spectrum-dci-p3-v1.bin`;
- `rgb-to-spectrum-rec2020-v1.bin`;
- `rgb-to-spectrum-aces2065-1-v1.bin`.

## Implementation Notes

- `RGBColorSpace` now defines canonical sRGB, DCI-P3, Rec.2020, and
  ACES2065-1 primaries and white points.
- The R schema-v2 descriptor layer accepts only canonical color-space names:
  `sRGB`, `DCI-P3`, `Rec.2020`, and `ACES2065-1`.
- RGB-to-spectrum table cache keys now include both path and expected
  color-space id, so a cached sRGB table cannot satisfy a DCI-P3 request.
- RGB table load counters record lookups, cache hits, and file loads.
- Non-sRGB spectral reconstruction throws a color-space-specific missing asset
  error until the matching binary table is packaged and manifest-listed.

## Performance Scope

The spectral renderer is still isolated from legacy transport and the direct
schema-v2 renderer bridge is not enabled. PR23 therefore does not run rendered
single-thread or multi-thread spectral benchmarks against legacy scenes. This is
a measured exception to Section 16: the package slice can validate table loading,
schema validation, asset lookup, and sanitizer behavior, but cannot yet measure
end-to-end Path performance or estimator bias for pbrt-equivalent scenes.

No spectral clamping, packet resampling, or wavelength-dependent light-selection
shortcut was introduced in PR23.

## Gate Evidence

Commands used for the PR23 package slice:

- `air format R/schema_v2_descriptors.R tests/testthat/test-schema-v2-descriptors.R`
- `Rscript tools/spectral-tests/run-pr5-rgb-spectrum-tests.R`
- `Rscript -e "devtools::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-schema-v2-descriptors.R')"`
- `Rscript tools/spectral-tests/check-spectral-assets.R --manifest docs/spectral/assets-manifest.csv --build-source`
- `R CMD INSTALL .`
- `R CMD build --no-build-vignettes --no-manual /Users/tyler/Desktop/R/rayrender`
- `R CMD check --no-manual --no-vignettes /private/tmp/rayrender_0.41.3.tar.gz`

The PR5 RGB spectrum command compiles and runs normal, installed-asset, and
sanitizer variants. The schema-v2 descriptor test passed with 144 assertions.
The built-tarball package check completed with 3 warnings and 1 note, all in
pre-existing Rd documentation files (`as_scene_compiler_input.Rd`,
`validate_ray_scene_v2.Rd`, and `schema_v2_materials.Rd`).
