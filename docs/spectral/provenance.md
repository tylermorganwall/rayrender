# Spectral Data and Source Provenance

This file tracks copied or generated spectral assets and source-derived implementation material. PR 0 adds the inventory and notice policy only; it does not copy pbrt source, pbrt data tables, RGB-to-spectrum tables, or measured spectral datasets into rayrender.

## Current Inventory

| Asset or source | Location in rayrender | Origin | License | Status | Required notice |
|---|---|---|---|---|---|
| pbrt-v4 reference source | `mmp/pbrt-v4` | `https://github.com/mmp/pbrt-v4.git` at `8c19f304558fd7681e2fef2c395a689d0106fb05` | Apache-2.0 | Reference only; not copied into package | Preserve copyright and Apache-2.0 notices before copying or adapting source |
| Spectral implementation plan | `docs/spectral/rayrender_spectral_rendering_codex_plan.md` | Project-local plan document | Project documentation | Present before PR 0 implementation | None beyond repository license policy |
| Spectral packaged assets | None | None | None | Not started | Future assets must be added to `docs/spectral/assets-manifest.csv` |

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

