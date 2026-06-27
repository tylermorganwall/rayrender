# Spectral Renderer Version Pins

This file records the immutable project baselines for the spectral renderer work.

## rayrender

- Package version: `0.41.3`
- Package date: `2026-06-15`
- Git commit: `0ab9b6a2de4b075f91daefee6d38487782cabd9e`
- Commit date: `2026-06-15T21:24:06-04:00`
- Commit subject: `rayrender v0.41.3: Update SIMD config for portability on windows (CRAN)`
- Source checksum: `53ab5e593dd457d3b13297e2e4a0559f74a5d1fdb279b4a3ef47dc43f2fa202f`
- Checksum command: `git archive --format=tar HEAD | shasum -a 256`

The source checksum is for the committed Git tree at the pin above. It intentionally excludes uncommitted workspace changes, build products, and ignored files.

## pbrt-v4

- Local reference path: `mmp/pbrt-v4`
- Remote: `https://github.com/mmp/pbrt-v4.git`
- Git commit: `8c19f304558fd7681e2fef2c395a689d0106fb05`
- Commit date: `2025-12-08T14:30:32-08:00`
- Commit subject: `Merge pull request #521 from gonsolo/noseven`
- License file: `mmp/pbrt-v4/LICENSE.txt`
- License: Apache-2.0

## Pin Policy

The pbrt commit above is the normative implementation reference for the spectral renderer. Do not update it silently. Any change to the pbrt pin requires:

1. an ADR describing why the pin changed;
2. updates to `docs/spectral/pbrt-conformance.md`;
3. regenerated reference fixtures where the executable behavior changed;
4. maintainer approval before the new pin is used for implementation work.

