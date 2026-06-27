# rayrender 0.41.3 Fully Spectral Renderer

## Codex implementation specification: pbrt-v4 conformance revision 3 3

**Repository baseline:** attached `rayrender` 0.41.3 source package, dated 2026-06-15  
**Primary behavioral reference:** a pinned commit of `mmp/pbrt-v4` (available in this directory at mmp/pbrt-v4), plus the online fourth edition of *Physically Based Rendering*

This document is written for Codex. Implement exactly one numbered pull request at a time. Do not skip prerequisite PRs, combine adjacent PRs, or opportunistically convert unrelated legacy code.

The objective is not merely to add wavelength-valued colors to the existing integrator. The objective is to add a pbrt-v4-style rendering core to rayrender: sampled spectra, pbrt wavelength sampling, role-aware RGB reconstruction, pbrt-style texture evaluation, per-hit material closures, BxDF/BSDF contracts, explicit lights and light samplers, Film and PixelSensor, path-integrator MIS, participating media, and wavelength-dependent cameras. The existing RGB renderer remains available until the new renderer reaches the acceptance criteria in this document.

Rayrender's nested-dielectric priority behavior is retained as an explicit extension. It is not implemented by keeping mutable `dielectric*` pointers on spectral rays. Instead, a path-local dielectric-region membership state resolves the effective optical media on the two sides of each boundary and supplies that resolved interface to a pbrt-style `DielectricBxDF` closure.

This plan adds normative public API instructions, `ray_scene` and `ray_material` schema changes, the functional object-decorator model for area lights and region boundaries, multi-primitive optical-region assemblies, and CSG-specific support and limitations. These additions are part of the implementation requirements, not optional user-interface polish.

## Contents

1. Codex usage and non-negotiable decisions
2. Definition of fully spectral and pbrt wavelength model
3. Current and target architecture
4. Core invariants and normative interfaces
5. R API and schema-v2 scene descriptors
6. Material closure mapping and PathIntegrator behavior
7. Nested and overlapping dielectric-region specification
8. CSG and implicit geometry integration
9. PR 0 through PR 24 implementation sequence
10. Unit, statistical, render, and pbrt comparison tests
11. Feature integration, performance, packaging, and compatibility
12. Risks, Definition of Done, Codex execution protocol, and pbrt reference map

---

# 1. How Codex must use this document

1. Note that this file is in the repository as `docs/spectral/rayrender_spectral_rendering_codex_plan.md`.
2. Read the repository `AGENT.md` before every PR and use `tools/codex/install-local.sh` for local package installation.
3. In PR 0, pin the exact rayrender and pbrt-v4 Git commits used for the project. Do not silently track later pbrt changes after that point.
4. Treat sections marked **MUST** as acceptance requirements. A deviation requires an ADR entry, a focused test, and explicit approval.
5. Implement one PR, run its gate, report results, and stop. The next PR starts only after the current one is accepted.
6. Preserve `render_mode = "rgb_legacy"` output throughout the migration unless a separately reviewed legacy bug fix is requested.
7. Use `=` for assignment in all R code. Use ASCII only in source-code text.
8. Do not copy pbrt source or generated tables without preserving provenance and completing the repository's license and notice work.

Normative terms in this plan use the meanings in RFC-style specifications:

- **MUST**: required for correctness or acceptance.
- **SHOULD**: expected unless a documented reason applies.
- **MAY**: optional.

---

# 2. Non-negotiable architectural decisions

## 2.1 pbrt-v4 is the normative model for the new renderer

The new spectral renderer MUST follow the pinned pbrt-v4 implementation for:

- the visible wavelength range and four-sample correlated wavelength packet;
- `SampledSpectrum` and `SampledWavelengths` semantics;
- secondary-wavelength termination;
- `Spectrum`, RGB spectrum roles, RGB color spaces, and RGB-to-spectrum table lookup;
- scalar and spectral texture separation;
- `Material::GetBSDF()` and per-hit concrete BxDF evaluation;
- BxDF flags, `BSDFSample`, `TransportMode`, BSDF coordinate conventions, and PDF measures;
- Fresnel, Trowbridge-Reitz, conductor, dielectric, thin dielectric, coated, and hair models where ported;
- explicit `Light`, `LightSampler`, area-light, and infinite-light behavior;
- Film wavelength sampling and PixelSensor integration;
- camera APIs that may mutate sampled wavelengths;
- PathIntegrator state, emitter-hit MIS, direct-light MIS, throughput updates, `etaScale`, regularization, and Russian roulette;
- Medium and VolPathIntegrator behavior when the volume stages are reached.

The goal is behavioral and mathematical conformity, not source-layout identity.

## 2.2 Allowed implementation-level deviations

The following deviations from pbrt are permitted and expected:

1. **CPU dispatch:** rayrender MAY use `std::variant` or an explicit enum-plus-pointer tagged handle rather than pbrt's pointer-bit `TaggedPointer`. The public contracts and evaluation order MUST remain pbrt-like.
2. **Memory allocation:** rayrender MUST use a per-worker scratch arena for per-hit closures, but its allocator implementation need not be pbrt's allocator.
3. **R scene descriptions:** rayrender uses versioned R descriptors instead of pbrt scene files.
4. **Package assets:** RGB-to-spectrum tables MAY be distributed as checked binary package assets rather than generated C++ arrays.
5. **Nested dielectric priority:** rayrender retains lower-number-wins overlapping dielectric regions. This is an explicit extension described in Section 12.
6. **Legacy coexistence:** the existing RGB renderer remains compiled and selectable during migration.
7. **Scheduling and preview:** rayrender retains its RcppThread, adaptive rendering, and preview infrastructure, after adapting them to Film.
8. **CSG support:** rayrender retains CSG as a rayrender-specific Shape extension. CSG must satisfy the same Shape, Primitive, Material, Light, Medium, and region-boundary contracts as other geometry, but it is not expected to exist in pbrt itself.

Every other intentional difference MUST be recorded in `docs/spectral/pbrt-conformance.md`.

## 2.3 Do not implement a generic closure graph first

In this plan, a material "closure" means the evaluated per-hit scattering state represented by a concrete BxDF allocated in scratch storage and wrapped in a BSDF. This is the pbrt-v4 model.

Do not begin with an unrestricted recursive closure DAG or a heap-allocated list of arbitrary lobes. Compound behavior SHOULD be implemented with explicit pbrt-style BxDFs such as:

- `DielectricBxDF`;
- `CoatedDiffuseBxDF`;
- `CoatedConductorBxDF`;
- `ThinDielectricBxDF`;
- `HairBxDF`;
- later, `LayeredBxDF` if required.

This keeps `f()`, `Sample_f()`, and `PDF()` mutually consistent and makes material sampling testable.

## 2.4 Keep the legacy renderer separate until parity

Do not incrementally change all existing `point3f` radiometry to spectra. Add a new renderer behind:

```r
render_mode = c("rgb_legacy", "spectral")
```

During migration:

- legacy material classes, `scatter_record`, `pdf`, `color.cpp`, and the current dielectric pointer stack remain available only to `rgb_legacy`;
- new spectral code lives in the `rayrender` namespace and uses new interfaces;
- scene schema v1 is adapted into schema v2 for spectral compilation;
- legacy seeded regression outputs remain unchanged.

## 2.5 Wavelengths do not belong in geometric ray intersection state

The new `Ray` MAY carry a non-owning current `Medium` handle, as pbrt does, but it MUST NOT carry a `SampledWavelengths` packet. The integrator owns wavelengths and passes them to camera, material, texture, light, medium, and film evaluations.

The nested-region membership set is also path state, not geometric intersection state.

## 2.6 Emission is a light property, not a material side effect

The spectral path MUST not call `material::emitted()`. A primitive may reference both a surface material and an area light. Infinite lights are evaluated on ray misses. This is required for pbrt-style light sampling and MIS.

## 2.7 The output remains tristimulus

"Fully spectral" means that transport quantities are spectral. Film normally stores three sensor channels, not hundreds of wavelength bins. `PixelSensor` integrates radiance samples against sensor response curves and divides by wavelength PDFs. Output RGB conversion occurs after transport.

---

# 3. Definition of fully spectral for this project

All quantities below MUST be spectral in `render_mode = "spectral"`:

- path radiance and throughput;
- emitted radiance and light intensity distributions;
- diffuse, conductor, dielectric, coated, hair, and phase-function weights;
- conductor `eta` and `k`;
- dielectric index of refraction and absorption;
- medium absorption, scattering, emission, and transmittance;
- RGB image textures after role-specific spectral reconstruction;
- environment maps after illuminant reconstruction;
- realistic-camera lens IOR where lens data permit dispersion;
- sensor response and final spectral-to-tristimulus integration.

The following are outside the initial completion target:

- fluorescence or phosphorescence;
- polarization;
- diffraction and wave optics;
- wavelength redistribution;
- splitting one wavelength packet into independently traced spectral rays;
- bidirectional path tracing, MLT, or SPPM;
- GPU wavefront rendering.

These exclusions do not justify shortcuts in the surface and volume path tracers.

---

# 4. Correct pbrt-v4 wavelength model

Use a compile-time sample count of four:

```cpp
inline constexpr int NSpectrumSamples = 4;
inline constexpr Float LambdaMin = 360.f;
inline constexpr Float LambdaMax = 830.f;
```

For each camera sample:

1. Draw one scalar wavelength sample `u`.
2. Generate four correlated visible wavelengths using shifted copies of `u`.
3. Store each wavelength and its density in `SampledWavelengths`.
4. Trace one geometric path while evaluating `SampledSpectrum` values at those wavelengths.
5. When a camera or scattering closure must select a wavelength-dependent direction, call `TerminateSecondary()` before selecting the direction.
6. Continue the same path with only `lambda[0]` active.

The implementation MUST match pbrt's visible-wavelength inverse CDF and PDF at the pinned commit. Do not substitute uniform wavelength sampling in production merely because it is simpler.

`TerminateSecondary()` MUST be idempotent and equivalent to:

```cpp
void SampledWavelengths::TerminateSecondary() {
    if (SecondaryTerminated()) {
        return;
    }
    for (int i = 1; i < NSpectrumSamples; ++i) {
        pdf[i] = 0;
    }
    pdf[0] /= NSpectrumSamples;
}
```

Wavelength-PDF division MUST occur exactly once, in PixelSensor/Film conversion. Do not divide path throughput by the wavelength PDF at every bounce.

---

# 5. Current rayrender architecture and required response

The following audit is based on the attached 0.41.3 source.

| Current source | Current behavior | Required spectral response |
|---|---|---|
| `src/core/ray.h` | `Ray` stores geometry plus a raw pointer to `std::vector<dielectric*>`. | Retain for legacy only. New path uses a value-semantic dielectric-region state; `Ray` may gain a pbrt-like `Medium` handle. |
| `src/materials/material.h` | One virtual `material` combines scattering, `f`, emission, and AOV albedo. | Add immutable pbrt-style Material handles whose `GetBxDF()` creates per-hit closures. Emission moves to Light. |
| `scatter_record` | Owns a raw `pdf*`, stores RGB attenuation, and switches between specular and sampled behavior. | Do not reuse in spectral mode. Replace with `BSDFSample` and BSDF sampling. |
| `src/core/color.cpp` | Three integrators use RGB `point3f`; BSDF/cosine conventions are inconsistent. | Add new pbrt-style integrators. BSDF `f` excludes cosine; the integrator applies `AbsDot()` once. |
| `src/core/integrator.cpp` | Allocates one mutable dielectric pointer vector per worker and writes directly to `RayMatrix`. | Add per-worker sampler clone, scratch buffer, Film tile, and path-local dielectric state. |
| `src/hitables/hitable.h` | Intersection, material binding, emission geometry, and light sampling are combined. | Add pbrt-style Shape/Primitive/Light separation, initially through adapters. |
| `src/hitables/infinite_area_light.*` | Environment is an emissive intersectable sphere with RGB-luminance importance sampling. | Add explicit InfiniteLight evaluated on misses; keep sphere only for legacy. |
| `src/materials/texture.*` | One RGB-like texture interface plus ad hoc scalar texture classes. | Add `FloatTexture`, `SpectrumTexture`, `TextureEvalContext`, and texture evaluators. |
| `src/materials/constant.*` | A medium is modeled as a stochastic hitable with isotropic material. | Replace in spectral mode with pbrt-style Medium/PhaseFunction and VolPathIntegrator. |
| `src/core/camera.*` | Camera APIs are not wavelength-aware; realistic lenses use scalar IOR. | Add camera methods accepting mutable `SampledWavelengths&` and spectral camera weights. |
| `src/core/adaptivesampler.*` | Accumulates RGB, normal, albedo, and alpha directly. | Move accumulation behind Film; adaptive decisions use sensor/output channels, never packet components. |
| `src/core/buildscene.*` | Parses positional numeric material payloads and hard-coded enums. | Add schema-v2 named descriptors and a SceneCompiler. Keep a v1 adapter. |
| `src/render_scene_rcpp.cpp`, `src/render_animation_rcpp.cpp` | Duplicate scene, environment, camera, and output setup. | Extract one RenderSession so spectral initialization is implemented once. |
| `tools/config/configure.R` | Knows a fixed set of source subdirectories. | Make source discovery recursive or explicitly register every new directory. |
| `R/materials.R` | `dielectric()` documents lower numeric priority as dominant and uses its stack for overlap/void behavior. | Preserve that semantic through `DielectricRegion` and `DielectricPathState`. |
| tests | Primarily seeded legacy render sums. | Preserve them and add unit, statistical, reference-render, and pbrt-comparison suites. |

Do not globally replace `point3f`: most occurrences are geometry. Introduce explicit color/radiometry types and migrate only semantically appropriate values.

---

# 6. Target architecture

```text
R constructors and scene objects
        |
        v
Schema-v2 named descriptors
  SpectrumDesc / TextureDesc / MaterialDesc / LightDesc / MediumDesc
        |
        v
SceneCompiler
  Shapes -> Primitives -> Aggregate
  Materials -> immutable parameter objects
  Lights -> explicit Light handles
  Dielectric solids -> DielectricRegion handles
  Media -> Medium handles
  Camera -> wavelength-aware Camera
  Film -> PixelSensor + filter + output color space
        |
        v
RenderSession
  shared by still and animation entry points
        |
        v
RayIntegrator::EvaluatePixelSample
  SampledWavelengths lambda
  CameraRay weight
  ScratchBuffer
  DielectricPathState
        |
        v
PathIntegrator / VolPathIntegrator
  SampledSpectrum L, beta
  BSDF closures
  explicit LightSampler and MIS
  current Medium on Ray
        |
        v
Film::AddSample(L, lambda, ...)
        |
        v
PixelSensor integration -> XYZ/sensor RGB -> output linear RGB
        |
        v
adaptive preview / denoising / tone mapping / file output
```

## 6.1 Ownership model

- Scene-level spectra, textures, materials, lights, media, shapes, and primitives are immutable after compilation.
- Handles are non-owning references into scene-owned storage.
- Per-hit BxDFs are allocated in a worker-local `ScratchBuffer` and are invalid after the buffer reset.
- `BSDF` is a lightweight wrapper around the local shading frame and one concrete BxDF handle.
- `DielectricPathState` is copied or moved with a path and never aliases another path's mutable state.
- Film tiles own temporary accumulation; Film owns final accumulation.
- Legacy `shared_ptr` ownership remains isolated to legacy adapters until retirement.

## 6.2 Suggested source layout

The exact layout may be adapted, but the contracts must remain separated.

```text
src/spectral/
  sampled_spectrum.h
  sampled_wavelengths.h
  spectrum.h
  spectra.h
  spectra.cpp
  spectral_data.h
  spectral_data.cpp
  rgb_to_spectrum.h
  rgb_to_spectrum.cpp

src/color/
  rgb.h
  xyz.h
  rgb_color_space.h
  rgb_color_space.cpp
  pixel_sensor.h
  pixel_sensor.cpp

src/textures/
  texture_eval_context.h
  texture.h
  textures.h
  textures.cpp
  texture_evaluator.h

src/reflection/
  bxdf.h
  bxdfs.h
  bxdfs.cpp
  bsdf.h
  bsdf.cpp
  fresnel.h
  microfacet.h

src/render/
  interaction.h
  surface_interaction.h
  primitive.h
  primitive.cpp
  scene.h
  scene.cpp
  film.h
  film.cpp
  camera.h
  camera.cpp
  integrator.h
  path_integrator.cpp
  render_session.h
  render_session.cpp

src/lights/
  light.h
  lights.h
  lights.cpp
  light_sampler.h
  light_sampler.cpp

src/media/
  medium.h
  media.h
  media.cpp
  phase_function.h
  phase_functions.cpp

src/materials/
  spectral_material.h
  spectral_materials.h
  spectral_materials.cpp
  dielectric_region.h
  dielectric_region.cpp

src/scene/
  scene_descriptor.h
  scene_compiler.h
  scene_compiler.cpp
  legacy_scene_adapter.cpp
```

All new code MUST be inside `namespace rayrender`.

---

# 7. Core invariants

Codex MUST preserve these invariants in every PR:

1. Legacy render results do not change unless the PR explicitly declares a legacy bug fix.
2. New geometric quantities are not represented by color types, and new color/radiometric quantities are not represented by `point3f`.
3. `SampledSpectrum` and `SampledWavelengths` have fixed-size inline storage and no dynamic allocation.
4. A sampled spectrum is only combined with spectra evaluated at the same wavelength packet.
5. Spectral-to-RGB conversion never occurs inside light transport.
6. Every RGB-to-spectrum conversion has an explicit `SpectrumType`: Albedo, Illuminant, or Unbounded.
7. Encoded color images are decoded once. Scalar/data maps are never gamma decoded.
8. BxDF `f()` excludes the cosine factor.
9. Directional PDFs are scalar and measured with respect to solid angle unless explicitly documented otherwise.
10. `BSDFSample::f`, `pdf`, flags, `eta`, and `pdfIsProportional` follow pbrt semantics.
11. Secondary wavelengths are terminated before a wavelength-dependent direction is selected.
12. A skipped nested-dielectric boundary does not terminate wavelengths, consume a scattering depth, alter prior MIS state, or change `etaScale`.
13. Dielectric-region membership changes only after transmission or a null/skipped crossing, never after reflection.
14. Region entry/exit uses the geometric normal. Bump and shading normals cannot change membership.
15. Segment medium attenuation/scattering is evaluated using the active region before the ending boundary.
16. Alpha rejection occurs before any dielectric-region transition is analyzed or committed.
17. Shadow/visibility traversal uses a local copy of dielectric state and may pass through null boundaries without mutating the camera path.
18. Light-selection probabilities do not depend on the random wavelength packet.
19. Wavelength PDF division occurs exactly once at the sensor.
20. New hot-path interfaces do not own raw pointers.
21. Random sample dimensions have stable documented meanings.
22. Every copied/generated spectral asset has provenance, checksum, and version metadata.
23. R code uses `=` assignment.
24. Source code uses ASCII characters only.


# 8. Normative core interfaces

The declarations in this section are interface sketches, not line-for-line requirements. Codex MAY adjust naming to fit the repository, but MUST preserve semantics, ownership, and evaluation order.

## 8.1 Numeric and color types

Add explicit types before implementing transport:

```cpp
// Reuse rayrender's existing Float alias and configured precision.
struct RGB {
    Float r = 0;
    Float g = 0;
    Float b = 0;
};

struct XYZ {
    Float x = 0;
    Float y = 0;
    Float z = 0;
};

struct RGBColorSpace;
struct RGBColorEncoding;
```

Do not overload geometry types as colors. Existing legacy code MAY continue to use `point3f` internally. New spectral code MUST use geometry-specific vector and point types for geometry and the explicit types above for tristimulus data.

The following conversions MUST be explicit:

- encoded RGB to linear RGB;
- linear RGB to XYZ;
- XYZ to output linear RGB;
- output linear RGB to encoded output RGB;
- RGB to spectral reconstruction, with an explicit `SpectrumType`.

No implicit constructor may silently gamma-decode or convert a color space.

## 8.2 SampledSpectrum

```cpp
class SampledSpectrum {
public:
    SampledSpectrum() = default;
    explicit SampledSpectrum(Float c);

    Float& operator[](int i);
    Float operator[](int i) const;

    bool IsPositive() const;
    bool HasNaNs() const;
    bool IsInf() const;
    Float MaxComponentValue() const;
    Float Average() const;

private:
    std::array<Float, NSpectrumSamples> values_{};
};
```

Provide componentwise arithmetic, `SafeDiv`, `Sqrt`, `Exp`, `ClampZero` (or check if they already exist and are in suitable form), and finiteness helpers. Avoid a broad implicit conversion surface. In particular:

- scalar multiplication and division are allowed;
- spectrum-by-spectrum multiplication and division are componentwise;
- comparisons MUST be named methods rather than ambiguous relational operators;
- negative values are permitted in intermediate sensor/color-transform calculations but radiometric source spectra SHOULD validate nonnegativity where physically required.

## 8.3 SampledWavelengths

```cpp
class SampledWavelengths {
public:
    static SampledWavelengths SampleVisible(Float u);
    static SampledWavelengths SampleUniform(Float u,
                                            Float lambdaMin,
                                            Float lambdaMax);

    Float operator[](int i) const;
    Float PDF(int i) const;

    void TerminateSecondary();
    bool SecondaryTerminated() const;

private:
    std::array<Float, NSpectrumSamples> lambda_{};
    std::array<Float, NSpectrumSamples> pdf_{};
};
```

`SampleUniform()` is for tests and diagnostic render modes only. Production Film sampling defaults to `SampleVisible()`.

Add invariant checks in debug builds:

- wavelengths are finite and within the declared interval;
- PDFs are nonnegative and finite;
- before termination, all PDFs are positive;
- after termination, only the first PDF is positive;
- repeated termination does not alter the packet.

## 8.4 Spectrum hierarchy

Use immutable scene-level spectrum objects. The exact dispatch mechanism MAY be a tagged handle or `std::variant`.

```cpp
enum class SpectrumType {
    Albedo,
    Illuminant,
    Unbounded
};

class Spectrum {
public:
    SampledSpectrum Sample(const SampledWavelengths& lambda) const;
    Float operator()(Float lambdaNm) const;
    Float MaxValue() const;
};
```

Required concrete spectrum forms:

- `ConstantSpectrum`;
- `PiecewiseLinearSpectrum`;
- `BlackbodySpectrum`;
- `DenselySampledSpectrum` for internal tables and named data;
- `RGBAlbedoSpectrum`;
- `RGBUnboundedSpectrum`;
- `RGBIlluminantSpectrum`;

The sampled-data constructor MUST:

- sort wavelengths or reject unsorted input according to an explicit policy;
- reject duplicate wavelengths unless duplicate handling is specified;
- reject nonfinite values;
- define interpolation explicitly;
- define zero, constant, or error extrapolation explicitly;
- record source units and normalization metadata when they matter.

For physical SPDs and optical constants, default extrapolation SHOULD be zero only when the declared dataset intentionally ends outside the renderer's wavelength interval. Optical-constant tables SHOULD use endpoint clamping only after an explicit warning or named policy.

## 8.5 Named spectral data

Add a registry with stable names and provenance for at least:

- CIE 1931 x, y, and z color matching functions;
- the standard illuminants required by supported RGB color spaces;
- a minimal set of measured conductor `eta` and `k` data matching pbrt's named spectra;
- camera sensor response data used by built-in sensors;
- standard glass data added in the camera stage.

Each packaged asset MUST include:

```text
name
semantic type
wavelength unit
value unit or normalization
source citation
source license
asset-generation script version
SHA-256 checksum
```

Do not hide unit conversion in the parser. Convert wavelengths to nanometers in one documented ingestion step.

## 8.6 RGB color spaces and RGB-to-spectrum tables

```cpp
struct RGBColorSpace {
    Point2f r;
    Point2f g;
    Point2f b;
    Point2f white;
    Spectrum illuminant;
    SquareMatrix<3> rgbToXYZ;
    SquareMatrix<3> xyzToRGB;
    RGBToSpectrumTable table;
};
```

The initial production color space is sRGB. Add DCI-P3, Rec.2020, and ACES2065-1 only after sRGB is complete and package-size behavior is measured.

The table representation MUST retain pbrt's resolution and coefficient layout unless the pinned commit differs:

```cpp
inline constexpr int RGBToSpectrumResolution = 64;
// Three max-component cases, a 64 x 64 x 64 grid, and three polynomial coefficients.
using CoefficientArray = /* [3][64][64][64][3] */;
```

Do not reduce table resolution, quantize coefficients, or change interpolation before full conformance and error measurements exist.

The table loader MUST:

- use a versioned binary format with a fixed magic number;
- declare scalar precision, endianness, dimensions, color-space identifier, and checksum;
- reject partial or malformed files;
- load once per process;
- expose no mutable table state after loading;
- work from an installed R package without network access;
- have a deterministic asset-generation script checked into the repository.

The lookup and sigmoid polynomial evaluation MUST numerically match the pinned pbrt implementation within the tolerance recorded in the conformance suite.

The role wrappers MUST match pbrt behavior:

- `RGBAlbedoSpectrum` directly reconstructs the bounded RGB value with the table's sigmoid polynomial.
- `RGBUnboundedSpectrum` uses `scale = 2 * max(rgb)`, looks up the normalized value `rgb / scale`, and multiplies the reconstructed shape by `scale`; black is handled without division.
- `RGBIlluminantSpectrum` uses the same scale/normalized lookup and multiplies the result by the RGB color space's standard illuminant spectrum.

Do not use albedo reconstruction for emitted radiance merely because it stays in `[0, 1]`.

## 8.7 Texture interfaces

Introduce distinct scalar and spectral texture handles:

```cpp
struct TextureEvalContext {
    Point3f p;
    Vector3f dpdx;
    Vector3f dpdy;
    Point2f uv;
    Float dudx = 0;
    Float dudy = 0;
    Float dvdx = 0;
    Float dvdy = 0;
    int faceIndex = -1;
};

class FloatTexture {
public:
    Float Evaluate(const TextureEvalContext& ctx) const;
};

class SpectrumTexture {
public:
    SampledSpectrum Evaluate(const TextureEvalContext& ctx,
                             const SampledWavelengths& lambda) const;
};
```

Implement pbrt-style texture evaluators:

```cpp
class UniversalTextureEvaluator {
public:
    Float operator()(const FloatTexture& texture,
                     const TextureEvalContext& ctx) const;

    SampledSpectrum operator()(const SpectrumTexture& texture,
                               const TextureEvalContext& ctx,
                               const SampledWavelengths& lambda) const;
};
```

A restricted/basic evaluator MAY be added later for constant and image textures, but correctness MUST not depend on it.

Image texture descriptors MUST distinguish:

- color image and its input encoding/color space;
- scalar/data image;
- normal map;
- bump/height map;
- alpha map;
- spectral image or spectral basis asset, when later supported.

Color-image decoding happens before filtering. Spectrum reconstruction happens after filtering the linear RGB texel value. Do not reconstruct every source texel into a dense spectrum at scene-load time unless profiling demonstrates a net benefit.

## 8.8 ScratchBuffer and closure lifetime

Add one `ScratchBuffer` per worker thread. It MUST:

- provide aligned allocation;
- construct non-owning per-hit BxDF objects;
- reset at the same lifetime boundary as pbrt's per-pixel-sample scratch use;
- never be shared concurrently;
- have deterministic high-water-mark instrumentation;
- avoid destructor-dependent resource ownership inside allocated BxDFs.

Per-hit BxDF objects may contain sampled spectral values, scalar roughness, local distributions, and non-owning references to immutable scene objects. They MUST NOT own heap allocations or survive the sample evaluation that created them.

## 8.9 Interactions and primitive separation

Add pbrt-style interaction data without initially rewriting every shape:

```cpp
struct Interaction {
    Point3fi pi;
    Normal3f n;
    Point2f uv;
    Vector3f wo;
    Float time = 0;
    MediumInterface mediumInterface;
};

struct SurfaceInteraction : public Interaction {
    Vector3f dpdu;
    Vector3f dpdv;
    Normal3f dndu;
    Normal3f dndv;

    struct {
        Normal3f n;
        Vector3f dpdu;
        Vector3f dpdv;
        Normal3f dndu;
        Normal3f dndv;
    } shading;

    Material material;
    AreaLight areaLight;
    Primitive primitive;
    RegionHandle dielectricRegion;
};
```

The exact pbrt interval-point type `Point3fi` MAY be deferred if it would force an unrelated robust-intersection rewrite. In that case, document the deviation and use existing rayrender hit-point error offsets until a later robustness PR.

Create adapters around existing `hitable` implementations first. The adapter MUST expose:

- intersection;
- intersection predicate;
- bounds;
- geometric normal;
- UV and differential geometry;
- shape area and direction sampling where supported;
- material, area-light, medium-interface, and dielectric-region bindings.

Do not let a shape call a material or emit light directly in spectral mode.

## 8.10 Ray and medium interface

The spectral ray SHOULD follow pbrt's conceptual state:

```cpp
struct Ray {
    Point3f o;
    Vector3f d;
    Float time = 0;
    Medium medium;
};
```

`RayDifferential` MAY derive from or contain `Ray` as appropriate. The ray does not own wavelengths, throughput, a dielectric-region set, or a sampler.

```cpp
struct MediumInterface {
    Medium inside;
    Medium outside;

    bool IsMediumTransition() const;
};
```

For ordinary pbrt primitives, spawning a ray across a surface updates the medium through `MediumInterface`. For priority dielectric regions, Section 12 determines the effective medium after a committed crossing.

## 8.11 BxDF and BSDF contracts

Use pbrt's flag structure:

```cpp
enum class BxDFFlags : uint32_t {
    Unset = 0,
    Reflection = 1 << 0,
    Transmission = 1 << 1,
    Diffuse = 1 << 2,
    Glossy = 1 << 3,
    Specular = 1 << 4
};

enum class TransportMode {
    Radiance,
    Importance
};

enum class BxDFReflTransFlags : uint32_t {
    Reflection = 1 << 0,
    Transmission = 1 << 1,
    All = (1 << 0) | (1 << 1)
};
```

`BSDFSample` MUST contain at least:

```cpp
struct BSDFSample {
    SampledSpectrum f;
    Vector3f wi;
    Float pdf = 0;
    BxDFFlags flags = BxDFFlags::Unset;
    Float eta = 1;
    bool pdfIsProportional = false;

    bool IsReflection() const;
    bool IsTransmission() const;
    bool IsSpecular() const;
};
```

Each concrete BxDF closure MUST implement the pbrt contract:

```cpp
BxDFFlags Flags() const;
SampledSpectrum f(Vector3f wo,
                  Vector3f wi,
                  TransportMode mode) const;
std::optional<BSDFSample> Sample_f(
    Vector3f wo,
    Float uc,
    Point2f u,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags) const;
Float PDF(Vector3f wo,
          Vector3f wi,
          TransportMode mode,
          BxDFReflTransFlags sampleFlags) const;
void Regularize();
```

The BSDF wrapper owns the shading frame and a non-owning tagged BxDF handle. It transforms between world and local coordinates, handles the geometric/shading normal relationship, and delegates `f`, `Sample_f`, and `PDF`.

Port pbrt's sample-remapping, flags, eta conventions, and roughness-to-alpha mapping exactly. Do not preserve the legacy scatter/PDF split in the new renderer.

## 8.12 Material closure interface

A scene-level Material is immutable and stores texture handles plus options. It does not store per-hit sampled values.

```cpp
struct MaterialEvalContext : public TextureEvalContext {
    Vector3f wo;
    Normal3f n;
    Normal3f ns;
    Vector3f dpdus;
    const ResolvedDielectricInterface* dielectric = nullptr;
};

class Material {
public:
    template <typename TextureEvaluator>
    BSDF GetBSDF(TextureEvaluator texEval,
                 const MaterialEvalContext& ctx,
                 SampledWavelengths& lambda,
                 ScratchBuffer& scratch) const;

    bool CanEvaluateTextures(const UniversalTextureEvaluator&) const;
};
```

Concrete material evaluation SHOULD mirror pbrt's pattern:

```cpp
class DiffuseMaterial {
public:
    using BxDF = DiffuseBxDF;

    template <typename TextureEvaluator>
    BxDF GetBxDF(TextureEvaluator texEval,
                 const MaterialEvalContext& ctx,
                 SampledWavelengths& lambda) const;
};
```

`Material::GetBSDF()` MUST perform the pbrt allocation pattern: dispatch to the concrete Material, obtain a concrete BxDF value from `GetBxDF()`, allocate storage for that concrete type in `ScratchBuffer`, move/copy the value into that storage, and construct a BSDF wrapper over the resulting non-owning BxDF handle. Concrete Material classes do not allocate their own closure storage.

The BxDF receives already evaluated values such as a `SampledSpectrum` reflectance or sampled scalar eta. Consequently, BxDF `f()` does not need a wavelength parameter.

For a nested dielectric material, `MaterialEvalContext::dielectric` provides the effective spectra on both sides after priority resolution. The material MUST NOT inspect or mutate the path's region state.

Material evaluation ordering MUST be:

1. resolve alpha and reject transparent hits;
2. compute geometric and shading interaction data;
3. analyze the dielectric transition if the primitive has a region;
4. construct `MaterialEvalContext` with resolved interface data;
5. evaluate textures and create a BxDF closure;
6. sample/evaluate the closure;
7. commit the transition only if the selected event crosses the boundary.

## 8.13 Initial BxDF set

Implement in this order, using the pinned pbrt formulas:

1. `DiffuseBxDF`;
2. `DiffuseTransmissionBxDF`, if exposed by rayrender's material API;
3. `ConductorBxDF`;
4. `DielectricBxDF`;
5. `ThinDielectricBxDF`;
6. `CoatedDiffuseBxDF`;
7. `CoatedConductorBxDF`;
8. `HairBxDF`;
9. other pbrt BxDFs needed by imported or legacy rayrender materials.

Do not port a BxDF until its standalone numerical and sampling tests exist.

## 8.14 Light interface

Use explicit pbrt-style lights:

```cpp
struct LightSampleContext {
    Point3f p;
    Normal3f n;
    Normal3f ns;
};

struct LightLiSample {
    SampledSpectrum L;
    Vector3f wi;
    Float pdf = 0;
    Interaction pLight;
};

class Light {
public:
    LightType Type() const;

    std::optional<LightLiSample> SampleLi(
        const LightSampleContext& ctx,
        Point2f u,
        const SampledWavelengths& lambda,
        LightSampleMode mode) const;

    Float PDF_Li(const LightSampleContext& ctx,
                 Vector3f wi,
                 LightSampleMode mode) const;

    SampledSpectrum Le(const Ray& ray,
                       const SampledWavelengths& lambda) const;

    SampledSpectrum Phi(const SampledWavelengths& lambda) const;
    void Preprocess(const Bounds3f& sceneBounds);
};
```

Area emission is represented by `AreaLight::L(surfaceInteraction, w, lambda)`. A Primitive binds an area light to geometry. The Material closure has no emission method.

Implement a `LightSampler` handle with at least:

- `UniformLightSampler`;
- `PowerLightSampler` after all light power estimates are correct;
- a spatial sampler only if needed later.

Selection probabilities MUST be scalar and wavelength independent. Power estimates SHOULD use a fixed photometric or standard spectral integral, never the current random packet.

## 8.15 Infinite lights

Replace the spectral environment sphere with an explicit `ImageInfiniteLight` or equivalent. It MUST:

- evaluate `Le(ray, lambda)` on a scene miss;
- sample directions using a wavelength-independent importance distribution;
- include the spherical Jacobian correctly;
- provide a PDF consistent with `SampleLi`;
- support transforms between light and render space;
- reconstruct RGB environment texels as illuminants;
- support directly sampled SPDs for procedural and constant environments;
- participate in both direct-light and emitter-hit MIS.

The internal environment representation SHOULD use the pinned pbrt equal-area square-to-sphere mapping. Existing equirectangular rayrender maps may be converted deterministically at scene load or handled by an explicit input-mapping adapter, but `SampleLi`, `PDF_Li`, and `Le` must all operate through one consistent internal mapping. Do not reuse an equirectangular luminance table with pbrt's equal-area PDF formulas.

The existing sphere implementation remains for legacy only. Do not include it in the spectral aggregate.

## 8.16 Film and PixelSensor

Film owns the wavelength-sampling distribution used by the camera sample:

```cpp
class Film {
public:
    SampledWavelengths SampleWavelengths(Float u) const;

    void AddSample(Point2i pFilm,
                   const SampledSpectrum& L,
                   const SampledWavelengths& lambda,
                   const VisibleSurface* visibleSurface,
                   Float sampleWeight);
};
```

As in pbrt, the caller multiplies the integrator result by `CameraRay::weight` before `AddSample()`. Wavelength-PDF division and sensor-response evaluation remain centralized in PixelSensor/Film.

```cpp
class PixelSensor {
public:
    RGB ToSensorRGB(const SampledSpectrum& L,
                    const SampledWavelengths& lambda) const;
};
```

`ToSensorRGB()` computes an unbiased Monte Carlo estimate of each sensor-channel integral by multiplying by sensor response and dividing each active component by its wavelength PDF. It then applies the sensor normalization used by the pinned pbrt implementation.

Required sensor modes:

- CIE 1931 XYZ-like sensor for predictable colorimetric output;
- pbrt-compatible RGB sensor constructed from a selected RGB color space;
- measured sensor curves in the camera stage.

Film accumulation MUST use enough precision and a concurrency-safe tile or atomic strategy. Packet components MUST never be treated as color channels.

## 8.17 Camera interface

```cpp
struct CameraRay {
    Ray ray;
    SampledSpectrum weight;
};

class Camera {
public:
    std::optional<CameraRay> GenerateRay(
        const CameraSample& sample,
        SampledWavelengths& lambda) const;
};
```

The mutable wavelength packet is required because a realistic camera may choose a wavelength-dependent lens trajectory and terminate secondary wavelengths. Pinhole and thin-lens cameras normally return a wavelength-independent weight and leave the packet intact.

The camera MUST set the ray's initial Medium. Dielectric-region initialization is separate and described in Section 12.13.

## 8.18 Media and phase functions

Introduce pbrt-style Medium and PhaseFunction handles:

```cpp
class Medium {
public:
    SampledSpectrum Tr(const Ray& ray,
                       Float tMax,
                       const SampledWavelengths& lambda,
                       Sampler& sampler) const;

    template <typename Callback>
    SampledSpectrum SampleTmaj(const Ray& ray,
                               Float tMax,
                               Float u,
                               const SampledWavelengths& lambda,
                               Callback callback) const;
};

class PhaseFunction {
public:
    Float p(Vector3f wo, Vector3f wi) const;
    std::optional<PhaseFunctionSample> Sample_p(Vector3f wo,
                                                 Point2f u) const;
};
```

The first Medium implementation is `HomogeneousMedium`, followed by the heterogeneous medium required to replace rayrender's current constant-density hitable behavior. The implementation MUST follow the pinned pbrt majorant/null-collision method rather than adapting RGB stochastic-surface logic.

## 8.19 Integrator interfaces

Use a pbrt-like hierarchy:

```cpp
class Integrator {
public:
    virtual void Render() = 0;
};

class RayIntegrator : public Integrator {
public:
    void Render() override;

protected:
    virtual SampledSpectrum Li(
        RayDifferential ray,
        SampledWavelengths& lambda,
        Sampler& sampler,
        ScratchBuffer& scratch,
        DielectricPathState& regions,
        VisibleSurface* visibleSurface) const = 0;
};
```

Required new integrators:

- `RandomWalkIntegrator`, used as the first correctness vertical slice;
- `PathIntegrator`, matching pbrt surface path behavior;
- `VolPathIntegrator`, after Medium support;
- optional `SimplePathIntegrator` only as a debugging aid.

Do not adapt the current `color()` recursion into the final spectral integrator. A fresh implementation is less risky than preserving inconsistent legacy conventions.

## 8.20 Shape capabilities, implicit shapes, and CSG hooks

The pbrt-like renderer needs a shape contract that advertises capabilities explicitly. This is mandatory because rayrender supports CSG and other geometry that can be intersected but may not support pbrt-style surface sampling, UVs, displacement, or watertight containment.

Add capability metadata to every spectral Shape adapter:

```cpp
struct ShapeCapabilities {
    bool canIntersect = true;
    bool canSampleArea = false;
    bool hasExactArea = false;
    bool hasImplicitInterior = false;
    bool hasReliableContainment = false;
    bool hasUV = false;
    bool supportsNormalMap = false;
    bool supportsDisplacement = false;
};
```

The base Shape contract SHOULD include optional hooks:

```cpp
class Shape {
public:
    Bounds3f Bounds() const;
    std::optional<ShapeIntersection> Intersect(const Ray& ray) const;
    bool IntersectP(const Ray& ray) const;

    ShapeCapabilities Capabilities() const;

    Float Area() const;
    std::optional<ShapeSample> Sample(Point2f u) const;
    std::optional<ShapeSample> Sample(
        const ShapeSampleContext& ctx,
        Point2f u) const;
    Float PDF(const Interaction& intr) const;
    Float PDF(const ShapeSampleContext& ctx, Vector3f wi) const;

    bool Contains(Point3f p) const;
    Float SignedDistance(Point3f p) const;
};
```

If a method is unsupported, it MUST fail through validation before rendering or return an invalid optional result in a controlled way. Do not let a CSG object report a fake `PDF()` of one or a fake random direction in spectral mode.

CSG becomes a custom Shape implementation:

```cpp
class CSGShape final : public Shape {
public:
    Bounds3f Bounds() const;
    std::optional<ShapeIntersection> Intersect(const Ray& ray) const;
    bool IntersectP(const Ray& ray) const;

    ShapeCapabilities Capabilities() const;
    bool Contains(Point3f p) const;
    Float SignedDistance(Point3f p) const;
};
```

Initial CSG capability defaults:

```text
canIntersect = true
canSampleArea = false
hasExactArea = false
hasImplicitInterior = expression-dependent
hasReliableContainment = expression-dependent
hasUV = false
supportsNormalMap = false
supportsDisplacement = false
```

If `csg_mesh_sampler(mode = "render_mesh")` is supplied, the SceneCompiler MAY tessellate the CSG to a triangle mesh and compile the mesh as the rendered Shape. That makes area sampling and rendered intersections refer to the same surface, which is required for unbiased area-light MIS. A sampling proxy for an exact SDF-rendered surface is an advanced diagnostic mode only and is outside strict correctness guarantees.


# 9. R API and schema-v2 scene descriptors

The R API must preserve existing scenes while making every spectral interpretation explicit in new scenes. The public interface is descriptor-based and functional. R users create serializable descriptors for spectra, textures, materials, lights, media, regions, cameras, films, samplers, and integrators. The C++ SceneCompiler turns those descriptors into pbrt-style immutable scene objects and per-hit material closures.

Do not expose C++ BxDF closures to R. R material constructors describe parameters. C++ material evaluation creates pbrt-style per-hit BxDFs in scratch storage.

## 9.1 API design rules

1. Existing rayrender constructors and `render_scene()` arguments remain valid.
2. New spectral controls are provided through typed option descriptors rather than a much larger flat `render_scene()` signature.
3. Every descriptor has a `schema_version`, a `type`, named fields, validation, and a stable class tag.
4. Functional decorators such as `with_light()` and `with_region_boundary()` return modified object rows. They do not mutate global state.
5. Area lights are attached to geometry. Free lights are stored in the scene light registry.
6. Optical regions are independent of materials and shapes. Multiple boundary primitives may reference the same region.
7. The R API uses clear physical terms: radiance, intensity, power, reflectance, eta, sigma_a, sigma_s, sensor, and medium.
8. Legacy convenience remains available through centralized adapters, not through branches in BSDF, Light, Film, or integrator code.

## 9.2 Top-level render API

The preferred spectral call form is:

```r
render_scene(
  scene,
  render_mode = "spectral",
  integrator = path_integrator(max_depth = 50),
  sampler = sobol_sampler(pixel_samples = 512),
  camera = perspective_camera(
    lookfrom = c(278, 278, -800),
    lookat = c(278, 278, 0),
    fov = 40,
    aperture = 0
  ),
  film = rgb_film(
    width = 800,
    height = 800,
    sensor = cie1931_sensor(),
    output_color_space = "sRGB",
    white_balance = "D65",
    save_linear = TRUE
  ),
  spectral = spectral_options(
    sampling = "pbrt_visible_4",
    input_color_space = "sRGB"
  ),
  validation = scene_validation("warn")
)
```

Keep the legacy flat arguments and `integrator_type` for compatibility, but implement a precedence rule:

1. A structured descriptor wins over legacy scalar arguments.
2. Legacy scalar arguments fill missing descriptor fields.
3. Conflicts warn in `scene_validation("warn")` and error in `scene_validation("strict")`.

Recommended additions:

```r
render_scene(
  scene,
  ...,
  render_mode = c("rgb_legacy", "spectral"),
  integrator = NULL,
  sampler = NULL,
  camera = NULL,
  film = NULL,
  spectral = spectral_options(),
  environment = NULL,
  validation = scene_validation(),
  return_result = FALSE
)
```

Legacy arguments such as `width`, `height`, `samples`, `sample_method`, `lookfrom`, `lookat`, `fov`, `aperture`, `environment_light`, `rotate_env`, and `intensity_env` compile into the corresponding descriptors when the structured descriptor is absent.

## 9.3 Render-option descriptors

Implement these S3 descriptor families:

```r
path_integrator = function(
  max_depth = 50,
  regularize = FALSE,
  light_sampler = uniform_light_sampler(),
  russian_roulette = TRUE
)

volpath_integrator = function(
  max_depth = 50,
  regularize = FALSE,
  light_sampler = uniform_light_sampler()
)

random_walk_integrator = function(max_depth = 50)

sobol_sampler = function(pixel_samples = 128, randomize = TRUE)
stratified_sampler = function(x_samples = 4, y_samples = 4, jitter = TRUE)
independent_sampler = function(pixel_samples = 128)

perspective_camera = function(
  lookfrom = c(0, 1, -10),
  lookat = c(0, 0, 0),
  up = c(0, 1, 0),
  fov = 20,
  aperture = 0,
  focal_distance = NULL,
  shutteropen = 0,
  shutterclose = 1,
  initial_regions = character()
)

orthographic_camera = function(...)
realistic_camera = function(...)

rgb_film = function(
  width = NULL,
  height = NULL,
  sensor = cie1931_sensor(),
  output_color_space = "sRGB",
  white_balance = NULL,
  filename = NULL,
  save_linear = TRUE,
  alpha = TRUE
)

cie1931_sensor = function(output_color_space = "sRGB", white_balance = NULL)
rgb_sensor = function(color_space = "sRGB", white_balance = NULL)
measured_sensor = function(red, green, blue, imaging_ratio = 1, white_balance = NULL)

spectral_options = function(
  sampling = "pbrt_visible_4",
  wavelength_range = c(360, 830),
  input_color_space = "sRGB",
  rgb_numeric_encoding = "legacy",
  rgb_table_resolution = 64,
  table_cache = TRUE
)

scene_validation = function(
  mode = c("warn", "strict", "advanced", "none"),
  check_watertight = TRUE,
  check_region_consistency = TRUE,
  check_spectrum_bounds = TRUE,
  check_light_geometry = TRUE,
  allow_open_regions = FALSE,
  allow_unsampled_emitters = FALSE
)
```

These descriptors map to pbrt concepts: Integrator, Sampler, Camera, Film, PixelSensor, and wavelength sampling. The Film/PixelSensor descriptor owns wavelength-to-tristimulus conversion. The integrator never converts spectra to RGB.

## 9.4 Spectrum constructors

Add immutable R spectrum descriptors:

```r
spectrum_constant = function(value)

spectrum_rgb = function(
  value,
  role = c("albedo", "illuminant", "unbounded"),
  color_space = "sRGB",
  encoding = NULL,
  scale = 1
)

spectrum_sampled = function(
  wavelength_nm,
  value,
  role = NULL,
  interpolation = "linear",
  extrapolation = "zero",
  normalize = "none",
  units = NULL,
  scale = 1
)

spectrum_blackbody = function(
  temperature,
  scale = 1,
  normalize = TRUE
)

spectrum_named = function(name, scale = 1)
spectrum_cauchy_ior = function(A, B = 0, C = 0)
spectrum_sellmeier_ior = function(B, C)
```

Use `role` exactly as pbrt's RGB spectrum classes require:

```r
spectrum_rgb("#884422", role = "albedo")
spectrum_rgb("#fff0d0", role = "illuminant")
spectrum_rgb(c(0.1, 0.02, 0.005), role = "unbounded")
```

Compatibility mappings:

```r
diffuse(color = "red")
```

compiles as:

```r
diffuse(reflectance = spectrum_rgb("red", role = "albedo"))
```

and legacy emissive material input compiles as an explicit area light with illuminant reconstruction.

## 9.5 Texture constructors

Use explicit scalar and spectral texture descriptors. Do not infer physical semantics from the number of image channels.

```r
texture_constant = function(value, value_type = c("auto", "float", "spectrum"))

texture_image = function(
  filename,
  value_type = c("spectrum", "float"),
  role = NULL,
  color_space = "sRGB",
  encoding = "auto",
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1,
  invert = FALSE
)

texture_image_color = function(
  filename,
  color_space = "sRGB",
  encoding = "srgb",
  role = c("albedo", "illuminant", "unbounded"),
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1
)

texture_image_scalar = function(
  filename,
  channel = "luminance",
  encoding = "linear",
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1
)

texture_checker = function(tex1, tex2, dimension = 2, uscale = 1, vscale = 1)
texture_mix = function(tex1, tex2, amount = 0.5)
texture_scale = function(tex, scale = 1)
texture_float = function(value, range = NULL, encoding = "linear")
```

CSG and other shapes without UVs may accept generated mappings such as:

```r
triplanar_mapping = function(scale = 1, blend = 0.2)
object_mapping = function(transform = identity_transform())
```

Strict validation MUST reject UV-only image textures, normal maps, and displacement maps on CSG unless a supported generated mapping or tessellated replacement is supplied.

## 9.6 Material constructors

Keep existing public names where possible, but make the spectral descriptors pbrt-shaped. Recommended primary constructors:

```r
diffuse = function(
  reflectance = spectrum_rgb(c(0.5, 0.5, 0.5), role = "albedo"),
  color = NULL,
  image_texture = NULL,
  displacement = NULL,
  normal_map = NULL
)

conductor = function(
  eta = spectrum_named("metal-Ag-eta"),
  k = spectrum_named("metal-Ag-k"),
  roughness = 0,
  u_roughness = NULL,
  v_roughness = NULL,
  remap_roughness = TRUE,
  displacement = NULL,
  normal_map = NULL
)

dielectric_interface = function(
  roughness = 0,
  u_roughness = NULL,
  v_roughness = NULL,
  remap_roughness = TRUE,
  displacement = NULL,
  normal_map = NULL
)

thin_dielectric = function(
  eta = spectrum_constant(1.5),
  displacement = NULL,
  normal_map = NULL
)

coated_diffuse = function(...)
coated_conductor = function(...)
diffuse_transmission = function(...)
hair = function(...)
measured_material = function(filename)
mix_material = function(mat1, mat2, amount = 0.5)
interface_material = function()
null_material = interface_material
```

Legacy aliases remain:

```r
metal(...)
dielectric(refraction = 1.5, attenuation = NULL, priority = 0, ...)
light(color = "white", intensity = 1, ...)
```

In `render_mode = "spectral"`, `light()` is a schema-v1 shorthand for an area light and SHOULD warn once per scene:

```text
`light()` as a material is a legacy shorthand. In spectral mode, prefer `with_light(area_light(...))`.
```

Eta ownership is unambiguous:

- `thin_dielectric()` owns its eta on the material and never binds a region.
- `dielectric_interface()` obtains effective eta from the region resolver when attached to a region boundary.
- A convenience `dielectric()` may create both `dielectric_interface()` and an automatic `optical_region()`, but the compiled spectral scene stores bulk eta on the region.
- Bulk absorption is a Medium property, not a surface-material color.

## 9.7 Medium and optical-region constructors

```r
homogeneous_medium = function(
  sigma_a = spectrum_constant(0),
  sigma_s = spectrum_constant(0),
  scale = 1,
  phase = henyey_greenstein_phase(0),
  emission = spectrum_constant(0)
)

henyey_greenstein_phase = function(g = 0)

optical_region = function(
  id,
  eta = spectrum_constant(1),
  priority = 0L,
  medium = NULL
)

region_boundary = function(
  region,
  side = c("negative_normal", "positive_normal", "inside")
)
```

`side = "inside"` is valid only for shapes with a reliable implicit or analytic containment predicate, including supported closed CSG shapes. Mesh boundary patches normally use `positive_normal` or `negative_normal`.

For compatibility, a convenience constructor MAY remain:

```r
dielectric_region = function(
  eta = spectrum_constant(1.5),
  priority = 0L,
  sigma_a = spectrum_constant(0),
  medium = NULL,
  id = NULL
)
```

If both `sigma_a` and `medium` are supplied, error unless one is exactly the default.

## 9.8 Light constructors

Lights are explicit descriptors. New spectral Material objects do not emit radiance.

Scene-level free lights:

```r
point_light = function(
  position,
  intensity = NULL,
  power = NULL,
  emission = spectrum_rgb("white", role = "illuminant"),
  scale = 1
)

spot_light = function(
  from,
  to,
  intensity = NULL,
  emission = spectrum_rgb("white", role = "illuminant"),
  cone_angle = 30,
  cone_delta_angle = 5,
  scale = 1
)

distant_light = function(
  direction,
  radiance = spectrum_rgb("white", role = "illuminant"),
  scale = 1
)

projection_light = function(...)
goniometric_light = function(...)
```

Infinite lights:

```r
uniform_infinite_light = function(
  radiance = spectrum_rgb("white", role = "illuminant"),
  scale = 1
)

image_infinite_light = function(
  filename,
  scale = 1,
  rotation = 0,
  encoding = "auto",
  color_space = "sRGB",
  importance_sample = TRUE
)
```

Area lights:

```r
area_light = function(
  emission = spectrum_rgb("white", role = "illuminant"),
  radiance = NULL,
  scale = 1,
  power = NULL,
  filename = NULL,
  two_sided = FALSE,
  visible = TRUE,
  sampling = c("auto", "sampled", "none")
)
```

The constructor documentation MUST state the supplied physical quantity. For area lights, the primary quantity is emitted radiance. For point lights, the primary quantity is spectral intensity or power. Do not use one generic `intensity` field in the new API unless its physical meaning is light-type specific and explicitly documented.

## 9.9 Object decorators: `with_light()` and `with_region_boundary()`

`with_light()` is an object decorator. It attaches an area-light descriptor to existing geometry rows and returns the modified object rows. It does not create a separate geometric object and it does not turn the material into an emitter.

Example:

```r
ceiling_panel = xz_rect(
  x = 278,
  y = 554,
  z = 278,
  xwidth = 130,
  zwidth = 105,
  material = interface_material()
) |>
  with_light(area_light(
    emission = spectrum_blackbody(temperature = 4000, scale = 20),
    two_sided = FALSE
  ))
```

Conceptually this changes the object row from:

```r
list(
  shape = "xz_rect",
  material = interface_material(),
  light = NULL
)
```

to:

```r
list(
  shape = "xz_rect",
  material = interface_material(),
  light = area_light(
    emission = spectrum_blackbody(temperature = 4000, scale = 20),
    two_sided = FALSE
  )
)
```

The spectral SceneCompiler then emits:

```text
Shape: XZ rectangle
Material: interface material
AreaLight: DiffuseAreaLight bound to the same primitive
Primitive: GeometricPrimitive(shape, material, areaLight, region metadata)
```

This is the R equivalent of pbrt's area light associated with shape geometry. pbrt uses scene-file attribute scoping; rayrender uses a functional decorator because object rows are explicit R values.

Required implementation sketch:

```r
with_light = function(object, light, replace = TRUE) {
  if (!inherits(light, "ray_light")) {
    stop("`light` must be a ray_light object.")
  }
  if (!isTRUE(attr(light, "area"))) {
    stop("`with_light()` only accepts area lights. Use `add_light()` for point, spot, distant, and infinite lights.")
  }

  object = ensure_scene_v2_columns(object)
  has_light = !vapply(object$light, is.null, logical(1))
  if (any(has_light) && !replace) {
    stop("Object already has an area light. Use `replace = TRUE` to replace it.")
  }

  object$light = rep(list(light), nrow(object))
  object
}
```

`with_region_boundary()` is the same pattern for optical-region metadata:

```r
with_region_boundary = function(object, boundary, append = TRUE) {
  object = ensure_scene_v2_columns(object)

  if (inherits(boundary, "ray_region_boundary")) {
    boundary = list(boundary)
  }

  if (!append) {
    object$region_boundaries = rep(list(boundary), nrow(object))
  } else {
    object$region_boundaries = Map(
      function(old) c(old, boundary),
      object$region_boundaries
    )
  }

  object
}
```

A single visible primitive may simultaneously carry a Material, an AreaLight, and one or more region-boundary attachments.

## 9.10 Scene-level registries and functional scene operations

Area lights use `with_light()`. Free lights and environments are scene-level registries:

```r
add_light = function(scene, light)
set_environment = function(scene, environment)
add_region = function(scene, region)
add_named_texture = function(scene, name, texture)
add_named_material = function(scene, name, material)
named_texture = function(name)
named_material = function(name)
```

Example:

```r
scene = ray_scene() |>
  add_region(optical_region(
    id = "glass_ball",
    eta = spectrum_named("glass-BK7"),
    priority = 1
  )) |>
  add_object(sphere(
    x = 0,
    y = 1,
    z = 0,
    radius = 1,
    material = dielectric_interface()
  ) |>
    with_region_boundary(region_boundary(
      region = "glass_ball",
      side = "negative_normal"
    ))) |>
  add_light(point_light(
    position = c(0, 5, -3),
    emission = spectrum_blackbody(temperature = 5500, scale = 1),
    power = 100
  )) |>
  set_environment(image_infinite_light(
    filename = "studio.exr",
    rotation = 30,
    scale = 0.5
  ))
```

Named registries are optional for small scenes but required for importers, large scenes, repeated resources, pbrt export, and stable diagnostics.

## 9.11 Multi-primitive optical-region example

A water volume may be bounded by separate meshes. Do not require a duplicate coincident seafloor cap.

```r
water_region = optical_region(
  id = "lake_water",
  priority = 10,
  eta = spectrum_named("water"),
  medium = homogeneous_medium(
    sigma_a = spectrum_sampled(
      wavelength_nm = c(400, 450, 500, 550, 600, 650, 700),
      value = c(0.01, 0.006, 0.004, 0.01, 0.04, 0.12, 0.30),
      role = "unbounded"
    )
  )
)

scene = ray_scene() |>
  add_region(water_region) |>
  add_object(
    mesh3d_model(
      terrain_mesh,
      material = diffuse(
        reflectance = texture_image(
          filename = "terrain_albedo.png",
          value_type = "spectrum",
          role = "albedo"
        )
      )
    ) |>
      with_region_boundary(region_boundary(
        region = "lake_water",
        side = "positive_normal"
      ))
  ) |>
  add_object(
    mesh3d_model(
      water_surface_mesh,
      material = dielectric_interface(roughness = 0)
    ) |>
      with_region_boundary(region_boundary(
        region = "lake_water",
        side = "negative_normal"
      ))
  )
```

The seafloor material remains the terrain material. Reflection from the seafloor leaves the ray in water. Transmission through a transmissive seafloor would atomically remove water and add whatever region is declared on the other side.

## 9.12 CSG API additions

CSG remains a rayrender-specific extension. Ordinary visible CSG works through the Shape/Primitive/Material contracts. Region-bearing CSG requires a coherent implicit inside. Area-light CSG requires a sampleable surface.

Extend `csg_object()` with metadata rather than changing CSG construction syntax wholesale:

```r
csg_object = function(
  object,
  x = 0,
  y = 0,
  z = 0,
  material = diffuse(),
  angle = c(0, 0, 0),
  order_rotation = c(1, 2, 3),
  flipped = FALSE,
  scale = c(1, 1, 1),
  topology = c("auto", "closed", "open", "unknown"),
  surface_sampler = NULL,
  texture_mapping = NULL,
  region_inside = TRUE
)

csg_mesh_sampler = function(
  resolution = 128,
  max_error = NULL,
  mode = c("render_mesh", "sampling_proxy"),
  bounds = NULL,
  normal_method = c("mesh", "sdf_gradient"),
  seed = 1
)
```

Default policy:

- CSG diffuse, conductor, dielectric-interface, and null-boundary surfaces are allowed when the shape can intersect robustly.
- CSG optical regions use `region_boundary(side = "inside")` and require a coherent inside.
- CSG as an area light is rejected in strict spectral mode unless `surface_sampler = csg_mesh_sampler(mode = "render_mesh")` or an equivalent exact sampleable surface is provided.
- CSG UV image textures, normal maps, and displacement are rejected unless generated mapping support is provided.
- Per-child CSG materials and per-child regions are out of scope for the first spectral implementation.

Example CSG glass region:

```r
glass_region = optical_region(
  id = "implicit_glass",
  eta = spectrum_named("glass-BK7"),
  priority = 1
)

scene = ray_scene() |>
  add_region(glass_region) |>
  add_object(
    csg_object(
      csg_combine(
        csg_sphere(radius = 1),
        csg_sphere(x = 0.4, radius = 0.5),
        operation = "subtract"
      ),
      material = dielectric_interface()
    ) |>
      with_region_boundary(region_boundary(
        region = "implicit_glass",
        side = "inside"
      ))
  )
```

Example CSG area light through tessellation:

```r
scene = ray_scene() |>
  add_object(
    csg_object(
      csg_round(
        csg_box(width = c(2, 0.25, 2)),
        radius = 0.08
      ),
      material = interface_material(),
      surface_sampler = csg_mesh_sampler(
        resolution = 192,
        mode = "render_mesh",
        max_error = 0.002
      )
    ) |>
      with_light(area_light(
        emission = spectrum_blackbody(temperature = 3200, scale = 15),
        two_sided = FALSE
      ))
  )
```

The `render_mesh` mode is the strict-mode default because the rendered emitting surface and sampled emitting surface are identical. A sampling proxy for a different rendered implicit surface is not pbrt-conformant unless a future proof and tests show unbiased MIS behavior.

## 9.13 `ray_scene` S3 changes

Do not replace `ray_scene` with a non-data-frame container in the first API migration. Existing workflows use object rows, `nrow(scene)`, `rbind()`, and pipeable constructors. Instead, implement a hybrid scene:

```r
class(scene) = c("ray_scene_v2", "ray_scene", "tbl_df", "tbl", "data.frame")
```

Required row-level list columns and metadata:

```r
x
y
z
shape
material
shape_info
transforms
animation_info
light
region_boundaries
medium_interface
shape_capabilities
object_id
object_name
visibility
user_data
```

Scene-level attributes:

```r
attr(scene, "ray_schema_version") = 2L
attr(scene, "regions") = list()
attr(scene, "lights") = list()
attr(scene, "environment") = NULL
attr(scene, "named_textures") = list()
attr(scene, "named_materials") = list()
attr(scene, "render_defaults") = list()
attr(scene, "validation") = scene_validation()
```

`add_object()` MUST preserve attributes and registries:

```r
add_object = function(scene, objects = NULL) {
  if (is.null(objects)) {
    return(scene)
  }

  scene = ensure_ray_scene_v2(scene)
  objects = ensure_ray_scene_v2(objects)

  scene_attrs = ray_scene_attrs(scene)
  object_attrs = ray_scene_attrs(objects)

  newscene = rbind(unclass(scene), unclass(objects))
  class(newscene) = class(scene)

  newscene = restore_ray_scene_attrs(
    newscene,
    merge_ray_scene_attrs(scene_attrs, object_attrs)
  )

  assign_missing_object_ids(newscene)
}
```

`print.ray_scene()` SHOULD report object count, area lights, free lights, infinite lights, regions, media, material summary, and schema version. It must not count lights by checking whether `material$type` is a legacy light enum.

## 9.14 `ray_material` S3 changes

The current flat material payload is too coupled to RGB-era C++ enum handling. Add a schema-v2 descriptor:

```r
new_ray_material = function(
  type,
  params = list(),
  textures = list(),
  normal_map = NULL,
  displacement = NULL,
  flags = list(),
  legacy = list()
) {
  vctrs::new_vctr(
    list(list(
      schema_version = 2L,
      type = type,
      params = params,
      textures = textures,
      normal_map = normal_map,
      displacement = displacement,
      flags = flags,
      legacy = legacy
    )),
    class = c("ray_material_v2", "ray_material")
  )
}
```

Store type strings on the R side, for example `"diffuse"`, `"conductor"`, and `"dielectric_interface"`. The C++ compiler may map them to compact enums later. Keeping strings in R makes validation, printing, serialized scenes, conversion reports, and importer diagnostics much safer.

R material objects remain serializable descriptors. They do not contain runtime BxDF closures or R closures. C++ handles this mapping:

```text
ray_material_v2 descriptor
  -> SpectralMaterial scene object
  -> Material::GetBSDF(textureEvaluator, ctx, lambda, scratch)
  -> concrete BxDF
  -> BSDF
```

Add two central adapters:

```r
material_to_spectral = function(material)
material_to_legacy = function(material)
```

No spectral BxDF, Light, Film, or integrator code may inspect legacy material payloads directly.

## 9.15 Schema-v1 adaptation

Create one adapter that converts legacy scene objects into schema v2 before spectral compilation. It MUST be centralized.

Important mappings:

| Legacy input | Spectral adaptation |
|---|---|
| surface RGB color | `SpectrumType::Albedo` unless the parameter is explicitly emissive |
| light RGB color | `SpectrumType::Illuminant` |
| dielectric scalar `refraction` | constant eta spectrum |
| dielectric RGB `attenuation` | unbounded spectral `sigma_a` compatibility mapping on a Medium |
| dielectric `priority` | `DielectricRegion::priority` |
| metal RGB color with fuzz | compatibility conductor closure with warning |
| environment RGB/image | illuminant reconstruction |
| scalar textures | FloatTexture without color decoding |
| legacy medium color/density | homogeneous sigma values through a documented mapping |
| CSG object | CSGShape with capability flags and optional implicit containment |

Legacy nonwhite dielectric surface color is not a pbrt dielectric parameter. For compatibility, map it to a separately named internal `LegacyTintedDielectricMaterial` and emit one warning per scene in spectral mode. New APIs MUST not expose this parameter.

## 9.16 Descriptor validation

Validate before entering C++ hot paths:

- finite numeric values;
- nonnegative radiometric spectra where required;
- albedo values in the allowed range;
- positive eta;
- valid roughness and anisotropy ranges;
- unique explicit region IDs;
- integer priority range;
- resolved region-boundary references;
- area-light attachment only to sampleable shapes or explicitly allowed advanced unsampled emitters;
- CSG optical-region use only for coherent implicit interiors;
- image file existence and supported encoding;
- consistent sampled-spectrum vector lengths;
- monotonic wavelength samples;
- supported source and output color spaces.

Error messages MUST identify the scene object and parameter path, for example:

```text
object[17].region["lake_water"].eta: all IOR samples must be positive
object[8].light: CSG area lights require csg_mesh_sampler(mode = "render_mesh") in strict spectral mode
```

Do not defer obvious schema errors until a worker thread is rendering.

# 10. Material migration and closure mapping

## 10.1 General migration rule

For each existing rayrender material:

1. identify whether its parameters represent surface scattering, emission, or bulk transport;
2. map it to an existing pbrt-v4 Material/BxDF when possible;
3. evaluate its textures once per surface hit into a concrete BxDF closure;
4. add independent `f`, `Sample_f`, PDF, reciprocity, and energy tests;
5. preserve the old implementation only in `rgb_legacy`;
6. document intentional behavior changes in spectral mode.

Do not first create a spectral version of the current `material::scatter()` interface and then replace it later. That would duplicate the riskiest part of the project.

## 10.2 Mapping table

| Current rayrender material | Spectral target | Notes |
|---|---|---|
| `lambertian` | pbrt `DiffuseMaterial` plus `DiffuseBxDF` | Reflectance is albedo-role spectrum. |
| `orennayar` | pbrt diffuse-rough model or a separately validated Oren-Nayar BxDF | Prefer the exact pinned pbrt model if present; otherwise retain only after convention tests. |
| `metal` | pbrt `ConductorMaterial`/`ConductorBxDF` | New API uses spectral eta/k. Legacy RGB plus fuzz uses a compatibility mapping. |
| `dielectric` | pbrt `DielectricMaterial`/`DielectricBxDF` plus region resolver | Eta comes from effective inside/outside region spectra. Bulk attenuation moves to Medium. |
| thin glass behavior, if any | pbrt `ThinDielectricMaterial`/`ThinDielectricBxDF` | Does not enter region membership. |
| `MicrofacetReflection` | pbrt conductor or coated closure | Remove duplicate standalone legacy model from spectral path. |
| `MicrofacetTransmission` | pbrt rough dielectric closure | Use one combined dielectric BxDF, not separately sampled reflection/transmission classes. |
| `glossy` | pbrt coated diffuse or coated conductor | Choose based on documented legacy semantics; do not invent a lobe mixture without tests. |
| `hair` | pbrt `HairMaterial`/`HairBxDF` | Evaluate spectral absorption and follow pbrt longitudinal/azimuthal sampling. |
| `diffuse_light` | `DiffuseAreaLight` attached to Primitive | Optional separate nonemissive Material for the same surface. |
| `spot_light` material | explicit `SpotLight` | Remove from Material in spectral scenes. |
| `isotropic` | phase function within Medium | Not a surface BxDF in VolPath. |
| null/interface material | pbrt interface/null boundary behavior | Used for medium and region boundaries; consumes no scattering depth. |

## 10.3 Diffuse closure

The diffuse vertical slice MUST implement:

```cpp
class DiffuseBxDF {
public:
    explicit DiffuseBxDF(SampledSpectrum R);

    BxDFFlags Flags() const;
    SampledSpectrum f(Vector3f wo,
                      Vector3f wi,
                      TransportMode mode) const;
    std::optional<BSDFSample> Sample_f(... ) const;
    Float PDF(... ) const;

private:
    SampledSpectrum R_;
};
```

Required behavior:

- `f = R / Pi` for same-hemisphere directions;
- cosine-weighted hemisphere sampling;
- `pdf = AbsCosTheta(wi) / Pi`;
- no cosine embedded in `f`;
- reflectance texture clamped/validated according to pbrt's policy;
- shading frame built consistently from the evaluated surface interaction.

## 10.4 Conductor closure

Use pbrt's conductor Fresnel and Trowbridge-Reitz implementation. The closure captures:

```cpp
SampledSpectrum eta;
SampledSpectrum k;
TrowbridgeReitzDistribution distribution;
```

Rules:

- reflection direction is shared by all packet wavelengths, so conductor dispersion does not terminate secondary wavelengths;
- Fresnel is evaluated componentwise at each sampled wavelength;
- smooth and rough limits match pbrt;
- eta/k datasets are sampled at current wavelengths during material evaluation;
- roughness remapping follows pbrt exactly;
- no arbitrary RGB multiplier is added in the new API.

## 10.5 Dielectric closure

Use pbrt's combined `DielectricBxDF` for reflection and transmission. The closure captures one scalar relative eta evaluated at `lambda[0]` plus a microfacet distribution.

For ordinary non-region primitives:

- evaluate the material's eta spectrum at `lambda[0]`;
- if eta is nonconstant, call `TerminateSecondary()` before constructing/sampling the closure;
- use pbrt's eta convention and Fresnel branch probabilities.

For a region boundary:

- obtain effective outside and inside eta spectra from `ResolvedDielectricInterface`;
- determine whether their relative ratio is wavelength dependent;
- terminate secondary wavelengths if required;
- evaluate `etaRelative = etaInside(lambda[0]) / etaOutside(lambda[0])` relative to the primitive geometric normal;
- construct the same pbrt `DielectricBxDF` with that scalar ratio;
- let the integrator commit or reject the region transition based on the returned event.

Bulk absorption, scattering, or emission are not captured by `DielectricBxDF`.

## 10.6 Thin dielectric closure

`ThinDielectricBxDF` models a thin sheet and MUST NOT:

- add or remove a dielectric region;
- change the ray's current Medium;
- apply finite-thickness absorption automatically;
- update a nested-region stack.

Follow the pinned pbrt-v4 `ThinDielectricMaterial` wavelength rule exactly: evaluate nonconstant eta at `lambda[0]` and call `TerminateSecondary()` before constructing the thin-dielectric closure. Constant eta preserves the wavelength packet. This termination is part of material evaluation only; it does not imply a region-state transition.

## 10.7 Coated closures

Port pbrt's dedicated coated BxDFs rather than mixing diffuse and glossy samples with an ad hoc probability. Preserve:

- layer depth sampling;
- interface Fresnel behavior;
- transport mode;
- PDF bookkeeping;
- regularization behavior;
- sample-count constants and any pbrt-specific approximations.

Expose physical parameter names in R. Keep legacy `glossy()` through an adapter with a documented mapping.

## 10.8 Hair closure

Port pbrt's Hair BxDF after core surface validation. Required spectral parameters include absorption `sigma_a` or one of pbrt's alternate parameterizations. The closure MUST capture the sampled absorption spectrum and scalar geometric parameters.

Add tests for:

- white-furnace energy behavior;
- azimuthal normalization;
- sampling/PDF agreement;
- melanin or color-to-absorption conversion if exposed;
- spectral absorption under narrow-band illumination.

## 10.9 Bump and displacement evaluation

Follow pbrt's material-evaluation order:

1. evaluate displacement/bump using a FloatTexture;
2. modify only shading geometry;
3. retain the original geometric normal for visibility, medium transitions, and dielectric region entry/exit;
4. build the BSDF frame from shading geometry;
5. apply pbrt's shading-normal correction through its standard conventions, not a custom spectral correction.

## 10.10 Emissive surfaces

A spectral emissive surface has:

- a Primitive with shape;
- an optional nonemissive Material for reflected light;
- a `DiffuseAreaLight` or another explicit area-light implementation.

The same sampled emission spectrum MUST be returned consistently by:

- direct-light sampling;
- `AreaLight::L()` when a path hits the emitter;
- power estimation for light selection;
- any AOV that reports emission.

# 11. Canonical pbrt-style path integration

## 11.1 Per-sample setup

The RenderSession and RayIntegrator perform the following for each pixel sample:

1. obtain filter, lens, time, wavelength, and auxiliary random samples from the worker sampler;
2. call `film.SampleWavelengths(uLambda)`;
3. pass the mutable packet to `camera.GenerateRay()`;
4. initialize `DielectricPathState` for the camera origin;
5. reset the worker ScratchBuffer;
6. evaluate the selected integrator;
7. multiply the integrator result by `CameraRay::weight`, then pass the result, wavelengths, filter weight, and visible-surface information to Film;
8. reset scratch storage and continue.

Random dimensions MUST be named in a developer document and stable across platforms for deterministic tests.

## 11.2 Path state

The surface PathIntegrator should mirror pbrt's state as closely as practical:

```cpp
SampledSpectrum L(0);
SampledSpectrum beta(1);
Float etaScale = 1;
bool specularBounce = false;
bool anyNonSpecularBounces = false;

Float pBsdf = 1;
LightSampleContext previousLightContext;

DielectricPathState regions;
```

The exact names may follow the pinned source. Do not collapse previous-light MIS state into a generic recursive return value.

## 11.3 Surface path loop

Implement the following order, with the pinned pbrt source open during coding:

1. Intersect the ray with the aggregate.
2. If the ray has no surface hit, accumulate infinite-light emission. Apply emitter-hit MIS unless this is the camera ray or the previous event was specular.
3. Resolve alpha masking before treating a hit as a real boundary.
4. Add area-light emission at the hit, with emitter-hit MIS under the same pbrt conditions.
5. Apply medium transmittance or medium sampling in VolPath; surface Path uses vacuum or absorption-only handling only as explicitly implemented.
6. Analyze any nested-dielectric region transition, but do not commit it.
7. Compute differentials, bump mapping, and the per-hit BSDF closure in ScratchBuffer.
8. If the boundary is null/skipped, commit its transition, spawn the unchanged-direction ray, and continue without consuming scattering depth or overwriting prior scattering MIS state.
9. If the BSDF has non-specular components, sample direct illumination through the LightSampler and apply MIS.
10. Sample the BSDF closure.
11. Update throughput with `f * AbsDot(wi, ns) / pdf` using the pinned shading-normal convention.
12. Update prior BSDF PDF and previous light context for possible emitter-hit MIS.
13. If the sampled event transmitted through a region boundary, commit its region transition and update the ray Medium.
14. If transmission occurred, update `etaScale` using `bs.eta` exactly as pbrt does.
15. Spawn the new ray using robust origin offsets.
16. Apply regularization if enabled and appropriate.
17. Apply Russian roulette using `beta * etaScale` and pbrt's depth policy.
18. Continue until maximum depth, zero throughput, invalid sample, or roulette termination.

## 11.4 Throughput convention

For a non-delta surface sample:

```cpp
beta *= bs.f * AbsDot(bs.wi, isect.shading.n) / bs.pdf;
```

Use the exact normal choice and correction from the pinned pbrt BSDF implementation. Never multiply by an additional cosine hidden in a legacy material `f()`.

For proportional PDFs, preserve pbrt's `pdfIsProportional` handling before using the sample PDF in MIS.

## 11.5 Direct-light MIS

For a selected non-delta light:

```text
pLight = light-selection PMF * light directional PDF
pBsdf = BSDF directional PDF
weight = PowerHeuristic(1, pLight, 1, pBsdf)
contribution = beta * f * Li * abs(dot(wi, ns)) * weight / pLight
```

For delta lights, do not evaluate a competing BSDF PDF. Visibility/transmittance is evaluated after the sample and uses local path state as described in Section 12.12.

Light-sampling code MUST reject:

- zero or nonfinite light PDF;
- zero or nonfinite selection PMF;
- black/nonfinite sampled radiance;
- a black BSDF value;
- blocked visibility.

Do not clamp valid high contributions in the integrator. Firefly controls belong in explicit optional policies, not correctness code.

## 11.6 Emitter-hit MIS

When a BSDF-sampled ray reaches an area light or infinite light:

- use weight 1 for the camera ray or after a specular event;
- otherwise compute the selected light's PMF and directional PDF from the previous light-sample context;
- combine the stored BSDF PDF and light PDF with pbrt's power heuristic;
- include every relevant infinite light on a miss according to the pinned implementation.

A null/skipped dielectric crossing does not replace the previous scattering event. Its traversal must preserve the prior BSDF PDF, prior light context, and specular flag.

## 11.7 Russian roulette

Match pbrt's roulette policy and `etaScale` use. The continuation probability is based on scaled throughput, not raw packet average alone. The implementation MUST be unbiased and tested statistically.

Do not use wavelength PDF in roulette probability; it belongs to the sensor estimator.

## 11.8 Regularization

Implement pbrt's optional BSDF regularization only after unregularized parity. Record whether a path has already had non-specular scattering and call `Regularize()` under the same conditions as the pinned integrator.

`regularize = FALSE` remains the correctness-reference default.

## 11.9 RandomWalkIntegrator

Before MIS and LightSampler work, implement pbrt's RandomWalkIntegrator using the same Film, Camera, Material, BSDF, and Light contracts. Its purpose is to validate:

- wavelength generation and sensor conversion;
- emission evaluation;
- BSDF sampling and throughput;
- closure lifetime;
- nested-region state transitions;
- deterministic sample streams.

Do not ship it as the recommended production integrator.

## 11.10 VolPathIntegrator

After the Medium PR, port pbrt's volumetric path logic rather than inserting medium samples into the surface integrator ad hoc. It MUST preserve:

- majorant transmittance sampling;
- real, absorption, and null collision handling;
- spectral throughput and majorant ratios;
- phase-function direct-light sampling and MIS;
- surface BSDF behavior;
- wavelength packet handling;
- emitter-hit MIS;
- nested dielectric region state across null boundaries.


# 12. Nested and overlapping dielectric regions

This section is normative. It preserves rayrender's lower-number-wins dielectric priority behavior while moving scattering to pbrt-style closures.

## 12.1 Design objective

A region-bearing closed primitive defines a volume with:

- a stable region identity;
- an index-of-refraction spectrum;
- an optional Medium for absorption, scattering, and emission;
- a numeric priority used only when multiple regions overlap.

The surface Material defines boundary roughness and scattering model. It does not define path membership.

At a boundary hit, the region resolver determines the effective optical region on the geometric outside and inside of that boundary. A dielectric Material uses the resulting relative eta to construct a normal pbrt `DielectricBxDF`. The integrator updates membership only if the sampled event crosses the boundary.

## 12.2 Region identity is not material identity

Do not store Material pointers in the region state. The same Material may be reused by multiple objects, and each object still represents a distinct spatial region.

```cpp
using RegionId = uint32_t;

struct DielectricRegion {
    RegionId id;
    int priority;
    Spectrum eta;
    Medium medium;
    std::string_view debugName;
};
```

Rules:

- each closed object instance receives a unique RegionId by default;
- all triangles belonging to one watertight mesh instance share one RegionId;
- transformed instances receive distinct RegionIds unless an advanced descriptor explicitly declares them to be one connected region;
- a region descriptor is immutable after scene compilation;
- RegionId zero is reserved for the exterior region;
- material sharing never causes region sharing implicitly.

The exterior descriptor is conceptually:

```cpp
DielectricRegion exterior {
    .id = 0,
    .priority = std::numeric_limits<int>::max(),
    .eta = ConstantSpectrum(1),
    .medium = nullptr
};
```

A future scene-level exterior medium or IOR MAY replace these defaults, but it must remain an explicit root descriptor.

## 12.3 Path-local membership state

Although the user-facing concept remains a dielectric stack, overlaps mean the implementation is a membership set, not necessarily a LIFO stack.

```cpp
class DielectricPathState {
public:
    bool Contains(RegionId id) const;
    RegionHandle Active() const;
    ResolvedDielectricTransition Analyze(
        RegionHandle crossed,
        bool entering) const;
    void Commit(const DielectricTransitionToken& token);
    Medium ActiveMedium() const;
    uint64_t Generation() const;

private:
    InlinedVector<RegionHandle, 4> members_;
    uint64_t generation_ = 0;
};
```

Start with a simple value-semantic container. Optimize only after profiling. Copying this state must produce independent membership. No raw pointer aliases mutable state between path branches or visibility rays.

`Active()` returns the contained region with the smallest numeric priority. Empty membership returns the exterior region.

For deterministic diagnostics, compare by `(priority, RegionId)`, but equal-priority overlapping regions are invalid by default. The deterministic RegionId tie-break exists only to avoid undefined behavior while reporting the error.

## 12.4 Geometric side states

For a boundary of region `C`, define two hypothetical membership states independent of the ray direction:

```text
outside state: current memberships with C absent
inside state:  current memberships with C present
```

Then:

```text
activeOutside = active(outside state)
activeInside  = active(inside state)
```

This construction is critical. It provides the eta values on the two geometric sides of the primitive normal, which is exactly what a pbrt `DielectricBxDF` expects.

Entry/exit uses only the geometric normal:

```cpp
bool entering = Dot(ray.d, si.n) < 0;
```

where `si.n` is oriented outward from the crossed region. Shading and bump normals do not participate.

For an entering ray:

```text
activeBefore = activeOutside
activeAfter  = activeInside
```

For an exiting ray:

```text
activeBefore = activeInside
activeAfter  = activeOutside
```

The resolver MUST verify that `activeBefore` matches the current path state's active region. A mismatch is a state-consistency diagnostic, not a reason to guess silently.

## 12.5 Transition object

```cpp
struct ResolvedIORRatio {
    Spectrum etaOutside;
    Spectrum etaInside;

    Float Evaluate(Float lambdaNm) const;
    bool IsConstant() const;
    bool IsUnity() const;
};

enum class DielectricBoundaryKind {
    PrioritySkipped,
    IndexMatchedNull,
    ScatteringInterface
};

struct DielectricTransitionToken {
    RegionId crossed;
    bool entering;
    uint64_t expectedGeneration;
};

struct ResolvedDielectricTransition {
    bool entering;
    DielectricBoundaryKind kind;

    RegionHandle crossed;
    RegionHandle activeOutside;
    RegionHandle activeInside;
    RegionHandle activeBefore;
    RegionHandle activeAfter;

    ResolvedIORRatio etaRatio;
    Medium mediumBefore;
    Medium mediumAfter;

    DielectricTransitionToken token;
};
```

Classification:

1. `PrioritySkipped`: `activeOutside == activeInside`. Crossing the object's boundary changes membership but not the effective optical region.
2. `IndexMatchedNull`: effective active region changes, but the relative IOR is identically one and the boundary Material has no separate reflective coating or rough-surface effect. The path crosses straight through and the Medium may change.
3. `ScatteringInterface`: the effective IOR ratio is not identically one and the pbrt dielectric closure must be evaluated.

A boundary with an explicit coated or otherwise non-dielectric Material is not automatically null merely because bulk eta is equal. The SceneCompiler must reject unsupported combinations rather than apply region logic ambiguously.

## 12.6 Constant-ratio detection

pbrt's dielectric Material terminates secondary wavelengths when its eta spectrum is nonconstant. For two effective regions, the relevant spectrum is their ratio:

```text
etaRatio(lambda) = etaInside(lambda) / etaOutside(lambda)
```

`ResolvedIORRatio::IsConstant()` MUST be conservative and deterministic:

- true if both eta spectra are constant;
- true if both handles refer to the same immutable spectrum;
- true for analytically recognized scaled forms with constant ratio;
- otherwise false unless a representation can prove constancy.

Do not decide constancy by testing only the four current wavelengths. That would make path topology depend unstably on the sampled packet.

`IsUnity()` follows the same proof-based rule. A numerically near-unity ratio may still use the dielectric closure; do not introduce an arbitrary physical threshold into transport. Exact index matching can use the null path.

## 12.7 Closure construction

For a scattering interface, the Material closure receives:

```cpp
struct ResolvedDielectricInterface {
    Spectrum etaOutside;
    Spectrum etaInside;
    bool ratioIsConstant;
    DielectricTransitionToken token;
};
```

The dielectric Material performs:

```cpp
if (!ctx.dielectric->ratioIsConstant) {
    lambda.TerminateSecondary();
}

Float etaOutside = ctx.dielectric->etaOutside(lambda[0]);
Float etaInside = ctx.dielectric->etaInside(lambda[0]);
Float eta = etaInside / etaOutside;

return DielectricBxDF(
    eta,
    TrowbridgeReitzDistribution(alphaU, alphaV));
```

This follows pbrt's convention: the BxDF's scalar `eta` is the IOR on the positive-normal inside side divided by the IOR on the outside side. The BxDF itself uses the sign of local `wo.z` to handle entry versus exit.

Secondary wavelengths MUST be terminated before the Fresnel branch is sampled, even if reflection is later selected. This deliberately follows pbrt-v4's behavior rather than applying a rayrender-specific optimization.

Do not terminate wavelengths for `PrioritySkipped` or `IndexMatchedNull` crossings because no wavelength-dependent direction is selected.

## 12.8 Transactional state updates

`Analyze()` is pure. It MUST NOT mutate the path state.

`Commit(token)` applies exactly one membership change:

```text
entering: add crossed region
exiting:  remove crossed region
```

It verifies:

- the state generation matches `expectedGeneration`;
- an entering region is not already present;
- an exiting region is present;
- the RegionId is valid.

After commit, increment the generation and assert that `Active()` equals the transition's recorded `activeAfter`.

The integrator commits according to event type:

| Event | Commit membership? | Update ray Medium? | Update `etaScale`? |
|---|---:|---:|---:|
| reflected dielectric sample | no | no | no |
| transmitted dielectric sample | yes | yes, to active-after Medium | yes, from `BSDFSample::eta` |
| total internal reflection | no | no | no |
| priority-skipped crossing | yes | only if active Medium changed, normally no | no |
| index-matched null crossing | yes | yes | no |
| alpha-rejected hit | no | no | no |
| thin dielectric reflection/transmission | no | no | use pbrt thin-sheet behavior only |

The closure never sees `DielectricPathState` and cannot commit it.

## 12.9 Null/skipped traversal

For `PrioritySkipped` and `IndexMatchedNull`:

1. commit the transition;
2. update `ray.medium` to the new active Medium;
3. spawn a ray in the unchanged direction using a robust offset across the geometric boundary;
4. do not increment scattering depth;
5. do not consume BSDF random dimensions;
6. do not terminate wavelengths;
7. do not set `specularBounce`;
8. do not overwrite the previous BSDF PDF or previous light context;
9. do not modify `etaScale`;
10. continue intersection traversal.

This mirrors pbrt's treatment of null material boundaries while retaining rayrender's priority semantics.

To avoid an infinite loop at a null boundary, the spawned ray origin MUST be offset to the transmitted side using the geometric normal and outgoing direction. Add a bounded null-boundary counter and fail diagnostically if a path traverses an implausibly large number without making progress.

## 12.10 Rough boundaries

Rough dielectric reflection and transmission use the same transaction rule as smooth glass:

- `BSDFSample::IsReflection()` leaves membership unchanged;
- `BSDFSample::IsTransmission()` commits the crossing;
- an invalid sample or TIR leaves membership unchanged;
- `bs.eta` and the BxDF's radiance-mode scaling are used exactly as in pbrt.

Do not infer crossing from the sign of the generated direction alone when a valid `BSDFSample` already reports reflection/transmission flags.

## 12.11 Segment media and absorption

The active region before an intersection owns the Medium through which the segment traveled. Therefore:

```text
segment medium = regions.ActiveMedium() at ray launch
```

This must agree with `ray.medium`. Add debug assertions after every committed transition.

Bulk Beer-Lambert attenuation is implemented by an absorption Medium:

```text
Tr(lambda) = exp(-sigma_a(lambda) * distance)
```

It is not applied in `DielectricMaterial::GetBxDF()` or when the ending surface is encountered.

Consequences:

- entering a region does not retroactively apply its absorption to the preceding segment;
- exiting a region applies its Medium over the segment before the exit through normal medium integration;
- skipped boundaries preserve the Medium if the active region does not change;
- active index-matched boundaries may change Medium while keeping direction unchanged;
- environment paths in a non-vacuum exterior Medium follow VolPath's medium handling.

For strict pbrt conformity, scenes with non-null media use `integrator_type = "volpath"`. The surface PathIntegrator MUST not contain a second ad hoc attenuation system. Before VolPath is complete, spectral compilation rejects nonzero region `sigma_a` or instructs the user to use an explicitly marked experimental absorption-only mode that is not part of acceptance.

## 12.12 Visibility and direct-light segments

A visibility query starts with a value copy of the current `DielectricPathState` and the current ray Medium. It MUST NOT mutate the camera path's state.

Traversal behavior:

- alpha-rejected surfaces are ignored without state changes;
- priority-skipped and index-matched null boundaries are committed locally and traversal continues;
- active dielectric scattering boundaries block the straight visibility segment;
- opaque surfaces block the segment;
- media contribute transmittance in VolPath;
- the sampled light endpoint is handled without crossing its boundary spuriously.

This is the appropriate estimator for an ordinary path tracer: it does not attempt to connect a shading point to a light through a refractive interface along a bent path. Such transport is reached through BSDF sampling.

The visibility routine SHOULD return spectral transmittance rather than a Boolean once media are implemented:

```cpp
SampledSpectrum Tr(
    const Interaction& p0,
    const Interaction& p1,
    const SampledWavelengths& lambda,
    Sampler& sampler,
    DielectricPathState regions) const;
```

The region state is passed by value intentionally.

## 12.13 Camera-origin initialization

A camera ray may begin inside one or more dielectric regions, so the initial membership cannot always be empty.

Add two mechanisms:

1. **Explicit initialization**, required for complete control:

```r
camera(..., initial_regions = c("outer_glass", "inner_liquid"))
```

2. **Automatic containment**, used when every relevant shape supplies a reliable point-containment query:

```r
camera(..., initial_regions = "auto")
```

The SceneCompiler assigns stable user-facing region names and validates explicit IDs. Automatic initialization evaluates the generated ray origin and time
against each region's complete boundary assembly or explicit containment
predicate.

Rules:

- analytic closed shapes MUST implement robust containment;
- watertight triangle meshes MAY use an acceleration structure or winding-number method;
- non-watertight meshes cannot participate in automatic initialization;
- animated camera or region transforms require initialization at the camera ray's sampled time;
- a static camera in static regions MAY cache the result;
- explicit initialization may be validated against geometry in strict/debug mode;
- ambiguity or unsupported containment produces an error, not an empty-stack fallback.

This corrects the legacy limitation where a camera placed inside a dielectric did not fully inherit surrounding priority state.


## 12.14 Region boundary assemblies and topology requirements

A `RegionId` identifies one optical volume instance, not one `Shape`, `Mesh`, or `Primitive`.

A region MAY be bounded by an oriented assembly of multiple primitives, including primitives from different meshes. An individual boundary primitive MAY be open. When closed-region validation is requested, the closure requirement applies to the oriented union of all boundary patches associated with the region.

Material assignment and region-boundary assignment MUST be independent. A primitive MAY:

* have a Material without defining a region transition;
* define a region boundary without having a Material;
* have both a Material and one or more region-boundary attachments;
* participate in the boundaries of multiple regions at a shared interface.

Each region-boundary attachment MUST specify which side of the final geometric normal contains the region:

```cpp
enum class RegionSide {
    NegativeNormal,
    PositiveNormal
};

struct RegionBoundaryAttachment {
    RegionHandle region;
    RegionSide containedSide;
};
```

The renderer MUST determine membership before and after an interaction from the incoming and outgoing geometric sides. It MUST NOT infer region membership solely from whether a mesh is nominally entered or exited.

Multiple primitives MAY reference the same `RegionId`. Such repeated references are required for boundaries assembled from multiple meshes. Duplicate region **definitions** with incompatible eta, Medium, priority, or identity properties are errors.

A visible primitive MAY serve simultaneously as shading geometry and as a region boundary. Do not add coincident duplicate geometry solely to close a region. In particular, an opaque or reflective seafloor primitive MAY form the lower boundary of a water region while retaining its terrain Material.

For a reflected surface event, region membership remains on the incident side. For a transmitted or null-boundary event, region membership is updated to match the outgoing side. If one primitive separates multiple tracked regions, all membership changes MUST be committed atomically.

For strict closed-region validation, the compiler SHOULD validate the complete boundary assembly after applying transforms and virtual seam welding. Diagnostics should detect when practical:

* unmatched boundary edges across the complete assembly;
* gaps between component meshes;
* inconsistent region-side declarations;
* inconsistent orientation;
* finite-area overlap between boundary patches;
* zero-volume boundary assemblies;
* self-intersections and cross-patch intersections;
* conflicting region definitions;
* an eta spectrum with nonpositive values;
* equal-priority overlapping active regions encountered at runtime.

Shared edges or seams between separate meshes are valid. Coincident surfaces overlapping over a finite area are not valid substitutes for a shared boundary.

Opaque primitives that cannot be crossed MAY omit explicit region-boundary metadata in pbrt-compatible transport mode. In that case, the incoming ray region is retained for reflected rays. However, such omitted boundaries cannot be used to prove region closure or perform automatic point containment.

Automatic camera-region initialization requires either:

* a closed analytic region;
* a closed and consistently oriented boundary assembly;
* an explicitly supported constructive or implicit volume representation; or
* explicit user-provided initial region membership.

A region with reachable open boundaries and differing optical state beyond those boundaries is outside correctness guarantees. The compiler SHOULD reject it in strict mode and issue a diagnostic in advanced mode.

## 12.15 Interaction with instancing and transforms

Region membership follows object instances, not source meshes. The intersection adapter must return the instance's RegionHandle after transformations.

Negative-determinant transforms can reverse orientation. The shape/instance layer MUST ensure that the final geometric normal still points toward the declared outside and that entry/exit tests remain correct.

Animated transforms require the same RegionId over time. A path cannot transfer membership between two different instances merely because they share source geometry.

Multiple component primitives within one region instance share one RegionId.
Instancing the complete assembly creates a new RegionId for that volume
instance.

## 12.16 Coincident and nearly coincident boundaries

Coincident active dielectric boundaries are order-sensitive and are not solved automatically by scalar priority alone. Establish the following policy:

- exact equal-priority overlaps are an error;
- coincident boundaries with distinct priorities emit a scene diagnostic;
- traversal is deterministic using aggregate intersection ordering and RegionId tie-breaks;
- acceptance scenes avoid mathematically coincident active interfaces unless a dedicated grouped-boundary algorithm is later implemented;
- a future grouped-boundary extension must collect all hits within a robust `t` interval and resolve one combined before/after state before constructing a closure.

Do not claim CSG robustness beyond the tested overlap cases.

A shared physical interface MUST normally be represented by one geometric
primitive with all relevant Material and region-side metadata. Do not create
one coincident primitive per region.

## 12.17 State-error policy

Add counters and optional path traces for:

- duplicate entry;
- exit from an absent region;
- active-before mismatch;
- equal-priority active overlap;
- commit generation mismatch;
- ray Medium disagreement;
- excessive null-boundary traversal;
- nonfinite eta ratio.

In development and test builds, these are assertions or hard errors. In release rendering, invalidate the affected sample, increment a counter, and report a summarized warning after rendering. Do not continue with an invented state.

## 12.18 Required nested-region unit tests

At minimum, test state transitions without rendering:

1. exterior -> glass -> exterior;
2. exterior -> outer glass -> inner liquid -> outer glass -> exterior;
3. enter lower-priority region while higher-priority region is active;
4. exit lower-priority region while higher-priority region remains active;
5. enter higher-priority region while lower-priority region is active;
6. exit higher-priority region to reveal lower-priority region;
7. three-way overlap with deterministic active selection;
8. reflection at every active boundary leaves state unchanged;
9. transmission commits exactly once;
10. TIR leaves state unchanged;
11. index-matched active transition changes Medium but not direction/depth;
12. skipped dispersive region does not terminate wavelengths;
13. active dispersive transition terminates wavelengths before sampling;
14. rough reflection does not commit, rough transmission does;
15. alpha-rejected dielectric hit makes no transition;
16. copied visibility state mutates independently;
17. explicit camera initialization;
18. automatic camera initialization;
19. equal-priority overlap diagnostic;
20. RegionId differs for two instances sharing one Material.

## 12.19 Required nested-region render tests

Create deterministic high-sample reference scenes:

- concentric air/glass/water spheres;
- a glass shell with an air bubble represented by priority regions;
- two partially overlapping dielectrics where the lower numeric priority dominates;
- the same scene with region declaration order reversed;
- nested constant-IOR objects compared with a scene using explicit nonoverlapping interfaces;
- a prism with a nested nondispersive inclusion;
- a dispersive outer glass and nondispersive inner liquid;
- a dispersive lower-priority object fully hidden by a higher-priority region;
- absorbing glass around a clear bubble under VolPath;
- camera inside one region and inside two overlapping regions;
- rough nested interfaces;
- a direct-light visibility segment crossing only skipped/null boundaries;
- a segment crossing an active glass interface, which must be treated as blocked for straight shadow visibility.

For symmetric scenes, include symmetry and region-declaration-order checks in addition to image RMSE.

## 12.20 Reference transition pseudocode

Codex should implement the equivalent of the following control flow:

```cpp
for (int depth = 0; ; ) {
    std::optional<ShapeIntersection> hit = Intersect(ray);

    if (!hit) {
        AccumulateInfiniteLights(...);
        break;
    }

    SurfaceInteraction& si = hit->intr;

    if (!AcceptAlpha(si, sampler)) {
        ray = si.SpawnRay(ray.d);
        continue;
    }

    std::optional<ResolvedDielectricTransition> transition;
    if (si.dielectricRegion) {
        bool entering = Dot(ray.d, si.n) < 0;
        transition = regions.Analyze(si.dielectricRegion, entering);
    }

    if (transition &&
        transition->kind != DielectricBoundaryKind::ScatteringInterface) {
        regions.Commit(transition->token);
        ray = si.SpawnRay(ray.d, regions.ActiveMedium());
        continue;
    }

    MaterialEvalContext ctx(si);
    ResolvedDielectricInterface resolved;
    if (transition) {
        resolved = MakeResolvedInterface(*transition);
        ctx.dielectric = &resolved;
    }

    BSDF bsdf = si.material.GetBSDF(
        textureEvaluator, ctx, lambda, scratch);

    // Direct-light sampling occurs here for eligible BSDFs.

    std::optional<BSDFSample> bs = bsdf.Sample_f(...);
    if (!bs) {
        break;
    }

    beta *= bs->f * AbsDot(bs->wi, si.shading.n) / bs->pdf;

    if (transition && bs->IsTransmission()) {
        regions.Commit(transition->token);
        etaScale *= Sqr(bs->eta);
    }

    Medium nextMedium = regions.ActiveMedium();
    ray = si.SpawnRay(bs->wi, nextMedium);
    ++depth;
}
```

The final implementation must also include emission, MIS, medium sampling, regularization, robust offsets, roulette, and diagnostics. The pseudocode only illustrates state ownership and commit timing.



## 12.21 CSG interaction with the spectral architecture

CSG is not part of pbrt, but it is a supported rayrender extension. The implementation MUST treat CSG as a Shape implementation with explicit capability flags, not as a special material or integrator case.

## 12.21.1 CSG shape capabilities

Add a Shape capability record used by validation and compilation:

```cpp
struct ShapeCapabilities {
    bool canIntersect = false;
    bool canSampleArea = false;
    bool hasExactArea = false;
    bool hasImplicitInterior = false;
    bool hasUV = false;
    bool supportsDisplacement = false;
    bool supportsNormalMap = false;
};
```

Initial CSG capabilities:

```text
canIntersect = true
canSampleArea = false unless tessellated or otherwise sampleable
hasExactArea = false unless tessellated or analytically known
hasImplicitInterior = true only for closed/coherent CSG fields
hasUV = false by default
supportsDisplacement = false by default
supportsNormalMap = false by default
```

The CSG adapter MUST expose bounds, intersection, normal, and, when valid, `Contains(p)` and `SignedDistance(p)`. Existing `getDistance()` logic may remain internally, but spectral region code should depend on a named signed-distance/containment interface.

```cpp
struct ImplicitEval {
    Float phi;
    Vector3f grad;
    bool validGrad;
};

class ImplicitShape {
public:
    virtual Float SignedDistance(const Point3f& p) const = 0;
    virtual ImplicitEval Eval(const Point3f& p) const;
    virtual bool Contains(const Point3f& p) const;
    virtual bool HasCoherentInterior() const = 0;
    virtual CSGTopologyFlags TopologyFlags() const = 0;
};
```

## 12.21.2 CSG as ordinary geometry

CSG with an ordinary surface Material is supported through the same pipeline as any other shape:

```text
CSGShape -> GeometricPrimitive -> Material -> Material::GetBSDF() -> BxDF closure
```

Constant materials and 3D procedural textures are allowed. UV image textures, normal maps, and displacement are rejected unless a generated mapping such as triplanar mapping is provided. Current CSG hit records that set placeholder UVs are not sufficient for spectral image-texture correctness.

## 12.21.3 CSG as an optical-region boundary

A closed CSG object may define a region boundary using:

```r
with_region_boundary(region_boundary(region = "glass_blob", side = "inside"))
```

For CSG, `inside` means `SignedDistance(p) < 0` after all object transforms. The region resolver determines membership before and after a transmission or null crossing using the signed field or robust offsets across the surface. This applies equally to dielectric regions and pure medium/null boundaries.

Open or ambiguous CSG expressions MUST NOT be accepted as optical-region boundaries in strict mode. Examples that are invalid by default include infinite planes, isolated triangles, and shape mixes that do not define Boolean occupancy.

Default topology expectations:

| CSG form | Region-boundary default |
|---|---:|
| sphere, box, torus, capsule, cylinder, ellipsoid | valid |
| cone, rounded cone, pyramid | valid if nondegenerate |
| plane, triangle | invalid |
| union/group | valid if all children are valid |
| subtract | valid if operands are closed and the result has coherent occupancy |
| intersection | valid if all children are valid |
| blend/subtractblend | valid with warning unless topology can be proven |
| mix | invalid for optical regions initially |
| onion | valid shell if child is valid and thickness is safely positive |
| round, elongate, translate, rotate, scale | inherit child validity, subject to nondegenerate transforms |

## 12.21.4 CSG as an area light

CSG does not satisfy pbrt's area-light contract unless it can provide area, surface sampling, and a PDF consistent with the rendered emitting surface. Therefore:

- CSG area lights are rejected in strict spectral mode unless the CSG shape is tessellated or otherwise made sampleable.
- `csg_mesh_sampler(mode = "render_mesh")` is the preferred strict-mode route because the rendered surface and sampled emitting surface are the same.
- `sampling_proxy` mode is advanced only and is outside conformance guarantees unless a future estimator proof and tests are added.
- Unsampled CSG emission may be allowed only in `scene_validation("advanced", allow_unsampled_emitters = TRUE)` and must be documented as noisy/non-pbrt-like.

`with_light()` remains generic: it attaches the descriptor. The SceneCompiler decides whether the decorated shape satisfies the Light contract.

## 12.21.5 Per-child CSG materials and regions

The first spectral implementation supports one Material, zero or one AreaLight, and zero or more region-boundary attachments per CSG object. It does not support per-child CSG materials or per-child optical regions.

Adding per-child support would require the CSG evaluator to return feature/material provenance in addition to distance and gradient:

```cpp
struct CSGEval {
    Float phi;
    Vector3f grad;
    int materialSlot;
    int regionSlot;
    int featureId;
};
```

Boolean provenance is tractable for sharp union/intersection/subtraction, but smooth blends make material and region provenance ambiguous. Do not include this in PRs 0-24.

## 12.21.6 CSG numerical robustness

CSG region boundaries are sensitive to missed intersections. Strict spectral mode requires additional diagnostics for region-bearing CSG:

- scale-aware intersection epsilon;
- root refinement after a near hit;
- warnings for shells thinner than the CSG marching threshold;
- tests for tangent rays, subtractive cavities, smooth blends, onion shells, and nonuniform transforms;
- a bounded null-boundary counter to catch no-progress loops;
- validation that `Contains()` changes consistently across a committed crossing.

A missed CSG boundary corrupts region membership and therefore eta, Medium, and path throughput. Treat such failures as correctness diagnostics, not merely visual artifacts.


# 13. Step-by-step pull request sequence

Each PR below is an independently reviewable changeset. The listed gate is mandatory. Codex must include a short implementation report containing changed files, tests run, benchmark delta, pbrt references used, and open deviations.

## PR 0: Pin baselines, create conformance records, and capture legacy output

### Goal

Freeze the reference points and establish reproducible baselines before architectural work.

### Required changes

1. Record the exact rayrender commit/source checksum and pbrt-v4 commit in `docs/spectral/versions.md`.
2. Add `docs/spectral/pbrt-conformance.md` with columns for subsystem, pbrt file/function, rayrender implementation, status, tests, and deviation.
3. Add `docs/spectral/adr/0001-spectral-renderer.md` containing the decisions in Sections 2 and 3.
4. Copy this plan into `docs/spectral/rayrender_spectral_rendering_codex_plan.md`.
5. Inventory copied/generated data licenses and expected notices in `docs/spectral/provenance.md`.
6. Add a developer test command that installs with `tools/codex/install-local.sh`, runs unit tests, and renders small regression scenes.
7. Capture legacy seeded outputs, hashes, compiler information, and render options for representative still and animation frames.
8. Capture CPU time and peak-memory baselines for at least three current scenes.
9. Add a CI job or local script that verifies the source package contains all generated/packaged spectral assets later added.

### Gate

- Existing package tests pass.
- Baseline scenes reproduce from a clean checkout.
- No production renderer code changes.
- The pinned pbrt commit is immutable in project documentation.

## PR 1: Generalize the build and extract RenderSession

### Goal

Create one scene/render/output setup path shared by still and animation without altering rendering mathematics.

### Required changes

1. Make `tools/config/configure.R` discover new C++ source directories recursively or consume one maintained manifest.
2. Add `RenderOptions`, `CompiledScene`, and `RenderSession` placeholders under new source directories.
3. Extract duplicated still/animation setup into shared functions:
   - scene compilation;
   - environment setup;
   - camera construction;
   - sampler construction;
   - output buffers;
   - progress and cancellation;
   - postprocessing and file output.
4. Keep all calls routed to the existing RGB integrators.
5. Make animation frame changes explicit inputs to `RenderSession::RenderFrame()` rather than rebuilding unrelated state.
6. Preserve Rcpp exception translation and thread interruption behavior.
7. Add source-layout documentation and prohibit new top-level renderer setup duplication.

### Gate

- Bitwise or existing-tolerance parity for legacy still renders.
- Bitwise or existing-tolerance parity for legacy animation frames.
- Package builds on all supported compilers.
- No spectral types introduced yet.

## PR 2: Add explicit base types, dispatch handles, and ScratchBuffer

### Goal

Build the low-level foundation needed by pbrt-style closures without changing rendered output.

### Required changes

1. Add explicit `RGB`, `XYZ`, color-encoding, and color-space matrix types.
2. Add safe componentwise math, finite checks, and matrix tests.
3. Add a tagged-handle strategy for new `Spectrum`, `Texture`, `Material`, `BxDF`, `Light`, `Medium`, and `PhaseFunction` types.
4. Add `ScratchBuffer` with alignment, reset, high-water instrumentation, and unit tests.
5. Add bitmask helpers for future BxDF flags.
6. Add `std::optional`-based sampling result conventions.
7. Keep these types unused by legacy transport except for isolated conversion tests.
8. Add compiler warnings or lint rules where feasible to prevent accidental implicit `RGB`/geometry conversion.

### Gate

- Legacy render baselines unchanged.
- Scratch allocations meet alignment requirements for every registered closure type.
- Sanitizer test demonstrates no use after reset and no leaks.
- Tagged dispatch tests cover every registered dummy type.

## PR 3: Implement SampledSpectrum and SampledWavelengths

### Goal

Port pbrt's fixed-width spectral packet and wavelength sampling exactly.

### Required changes

1. Add constants for four samples and the pinned visible interval.
2. Implement `SampledSpectrum` arithmetic and diagnostics.
3. Implement pbrt visible wavelength PDF and inverse sampling.
4. Implement four correlated wavelength generation.
5. Implement uniform sampling for tests only.
6. Implement idempotent `TerminateSecondary()`.
7. Add debug assertions for wavelength/PDF invariants.
8. Add serialization/debug printing only outside hot paths.
9. Update `pbrt-conformance.md` with exact referenced pbrt functions.

### Gate

- Dense numerical comparison to pbrt reference values.
- Monte Carlo histogram agrees with the intended visible PDF.
- Packet correlations match the shifted-sample construction.
- Termination tests verify zero secondary PDFs and first-PDF scaling.
- No renderer behavior change.

## PR 4: Implement Spectrum representations and named spectral assets

### Goal

Support sampled and analytic scene spectra independently of RGB reconstruction.

### Required changes

1. Implement `ConstantSpectrum`, `PiecewiseLinearSpectrum`, `DenselySampledSpectrum`, and `BlackbodySpectrum`.
2. Add interpolation, extrapolation, normalization, and validation policies.
3. Package CIE x/y/z and required illuminant data with provenance and checksums.
4. Add a named-spectrum registry.
5. Add spectrum integration helpers used only for setup, testing, sensor construction, and light power estimation.
6. Add R-independent C++ parsing for packaged binary/text assets.
7. Add deterministic scripts that regenerate the assets.
8. Add unit-aware comments and API validation; internal wavelength unit is nanometers.

### Gate

- CIE and illuminant samples match source data.
- Blackbody values match a high-precision implementation up to a common normalization.
- Interpolation and extrapolation edge cases are covered.
- Asset checksums are verified during tests.
- Installed-package lookup works without a source tree.

## PR 5: Add RGB color spaces and RGB-to-spectrum table assets

### Goal

Implement pbrt's color-space and RGB reconstruction machinery, beginning with sRGB.

### Required changes

1. Implement chromaticity conversion and RGB/XYZ matrices.
2. Implement linear and sRGB transfer functions explicitly.
3. Define the sRGB color space and its standard illuminant.
4. Add a versioned binary `RGBToSpectrumTable` format.
5. Add deterministic table-generation/import scripts tied to the pinned pbrt algorithm.
6. Add runtime loading, validation, one-time caching, and failure diagnostics.
7. Implement sigmoid polynomial evaluation and table interpolation.
8. Implement `RGBAlbedoSpectrum`, `RGBUnboundedSpectrum`, and `RGBIlluminantSpectrum` exactly.
9. Record table license/provenance and package-size impact.
10. Defer other color spaces until PR 23 unless they are needed for table-validation infrastructure.

### Gate

- RGB reconstruction coefficients match the pinned pbrt table within recorded tolerance.
- Dense RGB cube samples reconstruct to expected RGB under the reference illuminant.
- Albedo reconstruction remains bounded according to pbrt semantics.
- white illuminant reconstruction includes the standard illuminant rather than becoming equal-energy white.
- corrupt, missing, or wrong-version table files fail clearly.

## PR 6: Add schema-v2 spectral descriptors and legacy adaptation shell

### Goal

Expose validated descriptors and the new functional public API surface in R without rendering them yet.

### Required changes

1. Add the spectrum constructors in Section 9.4.
2. Add render-option descriptors: integrator, sampler, camera, film, sensor, spectral options, and scene validation.
3. Add schema-v2 descriptor classes and print methods for spectra, textures, materials, lights, media, regions, region boundaries, cameras, films, sensors, samplers, integrators, and validation.
4. Add `ray_scene_v2` as a tibble-compatible hybrid scene with row-level list columns and scene-level registries described in Section 9.13.
5. Add `ray_material_v2` as a pbrt-style named parameter dictionary described in Section 9.14.
6. Add `with_light()`, `with_region_boundary()`, `add_light()`, `set_environment()`, `add_region()`, named texture/material registries, and reference descriptors.
7. Add strict validation and stable error paths, including CSG capability validation hooks.
8. Add a centralized schema-v1 adapter with a conversion-report object.
9. Add warning aggregation so legacy conversions warn once per scene, not per hit.
10. Preserve all current R constructors and signatures.
11. Add serialization tests between R and the C++ SceneCompiler input layer.
12. Document RGB role selection, input color encoding, physical light units, and the difference between Material, Light, Medium, and optical Region.

### Gate

- [x] Existing R scene construction tests pass.
- [x] Every new descriptor round-trips through serialization.
- [x] invalid spectra, lights, regions, and CSG capability combinations fail before rendering.
- [x] schema-v1 adaptation is deterministic and produces no worker-thread warnings.
- [x] `ray_scene_v2` remains compatible with existing object-row workflows.
- [x] `with_light()` and `with_region_boundary()` preserve row counts and scene attributes.
- [x] R code uses `=` assignment throughout.

### Completion record

Completed in PR 6 by `R/schema_v2_descriptors.R`, `R/render_scene.R`,
`R/add_object.R`, `R/ray_scene.R`, `tests/testthat/test-schema-v2-descriptors.R`,
and generated export/documentation updates.

Gate evidence:

- `testthat::test_file("tests/testthat/test-schema-v2-descriptors.R")`
- `tools/codex/install-local.sh`

## PR 7: Add FloatTexture, SpectrumTexture, and image color management

### Goal

Port pbrt's texture evaluation split and make image semantics explicit.

### Required changes

1. Add `TextureEvalContext`, `FloatTexture`, `SpectrumTexture`, and tagged dispatch.
2. Implement constant, scale, mix, checker, procedural, and image forms required by existing rayrender scenes.
3. Add `UniversalTextureEvaluator`.
4. Decode encoded color images to linear RGB once.
5. Filter linear texels, then reconstruct spectra using the consuming texture's role.
6. Keep scalar/data maps linear and never route them through RGB reconstruction.
7. Add image cache keys that include filename, encoding, color space, semantic role, wrap mode, and filtering options.
8. Retain legacy texture classes only for the legacy renderer.
9. Add derivative and filtering adapters from existing hit records.
10. Add spectral texture debug evaluation utilities outside hot paths.

### Gate

- [x] Scalar images are unaffected by gamma settings.
- [x] sRGB color textures decode/filter/reconstruct in the specified order.
- [x] image cache does not alias different semantic interpretations.
- [x] analytic texture tests match expected scalar and spectral values.
- [x] no per-evaluation heap allocation.

### Completion record

Completed in PR 7 by `src/materials/spectral_texture.h`,
`src/materials/spectral_texture.cpp`,
`src/materials/spectral_texture_adapters.cpp`,
`src/materials/spectral_texture_test.cpp`,
`tools/spectral-tests/pr7-texture-tests.cpp`, and
`tools/spectral-tests/run-pr7-texture-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr7-texture-tests.R`
- `tools/codex/install-local.sh`
- focused schema-v2 descriptor test

## PR 8: Add PixelSensor and Film

### Goal

Establish the correct spectral-to-tristimulus estimator before tracing spectral paths.

### Required changes

1. Implement CIE-based and pbrt-compatible RGB PixelSensor modes.
2. Implement sensor-response integration with wavelength-PDF division.
3. Add Film sample accumulation, reconstruction filter support, and tile merging.
4. Store sensor channels plus sample weights in Film.
5. Add output-color conversion and encoding after Film accumulation.
6. Add visible-surface storage hooks for future AOVs.
7. Add deterministic single-thread mode for tests.
8. Add sensor calibration/normalization matching the pinned pbrt source.
9. Ensure preview conversion reads Film channels rather than packet components.
10. Keep legacy `RayMatrix` output path intact for `rgb_legacy`.

### Gate

- [x] Monte Carlo integration of known spectra agrees with deterministic CIE integration.
- [x] results are invariant, within statistical tolerance, under visible versus uniform wavelength sampling.
- [x] wavelength-PDF division happens once; a deliberate double-division test fails.
- [x] Film tile and single-thread accumulation agree.
- [x] white-point and output-space tests pass.

### Completion record

Completed in PR 8 by `src/render/spectral_film.h`,
`src/render/spectral_film.cpp`, `tools/spectral-tests/pr8-film-tests.cpp`,
and `tools/spectral-tests/run-pr8-film-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr8-film-tests.R`
- `tools/codex/install-local.sh`


## PR 9: Add wavelength-aware camera interfaces

### Goal

Move camera sampling onto the pbrt contract before implementing the spectral integrator.

### Required changes

1. Add `CameraSample`, `CameraRay`, and a camera handle used only by spectral mode.
2. Port/adapt orthographic, perspective, and thin-lens cameras to `GenerateRay(sample, lambda)`.
3. Return a `SampledSpectrum` camera weight.
4. Set the initial ray Medium.
5. Preserve ray differentials where currently supported.
6. Add stubs for camera-origin dielectric-region initialization.
7. Route Film wavelength samples into the camera API.
8. Ensure nondispersive cameras leave all four wavelengths active.
9. Keep current camera code as an adapter or legacy implementation until parity.
10. Defer realistic-lens dispersion to PR 21.

### Gate

- [x] geometric rays match legacy/reference cameras for nondispersive settings.
- [x] camera weight is spectrally neutral where expected.
- [x] wavelength packet remains unchanged for pinhole/thin-lens cameras.
- [x] camera ray/Film sample ordering is deterministic.

### Completion record

Completed in PR 9 by `src/render/spectral_camera.h`,
`src/render/spectral_camera.cpp`, `tools/spectral-tests/pr9-camera-tests.cpp`,
and `tools/spectral-tests/run-pr9-camera-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr9-camera-tests.R`
- `tools/codex/install-local.sh`

## PR 10: Add Interaction, Shape adapters, Primitive, and Scene

### Goal

Separate geometry from material and light binding sufficiently for pbrt-style integration.

### Required changes

1. Add Interaction and SurfaceInteraction types.
2. Add a Shape handle or adapter around existing hitables.
3. Add GeometricPrimitive, TransformedPrimitive, and Aggregate abstractions.
4. Bind Material, AreaLight, MediumInterface, and optional DielectricRegion at Primitive level.
5. Convert current hit records into SurfaceInteraction without losing UVs, derivatives, face IDs, transforms, and error bounds.
6. Preserve existing BVH traversal through an Aggregate adapter initially.
7. Add shape-area and directional-sampling adapters for emissive geometry.
8. Add robust ray spawning/offset functions and test them separately.
9. Ensure final geometric normal orientation is correct under instances and negative transforms.
10. Make spectral Scene own aggregate, lights, infinite lights, and immutable registries.
11. Add `ShapeCapabilities` and route validation through those capabilities.
12. Add a CSGShape adapter that exposes intersection, bounds, finite-difference normals, optional signed-distance containment, and topology flags.
13. Reject CSG UV texture, normal-map, displacement, area-light, and region-boundary uses that the shape capabilities do not support.

### Gate

- [x] intersection distances, UVs, normals, and instance transforms match legacy geometry tests.
- [x] Primitive bindings are correct for shared shapes/materials.
- [x] a shape can be used with a nonemissive Material and a separate AreaLight.
- [x] ray-origin offset tests avoid immediate self-intersections.

### Completion record

Completed in PR 10 by `src/render/spectral_scene.h`,
`src/render/spectral_scene.cpp`, `src/render/spectral_legacy_shape.cpp`,
`tools/spectral-tests/pr10-scene-tests.cpp`, and
`tools/spectral-tests/run-pr10-scene-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr10-scene-tests.R`
- `tools/codex/install-local.sh`

## PR 11: Add explicit Light, LightSampler, and InfiniteLight

### Goal

Port pbrt's light architecture before adding MIS.

### Required changes

1. Add Light flags/types and tagged dispatch.
2. Implement PointLight, SpotLight, DistantLight if currently supported, and DiffuseAreaLight.
3. Implement ImageInfiniteLight and constant infinite light.
4. Implement `SampleLi`, `PDF_Li`, `Le`, `L`, `Phi`, and `Preprocess` as appropriate.
5. Build infinite-map importance distributions in the pinned pbrt convention.
6. Reconstruct RGB light colors and environment maps using illuminant spectra.
7. Implement UniformLightSampler.
8. Add wavelength-independent power estimates and PowerLightSampler infrastructure.
9. Exclude the legacy emissive environment sphere from spectral Scene geometry.
10. Add visibility endpoints/interactions needed by direct lighting.
11. Convert legacy emissive materials and spot-light materials through the schema adapter.
12. Compile `with_light(area_light(...))` into a Primitive-bound DiffuseAreaLight.
13. Compile `add_light()` free lights and `set_environment()` infinite lights into scene-level Light registries.
14. Reject area lights attached to unsampleable shapes, including CSG without `csg_mesh_sampler(mode = "render_mesh")`, unless advanced unsampled-emitter mode is explicitly requested.
15. Verify that Material evaluation has no emission path in spectral mode.

### Gate

- [x] `SampleLi` histograms agree with each light's PDF.
- [x] area-light solid-angle PDFs agree with shape sampling.
- [x] infinite-light PDF includes the correct spherical Jacobian and transform.
- [x] sampled and directly evaluated emission agree.
- [x] light-selection PMFs sum to one and do not depend on the current wavelength packet.

### Completion record

Completed in PR 11 by `src/render/spectral_light.h`,
`src/render/spectral_light.cpp`, `tools/spectral-tests/pr11-light-tests.cpp`,
and `tools/spectral-tests/run-pr11-light-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr11-light-tests.R`
- `tools/codex/install-local.sh`

## PR 12: Add the BxDF/BSDF foundation

### Goal

Port pbrt's scattering contracts and coordinate conventions independent of Material parsing.

### Required changes

1. Add BxDF flags, reflection/transmission sample flags, TransportMode, and BSDFSample.
2. Add local reflection-coordinate helpers and the BSDF frame wrapper.
3. Port Fresnel dielectric and conductor functions.
4. Port Trowbridge-Reitz distribution, roughness remapping, and sampling.
5. Port `DiffuseBxDF` and any minimal null/interface BxDF required for tests.
6. Port `f`, `Sample_f`, PDF, rho, and regularization dispatch.
7. Implement pbrt's shading-normal orientation and world/local transformations.
8. Add standalone statistical sampling harnesses.
9. Do not connect to current Material or integrators yet.

### Gate

- [x] diffuse `f` and PDF integrate to expected reflectance.
- [x] microfacet distribution sampling/PDF agreement passes chi-square or equivalent tests.
- [x] Fresnel functions match high-precision references.
- [x] reciprocity/transport-mode tests pass where applicable.
- [x] no legacy cosine convention is imported.

### Completion record

Completed in PR 12 by `src/render/spectral_bsdf.h`,
`src/render/spectral_bsdf.cpp`, `tools/spectral-tests/pr12-bsdf-tests.cpp`,
and `tools/spectral-tests/run-pr12-bsdf-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr12-bsdf-tests.R`
- `tools/codex/install-local.sh`

## PR 13: Add pbrt-style Material closures and diffuse material

### Goal

Introduce immutable Material parameter objects and scratch-allocated per-hit BxDFs.

### Required changes

1. Add Material tagged dispatch and `Material::GetBSDF()`.
2. Add `MaterialEvalContext` derived from texture evaluation context.
3. Allocate concrete BxDFs from ScratchBuffer and return BSDF wrappers.
4. Implement DiffuseMaterial using SpectrumTexture reflectance.
5. Implement interface/null Material behavior.
6. Add bump mapping with FloatTexture and preserve geometric normals.
7. Add alpha evaluation hooks before Material closure construction.
8. Add material capability queries used by texture evaluators.
9. Ensure the closure captures evaluated spectral values, not texture handles where avoidable.
10. Add closure-lifetime and scratch-reset tests.

### Gate

- [x] Material-created diffuse BSDF equals directly constructed DiffuseBxDF.
- [x] texture evaluation happens once per closure construction under instrumentation.
- [x] no heap allocation occurs per diffuse hit.
- [x] reset invalidates closures only after the sample has completed.
- [x] alpha/bump tests preserve region and geometric-normal invariants.

### Completion record

Completed in PR 13 by `src/materials/spectral_material.h`,
`src/materials/spectral_material.cpp`, `tools/spectral-tests/pr13-material-tests.cpp`,
and `tools/spectral-tests/run-pr13-material-tests.R`.

Gate evidence:

- `Rscript tools/spectral-tests/run-pr13-material-tests.R`
- `tools/codex/install-local.sh`

## PR 14: Implement the first spectral RandomWalk vertical slice

### Goal

Render a complete spectral image with diffuse surfaces and explicit emitters before MIS complexity.

### Required changes

1. Add `render_mode = "spectral"` and `integrator_type = "randomwalk"` behind an experimental flag.
2. Compile diffuse Materials, area lights, and infinite lights into spectral Scene objects. Delta lights are deferred to the PathIntegrator because RandomWalk cannot reach them by ordinary path continuation.
3. Wire Film wavelength sampling, CameraRay, ScratchBuffer, RandomWalkIntegrator, and Film output.
4. Support vacuum only and no nested dielectric regions in this first slice.
5. Add deterministic sampler plumbing and per-worker state.
6. Add NaN/Inf/negative-radiance diagnostics.
7. Add minimal output and preview conversion through Film.
8. Keep adaptive sampling, denoising, and animation disabled for this experimental mode.

### Gate

- a diffuse Cornell-box-style scene renders correctly.
- narrow-band red, green, and blue emitters produce expected sensor responses.
- RGB albedo and RGB illuminant roles produce visibly and numerically distinct results.
- equivalent pbrt random-walk scenes agree within Monte Carlo confidence bounds.
- legacy mode remains unchanged.

## PR 15: Port the pbrt surface PathIntegrator with MIS

### Goal

Implement the production spectral surface path tracer for vacuum diffuse scenes.

### Required changes

1. Port the pinned pbrt PathIntegrator loop and state ordering.
2. Add direct-light sampling through LightSampler.
3. Add BSDF sampling and power-heuristic MIS.
4. Add area/infinite emitter-hit MIS.
5. Add specular-bounce and previous-light-context state.
6. Add Russian roulette and `etaScale` scaffolding.
7. Add depth accounting, cancellation, and robust failure checks.
8. Add UniformLightSampler first; enable PowerLightSampler only after validated.
9. Preserve random dimension assignments in documentation.
10. Add optional regularization only after unregularized parity.

### Gate

- direct-light and BSDF-sampled estimators are individually unbiased in test scenes.
- combined MIS matches high-sample reference images.
- environment hits receive correct MIS weights.
- roulette on/off images agree statistically.
- maximum-depth semantics match the pinned pbrt reference.

## PR 16: Port spectral conductor materials

### Goal

Add the first wavelength-dependent surface model while retaining all four path wavelengths.

### Required changes

1. Port ConductorBxDF and ConductorMaterial exactly.
2. Add named eta/k spectral datasets with provenance.
3. Add sampled eta/k R inputs.
4. Add isotropic and anisotropic roughness.
5. Add pbrt roughness remapping and smooth-limit behavior.
6. Add legacy metal adaptation from RGB/fuzz to a clearly named compatibility model.
7. Do not terminate secondary wavelengths for conductors.
8. Add conductor AOV reflectance estimates through BSDF rho or a documented approximation.

### Gate

- normal-incidence Fresnel matches analytic values at sampled wavelengths.
- rough and smooth conductor sampling/PDF tests pass.
- measured copper, silver, and gold scenes agree with pbrt references.
- packet remains unterminated after conductor interactions.
- legacy RGB-metal conversion emits one clear compatibility warning in spectral mode.


## PR 17: Add constant-IOR nested dielectric regions and smooth DielectricBxDF

### Goal

Introduce the region resolver and smooth nondispersive glass with exact transactional state behavior.

### Required changes

1. Add DielectricRegion, RegionId assignment, and exterior-region handling.
2. Add value-semantic DielectricPathState and transition analysis.
3. Add camera explicit-region initialization and analytic-shape automatic containment.
4. Port smooth DielectricBxDF and DielectricMaterial for constant eta.
5. Add `ResolvedDielectricInterface` to MaterialEvalContext.
6. Implement priority-skipped and index-matched null traversal.
7. Commit state only after transmission; reflection/TIR leave it unchanged.
8. Set ray Medium from active region, initially vacuum-only handles.
9. Add visibility traversal through skipped/null boundaries using copied state.
10. Add region diagnostics and debug counters.
11. Adapt legacy dielectric priority and scalar refraction into schema v2.
12. Do not add spectral absorption yet; reject nonzero `sigma_a` in production spectral scenes.
13. Support multi-primitive region-boundary assemblies from Section 12.14, including repeated references to one RegionId.
14. Support analytic and closed CSG `region_boundary(side = "inside")` for constant-IOR region tests where CSG capabilities indicate coherent containment.
15. Commit multi-region changes atomically when one interface declares multiple region boundaries.

### Gate

- all state-machine tests in Section 12.18 that do not require dispersion/media pass.
- concentric and overlapping constant-IOR reference scenes are declaration-order invariant.
- pbrt-equivalent single-interface glass scenes match pbrt.
- reflection, transmission, and TIR commit behavior is proven by unit tests.
- camera-inside analytic shapes initialize correctly.

## PR 18: Add rough and dispersive nested dielectrics

### Goal

Complete pbrt dielectric scattering and hero-wavelength termination for effective region interfaces.

### Required changes

1. Port full rough DielectricBxDF `f`, `Sample_f`, and PDF behavior.
2. Add spectral eta descriptors, sampled IOR, Cauchy, and Sellmeier forms.
3. Implement proof-based effective-ratio constancy detection.
4. Call `TerminateSecondary()` during closure evaluation for active nonconstant eta ratios.
5. Preserve all wavelengths at skipped and index-matched null boundaries.
6. Use eta inside/outside relative to geometric normal and pass their ratio to the pbrt closure.
7. Update `etaScale` from `BSDFSample::eta` exactly.
8. Add ThinDielectricMaterial separately; it never modifies region membership.
9. Add named glass spectra/data needed by tests.
10. Add path-level diagnostics showing wavelength termination location.
11. Keep pbrt's early termination even when the sampled branch reflects.

### Gate

- pbrt single-interface smooth and rough dielectric comparisons pass.
- dispersive prism/reference scenes agree with pbrt statistically and chromatically.
- skipped dispersive region tests prove no wavelength termination.
- active dispersive reflection still shows secondary termination, matching pbrt.
- rough reflection/transmission state commits pass.
- thin dielectric never changes region membership.

## PR 19: Port remaining surface materials and closure composition

### Goal

Replace remaining spectral surface behavior with dedicated pbrt-style closures.

### Required changes

1. Port CoatedDiffuseBxDF/Material.
2. Port CoatedConductorBxDF/Material.
3. Port HairBxDF/Material and spectral absorption parameterizations.
4. Port any pbrt diffuse-transmission or measured-material support required by rayrender scenes.
5. Map legacy `glossy` to the physically closest coated model through schema adaptation.
6. Decide Oren-Nayar handling:
   - port the pinned pbrt diffuse-rough model if available at the pinned commit; or
   - retain a separately validated Oren-Nayar BxDF with an ADR documenting the extension.
7. Add stochastic MixMaterial only if existing rayrender behavior requires material mixing; follow pbrt's material-selection semantics rather than constructing an unvalidated lobe graph.
8. Add regularization implementations for supported BxDFs.
9. Remove any spectral fallback to legacy `material::scatter()`.

### Gate

- each BxDF passes sampling/PDF and furnace tests.
- coated material images agree with pbrt references.
- hair sampling and absorption tests pass.
- every spectral schema-v2 surface Material resolves to a new closure or a clear unsupported error.
- no Material emits radiance directly.

## PR 20: Add Medium, PhaseFunction, absorption, and VolPathIntegrator

### Goal

Port pbrt's spectral volume architecture and make dielectric-region bulk properties physically correct.

### Required changes

1. Add Medium and PhaseFunction tagged handles.
2. Port Henyey-Greenstein and any pbrt phase functions required by current features.
3. Implement HomogeneousMedium with spectral sigma_a, sigma_s, emission, and scale.
4. Implement the heterogeneous/majorant medium needed to replace current constant-density hitables.
5. Port pbrt's majorant transmittance sampling and null-collision logic.
6. Port VolPathIntegrator rather than patching surface PathIntegrator.
7. Bind active dielectric regions to Medium handles.
8. Ensure segment medium is active-before and updated only after committed crossings.
9. Implement transmittance-aware visibility through null boundaries and media.
10. Map legacy dielectric attenuation to an absorption Medium using unbounded RGB reconstruction.
11. Map legacy constant-medium constructs through documented sigma values.
12. Add medium emission if supported by the pinned pbrt stage.
13. Add strict rejection of incompatible surface PathIntegrator plus non-null Medium scenes.

### Gate

- homogeneous Beer-Lambert transmittance matches analytic spectra.
- pure absorption nested glass scenes match analytic/reference output.
- homogeneous scattering scenes agree with pbrt VolPath references.
- null-collision estimator is unbiased under multiple majorants.
- local visibility-state copies traverse medium transitions correctly.
- current constant-medium reference scenes have documented spectral replacements.

## PR 21: Port realistic-camera spectral behavior and measured sensors

### Goal

Make camera optics and sensing spectral where the existing camera model supports it.

### Required changes

1. Port pbrt-compatible measured PixelSensor loading and calibration.
2. Add lens-element glass descriptors with constant, Cauchy, Sellmeier, or named-glass IOR.
3. Update realistic camera ray tracing to evaluate lens eta at `lambda[0]`.
4. Terminate secondary wavelengths before a wavelength-dependent lens direction is selected, following pbrt.
5. Preserve packet wavelengths for nondispersive lens systems.
6. Add wavelength-dependent camera weights and exit-pupil behavior as required.
7. Add sensor white-balance controls after sensor integration.
8. Add chromatic-aberration diagnostic scenes and ray tests.
9. Document any differences from pbrt due rayrender lens-file formats.
10. Ensure camera-origin region initialization occurs after final generated ray origin/time are known.

### Gate

- nondispersive realistic camera remains geometrically consistent with prior/reference behavior.
- dispersive lens rays match independently computed Snell-law tests.
- measured sensor integrations match deterministic spectral integrals.
- chromatic-aberration reference images are stable and physically ordered.
- camera wavelength termination is recorded correctly.

## PR 22: Complete importers, lights, and scene compilation

### Goal

Ensure all supported scene-ingestion paths produce correct spectral descriptors and explicit renderer objects.

### Required changes

1. Audit OBJ/MTL and other importer color assumptions.
2. Assign source color spaces/encodings explicitly or expose importer options.
3. Convert diffuse/base colors as albedo and emission as illuminant.
4. Map optical constants, transmission, opacity, and IOR where formats provide them.
5. Add remaining rayrender light types as explicit Light implementations.
6. Complete area-light attachment for all emissive shape types.
7. Ensure instances have unique RegionIds and correct transformed normals.
8. Add scene compilation caches for spectra, textures, Materials, Lights, Media, and shapes.
9. Add scene-level diagnostics summarizing legacy approximations and unsupported fields.
10. Reject positional schema data after the centralized adapter boundary.
11. Add deterministic compiler output hashes for tests.
12. Complete CSG scene compilation for ordinary spectral surfaces, CSG region boundaries, CSG null medium boundaries, and CSG-to-mesh area-light rendering.
13. Reject per-child CSG materials/regions, unsampled strict-mode CSG area lights, and unsupported CSG mappings with precise diagnostics.
14. Add importer hooks that attach `shape_capabilities` and generated mappings when imported geometry lacks UVs or normal support.

### Gate

- importer fixtures produce expected schema-v2 descriptors.
- emissive imported geometry participates correctly in direct and hit-light MIS.
- no imported scalar/data map is gamma decoded.
- shared Material plus multiple instances produces unique regions.
- compiler output is stable across runs and thread counts.

## PR 23: Add remaining RGB spaces, spectral assets, packaging, and performance work

### Goal

Finish color-space coverage and make the implementation practical as an R package without altering estimator behavior.

### Required changes

1. Add DCI-P3, Rec.2020, and ACES2065-1 color spaces and tables if package policy permits.
2. Add remaining named pbrt spectra needed by supported Materials/Lights/sensors.
3. Complete binary asset validation, lazy loading, and installed-package lookup.
4. Profile spectrum dispatch, texture reconstruction, closure allocation, LightSampler, Film merging, and region-state copying.
5. Add safe optimizations:
   - inline fixed-size spectrum arithmetic;
   - cache immutable spectrum/table handles;
   - use compact tagged dispatch;
   - reserve/inlined storage for region membership;
   - avoid repeated RGB reconstruction for constant descriptors;
   - batch Film tile merges.
6. Add per-subsystem counters and optional profiler labels.
7. Measure source package, installed package, and memory-map/load impact.
8. Add single- and multi-thread benchmarks against legacy and pbrt-equivalent scenes.
9. Do not introduce biased spectral clamping, packet resampling, or wavelength-dependent light-selection shortcuts.
10. Run sanitizers and undefined-behavior checks on the spectral test suite.

### Gate

- all color-space reconstruction tests pass.
- package checks pass without network access.
- no new estimator bias is detected by statistical tests.
- performance targets in Section 16 are met or a measured exception is documented.
- peak memory and asset sizes are recorded.

## PR 24: Restore feature parity, document, and define rollout

### Goal

Integrate spectral rendering with rayrender's user-facing workflow and decide when it becomes the default.

### Required changes

1. Route still and animation through the same spectral RenderSession.
2. Add adaptive sampling based on Film/sensor channels.
3. Add spectral-compatible alpha, albedo, normal, depth, emission, and variance AOVs.
4. Define denoiser inputs in output-linear RGB or another documented tristimulus space.
5. Restore preview, progress, cancellation, debug channels, and file-output parity.
6. Add spectral scene/material/light/camera documentation and migration examples.
7. Add a pbrt-conformance report with every intentional deviation.
8. Add a legacy compatibility guide, warning schedule, and numerical appearance-change examples.
9. Decide default policy only from acceptance data:
   - keep `rgb_legacy` default for one release; or
   - make `spectral` default with explicit release notes and a legacy escape hatch.
10. Do not remove legacy renderer code in the same PR that changes the default.
11. Open a separate future deprecation plan after at least one stable release.

### Gate

- all Section 15 test suites pass.
- still and animation frame zero agree under identical camera/scene/options.
- adaptive versus fixed sampling is statistically consistent.
- denoising/AOV pipelines contain no packet-component misuse.
- documentation builds and examples run.
- release candidate passes package checks on all supported platforms.


# 14. Test and validation program

Use the repository's existing `testthat` integration for R tests and embedded C++ unit tests where practical. Add a dedicated spectral reference-render harness under `tools/spectral-tests/`. Expensive statistical and image tests MAY run outside routine package checks, but every PR gate must have a fast deterministic subset.

## 14.1 Test tiers

1. **Tier A: deterministic unit tests**
   - run on every development test invocation;
   - no image comparison and no statistical flakiness;
   - target spectra, color, BSDF formulas, state transitions, parsing, and assets.
2. **Tier B: statistical estimator tests**
   - fixed sample counts and multiple deterministic seeds;
   - compare confidence intervals, means, PDFs, and known integrals;
   - run in CI where practical and always before a PR gate is accepted.
3. **Tier C: small reference renders**
   - low resolution, controlled sample count, stored high-sample references;
   - compare linear sensor/output channels with robust metrics;
   - include both pbrt and rayrender references where scenes are equivalent.
4. **Tier D: release validation renders**
   - larger/high-sample corpus;
   - performance, convergence, multithreading, animation, media, and realistic camera;
   - run before releases and major renderer changes.

## 14.2 Spectrum and wavelength unit tests

Required tests include:

- visible wavelength PDF integrates to one;
- inverse-CDF sampling matches reference quantiles;
- correlated packet samples match pbrt for fixed `u` values;
- uniform sampler bounds and PDFs;
- `TerminateSecondary()` idempotence;
- active PDF sum/weight behavior after termination;
- spectrum arithmetic and `SafeDiv` edge cases;
- blackbody shape at multiple temperatures;
- piecewise interpolation at knots and between knots;
- extrapolation policies;
- nonfinite/negative validation by semantic type;
- dense versus analytic spectrum sampling agreement;
- named asset checksums and lookup errors.

## 14.3 Color and RGB reconstruction tests

Required tests include:

- sRGB transfer-function round trip;
- RGB-to-XYZ and XYZ-to-RGB matrix round trip;
- whitepoint/primary matrix construction;
- table interpolation at cells, boundaries, and neutral axis;
- role-specific reconstruction for black, white, primaries, secondaries, gray ramp, and random RGB samples;
- reconstructed albedo boundedness;
- unbounded scale handling;
- illuminant multiplication by the color-space illuminant;
- deterministic integration of reconstructed spectra back to XYZ/RGB;
- cross-platform table checksum and byte-order tests;
- corrupt/truncated table diagnostics.

Use the pinned pbrt implementation or generated golden values as the numerical reference. Store the generation commit and tolerance rationale with each fixture.

## 14.4 Texture tests

Required tests include:

- FloatTexture and SpectrumTexture dispatch;
- constant, checker, mix, scale, and procedural values;
- UV/3D mapping transformations;
- derivative propagation;
- mip/filter behavior at known footprints;
- sRGB decode before filtering;
- scalar map never decoded as color;
- image cache key semantic separation;
- albedo versus illuminant reconstruction of the same filtered RGB texel;
- no heap allocation in repeated evaluations;
- thread-safe immutable texture-cache access.

## 14.5 BSDF numerical tests

For every BxDF:

1. evaluate flags for smooth/rough/reflection/transmission cases;
2. compare `f` to analytic or pbrt reference values at fixed direction pairs;
3. verify `Sample_f` output validity and event flags;
4. verify sampled directions against the reported PDF;
5. verify `PDF()` matches the sampling distribution for non-delta cases;
6. verify zero `f`/PDF behavior for delta distributions outside sampling;
7. test both TransportModes where non-symmetry matters;
8. test roughness extremes and grazing angles;
9. test finiteness and nonnegativity;
10. run white-furnace or directional-hemispherical energy tests.

Statistical sampling tests SHOULD use a robust goodness-of-fit method with predeclared significance correction. Do not accept or reject on one noisy seed.

## 14.6 Dielectric-specific numerical tests

Required tests include:

- dielectric Fresnel against analytic normal-incidence values;
- Brewster-angle and critical-angle behavior;
- smooth reflection/transmission probabilities;
- radiance-mode eta-squared correction;
- Importance versus Radiance relation;
- rough generalized half-vector equations;
- rough reflection/transmission Jacobians;
- TIR corner cases;
- eta equal to one;
- eta inside/outside inversion from the sign of `wo`;
- `BSDFSample::eta` convention;
- `etaScale` update sequence through enter/exit events;
- hero termination at active spectral interfaces;
- no termination at constant-ratio or skipped interfaces.

## 14.7 Light tests

For every light:

- sampled radiance equals direct evaluation at the sample;
- `SampleLi` distribution matches `PDF_Li`;
- power estimate is finite, nonnegative, and wavelength-packet independent;
- LightSampler selection PMFs normalize;
- area-light side and two-sided behavior;
- spot falloff boundaries;
- point/distant delta flags and MIS handling;
- infinite-map transform, pole behavior, Jacobian, and PDF normalization;
- RGB illuminant reconstruction versus directly sampled SPD;
- emitter-hit and direct-light emission consistency.

## 14.8 Film and sensor tests

Required tests include:

- deterministic sensor integral versus Monte Carlo packet estimator;
- invariance across wavelength-sampling distributions;
- inactive secondary PDFs after termination;
- filter normalization and edge pixels;
- tile merge versus serial accumulation;
- sensor white balance and output-space conversion;
- no negative clipping before the documented output stage;
- high-dynamic-range accumulation;
- sample-weight handling;
- visible-surface/AOV storage without affecting beauty output;
- repeatability across thread counts within the chosen accumulation policy.

## 14.9 Integrator estimator tests

Construct minimal analytic scenes for:

- constant environment plus Lambertian plane;
- one point light plus diffuse receiver;
- one area light plus diffuse receiver;
- direct-light sampling alone;
- BSDF sampling alone;
- combined MIS;
- specular mirror paths;
- specular glass paths;
- environment emitter-hit MIS;
- Russian roulette;
- maximum-depth boundary cases;
- null/interface boundaries;
- PowerLightSampler selection;
- regularization on/off.

For each stochastic estimator, compare the mean of multiple runs to an analytic or high-sample reference and verify the expected confidence interval coverage. Seed values are fixtures, not the correctness criterion by themselves.

## 14.10 Medium tests

Required tests include:

- pure absorption analytic transmittance;
- homogeneous single-scattering reference;
- zero-density limit;
- equal sigma_t channels;
- strongly spectrally varying sigma_a/sigma_s;
- phase-function normalization and sampling/PDF agreement;
- majorant validity and null-collision ratios;
- heterogeneous medium with a constant-density reduction;
- medium emission;
- visibility transmittance through multiple null boundaries;
- active region Medium before/after transitions;
- camera starting inside a Medium;
- environment path through exterior Medium if supported.

## 14.11 Render-reference scene corpus

Store scene generators, not opaque serialized binaries, where possible. Required scenes:

### Spectral/color

- equal-energy white source;
- standard illuminant over neutral diffuse patches;
- monochromatic/narrow-band emitters at multiple wavelengths;
- blackbody temperature sequence;
- RGB albedo versus the same RGB used as illuminant;
- metameric reflectance pair under two illuminants;
- color-checker-like measured reflectances if licensing permits.

### Surface materials

- diffuse Cornell box;
- smooth and rough measured conductors;
- smooth and rough nondispersive glass;
- dispersive prism;
- thin dielectric sheet;
- coated diffuse and coated conductor;
- hair under broadband and narrow-band light.

### Lighting

- point, spot, area, distant, and environment lights;
- bright small source for MIS stress;
- high-dynamic-range environment map;
- multiple-light power-sampling scene.

### Nested regions

All scenes listed in Section 12.19.

### CSG

- ordinary diffuse CSG object with constant reflectance;
- CSG conductor object with generated normals;
- CSG closed implicit glass sphere using `region_boundary(side = "inside")`;
- CSG subtractive cavity as a dielectric region;
- CSG onion shell as a nested dielectric shell;
- CSG pure medium boundary with `interface_material()`;
- CSG area light rejected without a sampler in strict mode;
- CSG area light accepted with `csg_mesh_sampler(mode = "render_mesh")`;
- CSG UV texture rejected unless generated mapping is provided;
- CSG smooth blend diagnostic when used as an optical region boundary.

### Media

- absorbing slab;
- homogeneous fog;
- colored/spectral fog;
- absorbing nested glass and bubble;
- heterogeneous density field;
- medium boundary with equal IOR.

### Camera

- pinhole/thin lens parity;
- realistic nondispersive lens;
- dispersive lens/chromatic aberration;
- measured sensor comparison;
- camera inside nested regions.

## 14.12 pbrt comparison protocol

For each equivalent reference scene:

1. pin the pbrt executable commit and build options;
2. document every scene-unit, camera, light, material, sensor, filter, sampler, and integrator mapping;
3. use the same resolution, sample count, maximum depth, and wavelength policy where possible;
4. disable denoising, tone mapping, adaptive sampling, and optional clamping;
5. compare linear XYZ or matched linear output RGB, not encoded display pixels;
6. crop or mask only documented non-equivalent regions;
7. use multiple seeds for stochastic comparisons;
8. report mean bias, RMSE, relative RMSE over nonblack pixels, 95th/99th error percentile, and structural diagnostics;
9. inspect convergence by increasing samples per pixel;
10. retain both images, metadata, and comparison report as versioned artifacts where repository size permits.

A low-sample image mismatch is not sufficient to declare a formula wrong. Check reference variance, sample streams, filter, and unit mappings first.

## 14.13 Tolerance policy

Do not put one global pixel tolerance on all tests.

- algebraic unit tests use tight precision-appropriate tolerances;
- table and asset tests use tolerances derived from storage precision;
- stochastic tests use confidence intervals or distribution tests;
- image tests use scene-specific thresholds established from repeated reference renders;
- cross-platform accumulation differences use a documented tolerance, not exact hashes;
- legacy tests retain their existing expected behavior unless intentionally changed.

Every tolerance above ordinary floating-point roundoff must include a comment explaining its source.

## 14.14 Failure diagnostics

A spectral debug build SHOULD be able to report, for a selected pixel/sample:

- sampled wavelengths and PDFs;
- camera weight;
- every hit and geometric/shading normal;
- active dielectric regions before and after each event;
- closure type and evaluated parameters;
- BSDF/light samples and PDFs;
- MIS weights;
- throughput after each update;
- Medium events and transmittance;
- wavelength termination point;
- roulette probability;
- final sensor contribution.

The trace must be deterministic, bounded, and disabled by default. It is a debugging aid, not a production logging path.


# 15. Feature integration requirements

## 15.1 Adaptive sampling

Adaptive decisions operate on accumulated Film/sensor channels. They MUST NOT treat the four wavelength packet values as four color channels.

Recommended procedure:

1. maintain per-pixel running means and variances in linear sensor or output-linear RGB;
2. use a luminance-like scalar plus an optional chromatic error criterion;
3. account for filter reconstruction consistently;
4. define minimum samples before adaptation;
5. keep wavelength sampling independent for every camera sample;
6. validate adaptive output against fixed-sample estimates.

The existing adaptive sampler may be reused only after its RGB accumulation assumptions are removed or isolated behind Film.

## 15.2 AOV semantics

Required spectral-mode AOVs:

| AOV | Definition |
|---|---|
| beauty | sensor-integrated transported radiance |
| albedo | pbrt-style visible-surface spectral albedo integrated through the selected sensor/output transform |
| normal | geometric or shading normal, explicitly named, in a documented coordinate system |
| depth | camera-space or ray distance with documented units |
| position | render/world coordinates if exposed |
| emission | sensor-integrated emitted contribution where supported |
| direct/indirect | optional sensor-integrated contribution split, with estimator definition documented |
| variance/sample count | Film estimator metadata |
| alpha | documented coverage/transmittance semantics, not spectral packet opacity |
| region ID/material ID | integer debug channels |

Albedo for denoising SHOULD follow pbrt's visible-surface/GBuffer logic rather than evaluating one arbitrary RGB texture field. For complex closures, use pbrt's rho or documented closure-specific estimates.

## 15.3 Denoising

Denoisers consume three-channel linear data. Convert beauty and albedo through the same sensor/output-linear basis before denoising. Normals remain geometric vectors.

Do not:

- denoise packet components;
- divide or multiply denoiser inputs by wavelength PDFs again;
- encode sRGB before denoising;
- hide negative sensor values by clipping before the denoiser unless the denoiser contract explicitly requires it and the policy is documented.

## 15.4 Alpha and transparent cutouts

Follow pbrt's stochastic alpha approach where practical:

- evaluate alpha before Material closure or region transition;
- an alpha-rejected hit continues the same ray direction;
- no region membership changes;
- no medium change;
- no depth increment;
- preserve prior MIS state;
- visibility rays use the same alpha policy with their local sampler/state.

A transmissive dielectric is not an alpha cutout.

## 15.5 Animation

Still and animation MUST share:

- scene schema/compiler;
- spectrum and image caches;
- RGB-to-spectrum tables;
- Film/PixelSensor implementation;
- Material/Light/Medium registries;
- integrator code;
- output color conversion;
- region initialization logic.

Frame-varying state includes camera transforms, animated object transforms, time-dependent lights/materials if supported, and per-frame output. Immutable spectral assets remain cached across frames.

For motion blur, camera-region initialization and animated Region intersections use the sampled ray time.

## 15.6 Units and radiometric semantics

New APIs MUST document physical quantity type even if rayrender retains relative-unit convenience:

- PointLight input is spectral intensity.
- DiffuseAreaLight input is spectral radiance.
- InfiniteLight input is spectral radiance by direction.
- `sigma_a` and `sigma_s` are inverse scene-length units before Medium scale.
- eta and reflectance are dimensionless.
- sensor response may be relatively normalized.

Legacy intensity multipliers remain adapters. The new renderer must avoid one generic `intensity` field whose physical meaning changes by light type without documentation.

## 15.7 Sampler and random dimensions

Add `docs/spectral/sample-dimensions.md` with stable assignments for:

- pixel filter sample;
- wavelength sample;
- lens sample;
- time sample;
- light selection;
- light direction/position;
- BSDF component and direction;
- alpha tests;
- medium distance/collision samples;
- phase-function samples;
- Russian roulette;
- layered-BxDF internal samples.

The exact sampler may continue to be rayrender's current implementation, but its API must support independent 1D/2D dimensions consistently. Any change to dimension order updates deterministic references intentionally.

## 15.8 Error handling and diagnostics

Worker threads MUST not call R APIs or emit unbounded warnings. Collect structured diagnostics per worker and merge them after rendering.

Report at least:

- invalid spectrum/texture/material counts;
- NaN/Inf path terminations;
- negative radiance/source-spectrum rejections;
- dielectric state inconsistencies;
- invalid BSDF/light PDFs;
- zero-throughput terminations;
- maximum null-boundary events;
- asset load failures;
- legacy approximation warnings.

A strict mode SHOULD turn all correctness diagnostics into errors after the first reproducible failing sample trace is saved.

# 16. Performance and packaging requirements

Correctness precedes optimization. Begin performance work only after the relevant conformance tests pass.

## 16.1 Hot-path requirements

The final spectral path MUST have:

- no heap allocation per bounce;
- no RGB-to-spectrum table file access after initialization;
- no repeated parsing of R objects during rendering;
- no virtual call chains that allocate temporary ownership objects;
- fixed-inline SampledSpectrum storage;
- scratch-allocated closures;
- immutable shared scene objects;
- per-worker samplers, scratch buffers, diagnostics, and Film tiles;
- inline or compact value-semantic region state for ordinary nesting depth;
- no global mutex on each sample or bounce.

## 16.2 Benchmark corpus

Benchmark at fixed resolution, sample count, and maximum depth:

1. diffuse Cornell box;
2. many-texture diffuse scene;
3. rough conductors;
4. nested smooth glass;
5. rough dispersive glass;
6. HDR infinite-light scene;
7. many-light scene;
8. homogeneous medium;
9. heterogeneous medium;
10. realistic camera.

Record:

- wall time;
- CPU time;
- rays and intersections per second;
- samples per second;
- peak resident memory;
- scratch high-water mark;
- spectrum/table/image cache size;
- region-state copy count and spill count;
- Film merge time;
- scaling by thread count.

## 16.3 Initial performance targets

Use the PR 0 machine and compiler baselines. After PR 23:

- diffuse spectral Path at equal camera samples SHOULD be no more than 3 times the legacy diffuse path wall time on the baseline CPU;
- typical region states of four or fewer memberships SHOULD not allocate dynamically;
- closure allocation SHOULD be less than 5 percent of diffuse render time;
- Film/sensor conversion SHOULD be less than 10 percent of diffuse render time;
- multithread parallel efficiency SHOULD not regress by more than 15 percentage points relative to legacy at the same thread count;
- table and named-spectrum initialization SHOULD be amortized and negligible after the first RenderSession;
- no benchmark may regress by more than 10 percent between accepted PRs without an explanation.

These are engineering targets, not permission to bias the estimator. If a target is missed, report the profile and retain correctness.

## 16.4 Binary assets and package installation

Do not generate large RGB-to-spectrum tables during ordinary package installation. Do not download them at first use.

Required packaging behavior:

- deterministic assets are part of the source or package-data workflow;
- installed-package lookup uses supported R package paths;
- source-package manifests include every asset;
- checksums are tested;
- binary format is portable or has a controlled conversion step;
- failures identify the missing/corrupt asset and color space;
- optional color-space tables may be split only if installation remains offline and reproducible.

Review source and installed package sizes before adding all color spaces. If size policy requires compromise, sRGB is mandatory and additional spaces may become optional packaged data with explicit installation documentation.

## 16.5 Licensing and provenance

pbrt-v4 source is under Apache-2.0, while rayrender declares GPL-3. Before copying source, data, or generated tables:

1. have the project maintainer review license compatibility and notice requirements;
2. preserve required copyright and SPDX notices in copied/adapted files;
3. identify whether generated tables are derivative assets and record their generator/source;
4. update `inst/COPYRIGHTS`, license notices, and provenance records as required;
5. prefer independent implementation from published formulas when provenance is unclear;
6. never strip attribution to reduce package size.

This plan is not legal advice; the repository must complete its normal licensing review.

# 17. Compatibility and rollout policy

## 17.1 Renderer coexistence

During PRs 0-24:

```text
rgb_legacy
  current Ray/material/scatter_record/pdf/integrators/environment sphere

spectral
  new Ray/Interaction/Primitive/Material closure/Light/Film/integrators
```

The two paths may share geometry math, shape intersection, image decoding infrastructure, samplers, progress reporting, and output libraries when semantics are compatible. They must not share ambiguous radiometric types or stateful material interfaces.

## 17.2 Legacy scene behavior

For `render_mode = "rgb_legacy"`:

- preserve current outputs and defaults;
- preserve current material parameter interpretation;
- preserve current environment behavior;
- preserve current dielectric pointer-stack implementation until separate retirement;
- do not route through schema-v2 spectral reconstruction unless explicitly tested as behavior preserving.

For `render_mode = "spectral"` with legacy constructors:

- adapt once to schema v2;
- produce an inspectable conversion report;
- warn for physically ambiguous mappings;
- never silently treat an emissive RGB as albedo or vice versa;
- preserve nested priority numbers and assign per-instance RegionIds.

## 17.3 Appearance changes

Spectral rendering is not expected to reproduce legacy RGB images exactly. Differences can result from:

- RGB-to-spectrum reconstruction roles;
- spectral Fresnel and optical constants;
- sensor/color-space conversion;
- corrected BSDF/cosine conventions;
- explicit MIS and infinite-light treatment;
- physically modeled absorption and dispersion;
- changed material mappings.

Document representative comparisons. Do not add arbitrary correction factors merely to force spectral mode toward legacy images.

## 17.4 Serialization/versioning

Any saved scene/cache format must include:

- schema version;
- renderer mode;
- color-space/encoding identifiers;
- RGB table version/checksum;
- named-spectrum asset version;
- pbrt conformance commit;
- material/light/medium descriptor versions.

Reject incompatible cached compiled scenes. Do not deserialize raw pointers, tagged handles, or scratch closures.

## 17.5 Default switch

The spectral renderer may become default only after:

- the Definition of Done is met;
- one release exposes spectral mode as opt-in or preview;
- package checks pass on supported systems;
- core legacy constructors have documented spectral mappings;
- performance and memory are acceptable;
- maintainers approve the appearance-change policy.

Retain `render_mode = "rgb_legacy"` for at least one stable release after the switch.

# 18. Risk register

## 18.1 Mixed BSDF conventions

**Risk:** legacy materials include inconsistent cosine/PDF behavior and contaminate the new integrator.

**Mitigation:** no spectral adapter may call legacy `scatter()` or `f()`. Port BxDFs from pbrt with standalone tests before integration.

## 18.2 Wavelength PDF misuse

**Risk:** division occurs in both throughput and Film, or disappears after hero termination.

**Mitigation:** centralize sensor integration, add estimator invariance tests, and instrument one debug assertion path showing exactly one division.

## 18.3 Incorrect RGB role

**Risk:** albedo, optical coefficients, and emission all use one RGB reconstruction.

**Mitigation:** require `SpectrumType` at every RGB conversion, make descriptors role-aware, and disallow default roles in ambiguous low-level constructors.

## 18.4 Asset provenance and package size

**Risk:** RGB tables or measured data create licensing, source-package, or platform problems.

**Mitigation:** pin provenance, use deterministic binary assets, validate checksums, begin with sRGB, and review notices before merge.

## 18.5 Material closure lifetime

**Risk:** a BSDF references scratch memory after reset or closures acquire heap-owned resources.

**Mitigation:** reset only after complete pixel-sample evaluation, use non-owning immutable references, test lifetime under sanitizers, and forbid closure escape in APIs.

## 18.6 Nested region state corruption

**Risk:** reflection mutates membership, copied rays alias state, or entry/exit uses bump normals.

**Mitigation:** pure Analyze plus checked Commit, value semantics, geometric-normal invariants, generation tokens, and exhaustive state-machine tests.

## 18.7 Shared-material region identity

**Risk:** two objects using one dielectric Material collapse into one region.

**Mitigation:** bind RegionId at Primitive/instance level and test shared Material with distinct instances.

## 18.8 Camera inside regions

**Risk:** initial state is empty, giving wrong eta and Medium.

**Mitigation:** explicit `initial_regions`, robust containment where supported, time-aware initialization, and hard failure on ambiguity.

## 18.9 Dispersive priority-skipped boundaries

**Risk:** secondary wavelengths are terminated merely because the crossed object's eta is spectral, even though the active optical region does not change.

**Mitigation:** resolve effective inside/outside regions before closure construction and terminate only for an active scattering interface with nonconstant effective ratio.

## 18.10 Rough dielectric commit errors

**Risk:** membership is inferred from direction or precommitted before rough reflection/transmission is selected.

**Mitigation:** use BSDFSample flags and commit only after a valid transmission sample.

## 18.11 Medium/surface conflation

**Risk:** absorption remains attached to the surface and is applied to the wrong segment or omitted on visibility paths.

**Mitigation:** map attenuation to Medium, use active-before segment state, and require VolPath for non-null media.

## 18.12 Environment MIS mismatch

**Risk:** sampled and hit-environment PDFs use different parameterization or transforms.

**Mitigation:** one explicit InfiniteLight implements both evaluation and PDF; test normalization and emitter-hit MIS independently.

## 18.13 Shape/light coupling

**Risk:** current hitables cannot expose pbrt-compatible area sampling or Primitive bindings cleanly.

**Mitigation:** adapters first, explicit capability diagnostics, then targeted shape refactors. Do not rewrite every intersection routine before the first vertical slice.

## 18.14 Image color management

**Risk:** encoded values are filtered directly, data maps are decoded, or cache entries alias roles.

**Mitigation:** explicit image descriptors, decode-before-filter rule, semantic cache keys, and texture tests.

## 18.15 Statistical test flakiness

**Risk:** fixed-seed pixel hashes either hide bias or fail on harmless platform differences.

**Mitigation:** separate deterministic algebra tests from multi-seed estimator tests and use predeclared confidence/tolerance policies.

## 18.16 Scope expansion

**Risk:** fluorescence, polarization, GPU work, or a generic material graph delays a correct surface spectral renderer.

**Mitigation:** enforce the exclusions in Section 3 and one-PR-at-a-time gates.

## 18.17 Still/animation divergence

**Risk:** wavelength, sensor, cache, or region initialization differs by entry point.

**Mitigation:** shared RenderSession in PR 1 and parity tests in PR 24.

## 18.18 Performance regression hidden by architecture work

**Risk:** dispatch, spectrum reconstruction, or state copying becomes prohibitively expensive.

**Mitigation:** preserve PR 0 benchmarks, instrument high-water/counters early, and optimize only after correctness gates.

## 18.19 Public API drift

**Risk:** the spectral renderer accumulates many flat `render_scene()` arguments, ad hoc list fields, and incompatible material aliases.

**Mitigation:** use schema-v2 descriptors, functional decorators, centralized adapters, and validation. Keep structured descriptors as the canonical API and legacy scalar arguments as adapters.

## 18.20 CSG contract mismatch

**Risk:** CSG can intersect rays but cannot provide pbrt-compatible area sampling, UVs, exact area, or reliable containment for every expression.

**Mitigation:** attach ShapeCapabilities to every compiled shape, support ordinary CSG and closed implicit CSG region boundaries early, reject unsupported CSG area lights/textures/regions in strict mode, and require tessellation with `mode = "render_mesh"` for CSG area lights.

## 18.21 CSG numerical boundary misses

**Risk:** sphere tracing or SDF thresholds miss a thin shell or tangent interface, corrupting dielectric-region state.

**Mitigation:** add scale-aware epsilon, root refinement, shell-thickness diagnostics, signed-distance containment checks, CSG-specific region tests, and bounded null-boundary traversal diagnostics.


# 19. Definition of Done

The spectral renderer is complete only when all of the following are true.

## 19.1 Architecture

- [ ] New spectral transport uses SampledSpectrum and SampledWavelengths throughout.
- [ ] Geometric vectors/points are not used as radiometric values in new code.
- [x] Material evaluation produces scratch-allocated per-hit BxDF closures.
- [ ] Spectral integrators never call legacy `material::scatter()`, `material::f()`, or `material::emitted()`.
- [ ] Primitive separates Shape, Material, AreaLight, MediumInterface, and DielectricRegion.
- [ ] explicit Light and LightSampler implementations drive direct lighting and MIS.
- [ ] Infinite lights are evaluated on ray misses, not through spectral environment geometry.
- [ ] Film and PixelSensor own spectral-to-tristimulus conversion.
- [ ] still and animation use one RenderSession.
- [x] hot-path per-hit closure allocation is scratch based and leak free.

## 19.2 Spectral representation

- [ ] pbrt visible wavelength sampling and four-sample packet match the pinned reference.
- [ ] secondary termination matches pbrt and is idempotent.
- [ ] wavelength PDF division occurs exactly once.
- [ ] CIE/illuminant/named spectrum assets are versioned and verified.
- [ ] sRGB albedo, unbounded, and illuminant reconstruction match pbrt.
- [ ] supported input encodings and color spaces are explicit.
- [ ] output conversion occurs after transport.

## 19.3 Materials and BSDFs

- [ ] Diffuse, conductor, dielectric, thin dielectric, coated, and hair models required by rayrender are ported or explicitly unsupported.
- [x] BxDF conventions, flags, PDFs, eta, and TransportMode match pbrt.
- [ ] conductor eta/k and dielectric eta are spectral.
- [x] roughness remapping and microfacet sampling match pbrt.
- [x] no new spectral Material combines unrelated lobes through unvalidated ad hoc probabilities.
- [ ] emissive behavior is represented by Light.
- [x] bump/normal mapping cannot alter geometric medium/region transitions.

## 19.4 Nested dielectrics

- [ ] each closed object instance has a stable RegionId distinct from Material identity.
- [ ] lower numeric priority wins in overlapping membership.
- [ ] Analyze is pure and Commit is transactional.
- [ ] reflection/TIR never commit; transmission commits once.
- [ ] skipped/null boundaries do not consume depth or disturb MIS state.
- [ ] effective outside/inside eta spectra are resolved relative to geometric normal.
- [ ] active dispersive ratios terminate wavelengths before closure sampling.
- [ ] skipped dispersive boundaries do not terminate wavelengths.
- [ ] ray Medium and active region Medium remain consistent.
- [ ] camera rays can initialize inside nested regions.
- [ ] visibility rays use independent copied region state.
- [ ] all Section 12 unit and render tests pass.

## 19.5 Lighting and integration

- [ ] PathIntegrator direct-light and emitter-hit MIS match pbrt behavior.
- [ ] all Light sample/PDF pairs are consistent.
- [ ] LightSampler PMFs are scalar and wavelength independent.
- [ ] Russian roulette and etaScale match pbrt.
- [ ] infinite lights participate correctly in both MIS paths.
- [ ] alpha/null traversal preserves prior scattering state.
- [ ] RandomWalk and Path vertical slices agree with reference scenes.

## 19.6 Media and cameras

- [ ] Medium/PhaseFunction interfaces and VolPath follow the pinned pbrt implementation.
- [ ] spectral absorption is applied along segments, not at ending surfaces.
- [ ] visibility includes medium transmittance where appropriate.
- [ ] realistic camera can evaluate wavelength-dependent lens IOR.
- [ ] measured sensor curves integrate correctly.
- [ ] camera/lens wavelength termination follows pbrt.

## 19.7 Product integration

- [ ] schema-v2 descriptors are named, versioned, validated, and documented.
- [ ] `ray_scene_v2` remains tibble-compatible while preserving scene-level registries.
- [ ] `ray_material_v2` stores pbrt-style named parameter descriptors, not legacy flat C++ payloads.
- [ ] `with_light()`, `with_region_boundary()`, `add_light()`, `set_environment()`, and `add_region()` work as functional scene/object operations.
- [ ] schema-v1 adaptation is centralized and reports approximations.
- [ ] existing scenes still render in `rgb_legacy` with preserved baselines.
- [ ] adaptive sampling, AOVs, denoising inputs, alpha, preview, and animation work through Film.
- [ ] package installation is offline and reproducible.
- [ ] spectral assets include provenance and required notices.
- [ ] package checks and sanitizers pass.
- [ ] performance and memory reports meet or explain Section 16 targets.
- [ ] pbrt comparison reports are checked in or archived reproducibly.
- [ ] user documentation explains physical parameters, scene descriptors, decorators, CSG limitations, and appearance changes.

## 19.8 CSG extension

- [ ] CSG is implemented as a Shape with capability flags, not as a special material/integrator path.
- [ ] ordinary CSG surfaces work with supported spectral Materials.
- [ ] closed/coherent CSG regions support `region_boundary(side = "inside")`.
- [ ] pure CSG medium boundaries work through `interface_material()` and VolPath.
- [ ] CSG area lights are rejected unless the emitting surface is sampleable, normally through `csg_mesh_sampler(mode = "render_mesh")`.
- [ ] unsupported CSG UV textures, normal maps, displacement, per-child materials, and per-child regions are rejected clearly.
- [ ] CSG-specific numerical and render tests pass.

# 20. Codex execution protocol

## 20.1 Before starting any PR

Codex must:

1. read `AGENT.md` and this plan;
2. read the current `docs/spectral/pbrt-conformance.md`;
3. inspect the exact rayrender files to be modified;
4. inspect the corresponding pinned pbrt source and book sections;
5. state the PR number and goal in the working notes;
6. list expected files and tests;
7. confirm that prerequisite PR gates are present in the branch;
8. avoid changing files unrelated to the PR unless a build/test necessity is documented.

## 20.2 Implementation rules

1. Follow the pinned pbrt implementation over memory or secondary summaries.
2. Preserve pbrt evaluation order where it affects random samples, wavelength termination, MIS state, eta scaling, or medium handling.
3. Add a short source comment naming the pbrt class/function for formulas whose conventions are easy to misread. Do not paste lengthy book text.
4. Use existing rayrender geometry/intersection code through adapters unless the PR explicitly calls for a rewrite.
5. Do not create temporary spectral versions of legacy `scatter_record`, raw `pdf*`, or RGB radiometry.
6. Do not put wavelengths or mutable region state on geometric Ray.
7. Do not let Material closures mutate integrator state.
8. Do not convert to RGB inside transport or use packet lanes as colors.
9. Do not use current packet values to choose LightSampler PMFs.
10. Do not add heap allocation to a per-hit or per-bounce path without a benchmark and approval.
11. Use `=` for all R assignments.
12. Use ASCII in source code, tests, generated code, and code snippets.
13. Preserve copied-source notices and update provenance in the same PR.
14. Remove debugging output before the PR gate unless it is behind the documented spectral trace option.
15. Update the conformance table and this plan's completion checklist in every PR.

## 20.3 How to handle discrepancies

When rayrender and pbrt disagree:

- preserve rayrender behavior in `rgb_legacy`;
- use pinned pbrt behavior in `spectral`;
- preserve rayrender nested-priority semantics only through the extension in Section 12;
- document any additional deliberate extension in an ADR;
- add a focused test demonstrating both the pbrt behavior and the extension;
- do not silently blend formulas or conventions.

When the pbrt book and pinned source appear to disagree, use the pinned source as executable reference and record the discrepancy. Verify whether a later correction exists, but do not silently change the pin.

## 20.4 Required commands and checks

At minimum, for each PR:

```sh
tools/codex/install-local.sh
```

Then run:

- package unit tests;
- the PR-specific C++ and R tests;
- the fast spectral conformance suite;
- relevant legacy render baselines;
- relevant reference renders;
- compiler warnings for the modified code;
- formatter/linter checks used by the repository;
- sanitizers when the PR changes ownership, scratch allocation, tagged dispatch, or path state.

Use repository-supported commands rather than inventing a parallel build system. Add helper scripts under `tools/spectral-tests/` when a reproducible command is missing.

## 20.5 PR report template

Each Codex PR response must use this structure:

```text
PR N: <title>

Implemented
- <major change>
- <major change>

pbrt references
- <pinned source file and symbol>
- <book section>

rayrender files changed
- <file>
- <file>

Tests
- <command>: PASS/FAIL
- <command>: PASS/FAIL

Legacy regression
- <scene/hash/result>

Spectral/reference result
- <metric/result>

Performance
- before: <value>
- after: <value>
- delta: <value>

Conformance deviations
- none
or
- <documented deviation and ADR>

Remaining work
- only items assigned to later PRs
```

Do not claim a PR is complete when a mandatory gate fails. Fix it within the same task or report the exact failing evidence and leave the completion checkbox unset.

## 20.6 Commit discipline

Each PR SHOULD contain logically separated commits:

1. interfaces/data assets;
2. implementation;
3. tests/reference fixtures;
4. documentation/provenance.

Generated assets and their generator changes belong in the same PR. Do not commit build products, local installed packages, or uncompressed high-sample images unless explicitly part of the reference-fixture policy.

## 20.7 Suggested first Codex prompt

Use this to start PR 0:

```text
Implement PR 0 from docs/spectral/rayrender_spectral_rendering_codex_plan.md.

Work only on the baseline, pinning, provenance, conformance-record, and reproducible-test tasks listed in PR 0. Do not modify rendering behavior. Read AGENT.md first and use tools/codex/install-local.sh. Inspect the attached/current rayrender source and pin the exact pbrt-v4 commit that will be used as the normative source reference.

Preserve all legacy tests and capture the required still/animation baselines and benchmarks. Add docs/spectral/versions.md, docs/spectral/pbrt-conformance.md, docs/spectral/provenance.md, and ADR 0001. Copy the plan into the repository. Run the PR gate and return the report in Section 20.5 format. Stop after PR 0.
```

## 20.8 Suggested nested-dielectric PR prompt

Use this to start PR 17 after all prerequisites pass:

```text
Implement PR 17 from docs/spectral/rayrender_spectral_rendering_codex_plan.md.

Follow Section 12 exactly. Preserve lower-number-wins priority behavior, but do not reuse the legacy vector<dielectric*> state in spectral mode. Add per-instance RegionId, immutable DielectricRegion, value-semantic DielectricPathState, pure Analyze(), checked Commit(), camera-origin initialization for supported analytic shapes, and pbrt smooth DielectricBxDF closure construction from resolved geometric-side eta values.

The Material closure must not mutate region state. Reflection and TIR leave state unchanged; transmission commits exactly once. Priority-skipped and exact index-matched boundaries pass straight through without consuming depth, wavelength termination, etaScale changes, or prior MIS-state replacement. Visibility traversal must use a copied state. Use geometric normals for entry/exit and retain shading normals only for BSDF evaluation.

This PR supports constant eta and vacuum media only. Reject nonzero absorption and defer spectral dispersion/roughness to PR 18. Add every applicable Section 12.18 state-machine test and the PR 17 reference renders. Run the full PR gate and return the Section 20.5 report. Stop after PR 17.
```

## 20.9 Suggested dispersive-dielectric PR prompt

Use this to start PR 18:

```text
Implement PR 18 from docs/spectral/rayrender_spectral_rendering_codex_plan.md.

Port the pinned pbrt-v4 rough DielectricBxDF and dielectric Material behavior, including Trowbridge-Reitz sampling, TransportMode scaling, BSDFSample::eta, and early SampledWavelengths::TerminateSecondary() for a nonconstant effective eta ratio. The ratio is etaInside(lambda) / etaOutside(lambda), with inside/outside defined by the primitive geometric normal after nested-region priority resolution.

Do not terminate wavelengths for priority-skipped or exact index-matched null boundaries. Do not optimize by postponing termination until transmission; match pbrt even when reflection is selected. Commit region state only for a valid transmission BSDFSample. Add sampled, Cauchy, Sellmeier, and named-glass eta descriptors plus ThinDielectricMaterial, which never changes region membership.

Run all dielectric formula, sampling/PDF, region-state, dispersion, and pbrt comparison gates. Return the Section 20.5 report and stop after PR 18.
```

## 20.10 Suggested public API PR prompt

Use this during PR 6:

```text
Implement the Section 9 public API and schema-v2 descriptor work assigned to PR 6.

Add typed descriptors for spectra, textures, materials, lights, media, optical regions, region boundaries, cameras, films, sensors, samplers, integrators, spectral options, and validation. Preserve existing constructors and add centralized legacy-to-schema-v2 adapters. Implement `ray_scene_v2` as a tibble-compatible scene with row-level list columns for `light`, `region_boundaries`, `medium_interface`, `shape_capabilities`, and stable object IDs. Implement scene registries for free lights, environments, regions, named textures, and named materials.

Implement `with_light()` as a functional area-light object decorator and `with_region_boundary()` as a functional region-boundary object decorator. Implement `add_light()`, `set_environment()`, and `add_region()` for scene-level registries. Do not expose C++ BxDF closures in R. `ray_material_v2` must be a named pbrt-style parameter descriptor, not a legacy flat material payload.

Add validation and tests showing that descriptors serialize deterministically, old scenes still construct, decorators preserve object rows and scene attributes, invalid descriptors fail before rendering, and R code uses `=` assignment. Do not modify renderer output in this PR.
```

## 20.11 Suggested CSG PR prompt

Use this during PR 10 or PR 22, according to the PR sequence:

```text
Implement the CSG spectral-interface requirements assigned to this PR.

Treat CSG as a rayrender-specific Shape implementation with explicit ShapeCapabilities. Add CSGShape adapters for bounds, intersection, normals, optional SignedDistance/Contains, topology flags, and validation. Ordinary CSG spectral surfaces are allowed. Closed/coherent CSG shapes may define optical regions through `region_boundary(side = "inside")`. Pure CSG medium boundaries use `interface_material()` and region metadata.

Reject CSG area lights unless the CSG is tessellated or otherwise made sampleable; the strict-mode supported path is `csg_mesh_sampler(mode = "render_mesh")`, where the rendered and sampled surfaces are identical. Reject UV image textures, normal maps, and displacement on CSG unless a supported generated mapping is supplied. Reject per-child CSG materials and regions.

Add diagnostics for thin shells, tangent hits, noncoherent topology, unsupported mappings, unsampled emitters, and no-progress null-boundary loops. Add tests for ordinary CSG surfaces, CSG glass regions, subtractive cavities, onion shells, pure medium boundaries, CSG area-light rejection, CSG area-light tessellation, and unsupported texture mapping.
```

# 21. Primary pbrt-v4 reference map

Pin exact source lines/commit in `docs/spectral/versions.md`; paths below describe the current pbrt-v4 organization and may move after the plan date.

| rayrender subsystem | pbrt-v4 source areas | pbrt book areas |
|---|---|---|
| SampledSpectrum/SampledWavelengths | `src/pbrt/util/spectrum.h`, `src/pbrt/util/spectrum.cpp` | Radiometry, Spectra, and Color; Representing Spectral Distributions |
| RGB color spaces and reconstruction | `src/pbrt/util/color.*`, `src/pbrt/util/colorspace.*`, RGB-to-spectrum table code/assets | Color and RGB-to-Spectrum Conversion |
| Film and PixelSensor | `src/pbrt/film.*` | Cameras and Film; Film and Imaging |
| Cameras | `src/pbrt/cameras.*`, base camera interfaces | Cameras and Film |
| Texture interfaces/evaluators | `src/pbrt/textures.*`, base texture interfaces | Textures and Materials |
| Material closure dispatch | `src/pbrt/materials.*`, base material interfaces | Material Interface and Implementations |
| BxDF/BSDF | `src/pbrt/bxdfs.*`, `src/pbrt/bsdf.*` | Reflection Models; BSDF Representation |
| DielectricBxDF | dielectric sections of `src/pbrt/bxdfs.*` and `src/pbrt/materials.*` | Dielectric BSDF; Material Implementations |
| Conductor/coated/hair | corresponding BxDF and Material implementations | Reflection Models; Materials |
| Interaction/Primitive | `src/pbrt/interaction.*`, `src/pbrt/base/primitive.h`, primitive implementations | Primitives and Intersection Acceleration |
| Lights | `src/pbrt/base/light.h`, `src/pbrt/lights.*` | Light Sources |
| LightSampler | light sampler source and interfaces | Light Sampling |
| Path/RandomWalk/VolPath integrators | `src/pbrt/cpu/integrators.*` | Light Transport I and II; Volume Scattering |
| Medium/PhaseFunction | `src/pbrt/base/medium.h`, `src/pbrt/media.*` | Volume Scattering; Media |
| Memory/scratch allocation | pbrt memory utilities and Material::GetBSDF implementation | System Overview; Material Interface |
| GBuffer/visible surface | Film/GBuffer implementation | Film and Imaging |
| CSG extension | no direct pbrt counterpart; must satisfy pbrt Shape/Primitive/Light contracts | Shapes; Light Sources; Volume Scattering |

For every ported function, use the pinned source rather than this table alone. Update the conformance document with exact symbol names and test fixtures as implementation proceeds. For CSG, record that it is a rayrender extension and cite the pbrt contract it is required to satisfy rather than a pbrt CSG implementation.

# 22. Final implementation stance

The shortest safe route is not to turn each existing RGB material into a four-float version of itself. The safe route is to preserve the current renderer as legacy and build a pbrt-style spectral core with explicit contracts.

The intentional architectural extensions are rayrender's overlapping nested-dielectric priority model and rayrender's CSG geometry. Nested regions remain a path-local spatial-state resolver, not part of Material or Ray geometry. They supply effective eta and Medium values to otherwise standard pbrt closures and integrators. CSG remains a Shape extension that must satisfy explicit pbrt-style capabilities before it can be used as a light or optical-region boundary. These boundaries keep the extensions testable and prevent them from changing pbrt BSDF mathematics.

The user-facing API remains functional and pipeable. New spectral descriptors are schema-versioned R objects, object decorators attach area lights and region boundaries, and scene-level registries hold free lights, environments, regions, and named resources.

Implement the PRs in order. Do not begin optimization, broad legacy deletion, or a generic closure graph before the corresponding correctness gates are complete.
