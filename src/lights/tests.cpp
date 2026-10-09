#ifdef NOT_CRAN
#include "infinite.h"
#include "../hitables/infinite_area_light.h"
#include <testthat.h>

namespace {
class DirectionalTexture final : public texture {
public:
  point3f value(Float u, Float v, const point3f &) const override {
    Float s = std::max(Float(0), Float(std::cos(2 * M_PI * u)));
    return point3f(.02 + 4 * s * s * v, .03 + v, .1);
  }
};
vec3f direction(Float z, Float phi) {
  Float r = std::sqrt(std::max(Float(0), 1 - z * z));
  return vec3f(r * std::cos(phi), z, r * std::sin(phi));
}
}

context("Image infinite lights") {
  test_that("latlong reconstruction preserves irradiance from bright seam texels") {
    const int width = 512, height = 256;
    std::vector<Float> pixels(size_t(width) * height * 3, .1f);
    for (int y = 129; y < 131; ++y)
      for (int x : {0, width - 1})
        for (int c = 0; c < 3; ++c) pixels[3 * (x + width * y) + c] = 5000;
    auto texture = std::make_shared<latlong_image_texture>(pixels.data(), width, height, 3);
    ImageInfiniteLight light(texture, width, height, 0);
    double reference = 0, estimate = 0;
    // Independent source-texel quadrature, not the texture's own lookup. The
    // old endpoint interpolation loses half of the bright seam texels' area.
    for (int y = 0; y < height; ++y) {
      const double theta = M_PI * (y + .5) / height;
      for (int x = 0; x < width; ++x) {
        const double phi = -2 * M_PI * (x + .5) / width;
        reference += pixels[3 * (x + width * y)] * std::max(0.0, std::sin(theta) * std::cos(phi)) * std::sin(theta);
      }
    }
    reference *= 2 * M_PI * M_PI / double(width * height);
    // Resolve the bilinear footprint directly. Importance samples would add
    // unrelated variance from the very dim neighbors of these bright texels.
    const int columns = width * 4, rows = height * 4;
    for (int y = 0; y < rows; ++y) {
      const double theta = M_PI * (y + .5) / rows;
      for (int x = 0; x < columns; ++x) {
        const double phi = -2 * M_PI * (x + .5) / columns;
        vec3f wi(std::sin(theta) * std::sin(phi), std::cos(theta),
                 std::sin(theta) * std::cos(phi));
        estimate += light.Radiance(point3f(0), wi, 0)[0] *
                    std::max(Float(0), wi[2]) * std::sin(theta);
      }
    }
    estimate *= 2 * M_PI * M_PI / double(columns * rows);
    expect_true(std::abs(estimate / reference - 1) < .02);
    // Sampling must also resolve the reconstructed footprint, including the
    // bright spill into neighboring cells whose center radiance is near zero.
    double sampled = 0;
    const int samples = 1024;
    for (int y = 0; y < samples; ++y) for (int x = 0; x < samples; ++x) {
      vec3f wi = light.Sample(point3f(0), vec2f((x + .5f) / samples, (y + .5f) / samples), 0);
      sampled += light.Radiance(point3f(0), wi, 0)[0] * std::max(Float(0), wi[2]) /
                 light.Pdf(point3f(0), wi, 0);
    }
    sampled /= double(samples * samples);
    expect_true(std::abs(sampled / reference - 1) < .02);
  }
  test_that("latlong texels wrap in longitude and clamp at the poles") {
    Float pixels[] = {2,2,2, 6,6,6, 10,10,10, 14,14,14};
    latlong_image_texture texture(pixels, 2, 2, 3);
    expect_true(texture.value(0, .75, point3f(0))[0] == Approx(4));
    expect_true(texture.value(1, .75, point3f(0))[0] == Approx(4));
    expect_true(texture.value(.25, .75, point3f(0))[0] == Approx(2));
    expect_true(texture.value(.25, 1, point3f(0))[0] == Approx(2));
    expect_true(texture.value(.25, 0, point3f(0))[0] == Approx(10));
    expect_true(texture.value(-.25, .25, point3f(0))[0] == Approx(14));
  }
  test_that("rotation transforms both radiance and directional PDF") {
    auto tex = std::make_shared<DirectionalTexture>();
    ImageInfiniteLight a(tex, 64, 32, 0), b(tex, 64, 32, 63);
    vec3f w = unit_vector(vec3f(.7, .4, .8)), rotated = RotateY(63)(w);
    point3f la = a.Radiance(point3f(0), w, 0), lb = b.Radiance(point3f(0), rotated, 0);
    expect_true((la - lb).length() < 1e-5);
    expect_true(std::abs(a.Pdf(point3f(0), w, 0) - b.Pdf(point3f(0), rotated, 0)) < 1e-5);
  }
  test_that("mixtures add radiance and evaluate every proposal density") {
    auto a = std::make_shared<ImageInfiniteLight>(std::make_shared<DirectionalTexture>(), 64, 32, 0);
    auto b = std::make_shared<ImageInfiniteLight>(std::make_shared<constant_texture>(point3f(0, 0, 3)), 32, 16, 110);
    InfiniteLightMixture mixture({a, b});
    vec3f w = unit_vector(vec3f(.3, -.4, .5));
    point3f p(2, 3, 4);
    expect_true((mixture.Radiance(p, w, 0) - (a->Radiance(p, w, 0) + b->Radiance(p, w, 0))).length() < 1e-6);
    double t = a->SamplingWeight() / (a->SamplingWeight() + b->SamplingWeight());
    double pdf = t * a->Pdf(p, w, 0) + (1 - t) * b->Pdf(p, w, 0);
    expect_true(std::abs(mixture.Pdf(p, w, 0) - pdf) < 1e-6);
    // Compare the sampled integral with independent equal-solid-angle quadrature.
    double estimated = 0, reference = 0, integrated_pdf = 0;
    const int n = 256;
    for (int y = 0; y < n; ++y) {
      for (int x = 0; x < n; ++x) {
        vec2f u(Float(x + .5) / n, Float(y + .5) / n);
        vec3f sample = mixture.Sample(p, u, 0);
        double q = mixture.Pdf(p, sample, 0);
        estimated += mixture.Radiance(p, sample, 0)[0] / q;
        vec3f uniform = direction(1 - 2 * u[0], 2 * M_PI * u[1]);
        reference += mixture.Radiance(p, uniform, 0)[0] * 4 * M_PI;
        integrated_pdf += mixture.Pdf(p, uniform, 0) * 4 * M_PI;
      }
    }
    expect_true(std::abs(integrated_pdf / (n * n) - 1) < .01);
    expect_true(std::abs(estimated / reference - 1) < .01);
  }
  test_that("black environments have finite sampling and zero emitted radiance") {
    auto black = std::make_shared<ImageInfiniteLight>(std::make_shared<constant_texture>(point3f(0)), 8, 4, 0);
    auto mix = std::make_shared<InfiniteLightMixture>(std::vector<std::shared_ptr<InfiniteLight>>{black, black});
    Transform identity;
    InfiniteAreaLight endpoint(mix, 100, point3f(0), &identity, &identity);
    random_gen rng(3);
    vec3f w = endpoint.random(point3f(0), rng);
    expect_true(std::isfinite(endpoint.pdf_value(point3f(0), w, rng)));
    expect_true(endpoint.pdf_value(point3f(0), w, rng) == Approx(1 / (4 * M_PI)));
    hit_record hit;
    Ray ray(point3f(0), w);
    expect_true(endpoint.hit(ray, 0, 1000, hit, rng));
    bool invisible;
    expect_true(hit.mat_ptr->emitted(ray, hit, hit.u, hit.v, hit.p, invisible).length() == 0);
  }
}
context("Celestial disk lights") {
  test_that("uniform disks have exact power and irradiance across angular sizes") {
    for (double diameter : {.53, 30., 170.}) {
      DiskInfiniteLight disk(point3f(2), vec3f(0, 1, 0), diameter, 0);
      double radius = diameter * M_PI / 360;
      double solid_angle = 4 * M_PI * std::pow(std::sin(radius / 2), 2);
      expect_true(std::abs(disk.SamplingWeight() / (2 * solid_angle) - 1) < 1e-6);
      double irradiance = 0;
      const int n = 128;
      for (int i = 0; i < n; ++i) {
        vec3f wi = disk.Sample(point3f(0), vec2f((i + .5f) / n, .37), 0);
        double pdf = disk.Pdf(point3f(0), wi, 0);
        expect_true(std::abs(pdf * solid_angle - 1) < 1e-6);
        irradiance += disk.Radiance(point3f(0), wi, 0)[0] * wi[1] / pdf;
      }
      expect_true(std::abs(irradiance / n / (2 * M_PI * std::pow(std::sin(radius), 2)) - 1) < 1e-5);
    }
    DiskInfiniteLight below(point3f(1), vec3f(0, -1, 0), 10, 0);
    expect_true(below.Radiance(point3f(0), vec3f(0, -1, 0), 0)[0] == 1);
  }
  test_that("disk proposals integrate irradiance and remain finite at lunar angles") {
    auto white = std::make_shared<constant_texture>(point3f(1));
    const double diameter = .53, radius = diameter * M_PI / 360;
    DiskInfiniteLight disk(white, 128, 128, vec3f(0, 1, 0), diameter, 0);
    double irradiance = 0;
    bool valid = true;
    const int n = 128;
    for (int y = 0; y < n; ++y) {
      for (int x = 0; x < n; ++x) {
        vec3f wi = disk.Sample(point3f(2), vec2f((x + .5f) / n, (y + .5f) / n), 0);
        Float pdf = disk.Pdf(point3f(2), wi, 0);
        valid &= std::isfinite(pdf) && pdf > 0;
        if (pdf > 0) irradiance += disk.Radiance(point3f(2), wi, 0)[0] * wi[1] / pdf;
      }
    }
    expect_true(valid);
    expect_true(std::abs(irradiance / (n * n) / (M_PI * std::pow(std::sin(radius), 2)) - 1) < 1e-5);
    expect_true(std::abs(disk.SamplingWeight() / (4 * M_PI * std::pow(std::sin(radius / 2), 2)) - 1) < .005);
    DiskInfiniteLight tiny(white, 16, 16, vec3f(0, 1, 0), .0001, 0);
    vec3f wi = tiny.Sample(point3f(0), vec2f(.5, .8), 0);
    expect_true((std::isfinite(tiny.Pdf(point3f(0), wi, 0)) && tiny.Pdf(point3f(0), wi, 0) > 0));
  }
  test_that("disks preserve texture orientation, rotations, and finite angular support") {
    auto tex = std::make_shared<DirectionalTexture>();
    DiskInfiniteLight disk(tex, 64, 64, vec3f(0, 0, 1), 20, 0, false);
    double t = std::tan(10 * M_PI / 180);
    vec3f wi(-.5 * t, .25 * t, 1);
    point3f expected = tex->value(.75, .625, point3f(0));
    expect_true((disk.Radiance(point3f(0), wi, 0) - expected).length() < 1e-6);
    expect_true(disk.Pdf(point3f(0), vec3f(1, 0, 0), 0) == 0);
    expect_true(disk.Radiance(point3f(0), vec3f(0, 0, -1), 0).length() == 0);
    DiskInfiniteLight rotated(tex, 64, 64, vec3f(0, 0, 1), 20, 47, false);
    expect_true((rotated.Radiance(point3f(100), RotateY(47)(wi), 1) - expected).length() < 1e-6);
    expect_true(rotated.Pdf(point3f(100), RotateY(47)(wi), 1) == disk.Pdf(point3f(0), wi, 0));
    DiskInfiniteLight horizon(tex, 64, 64, vec3f(0, 0, 1), 20, 0);
    expect_true(horizon.Radiance(point3f(0), vec3f(0, -.01, 1), 0).length() == 0);
    expect_true(horizon.Pdf(point3f(0), vec3f(0, -.01, 1), 0) > 0);
    DiskInfiniteLight below(tex, 64, 64, vec3f(0, -1, 0), .53, 0);
    expect_true(below.SamplingWeight() == 0);
    expect_true(below.Radiance(point3f(0), vec3f(0, -1, 0), 0).length() == 0);
  }
  test_that("mixed sky and disk proposals use compatible solid-angle power weights") {
    auto white = std::make_shared<constant_texture>(point3f(1));
    auto sky = std::make_shared<ImageInfiniteLight>(white, 128, 64, 0);
    auto disk = std::make_shared<DiskInfiniteLight>(white, 128, 128, vec3f(0, 1, 0), .53, 0);
    InfiniteLightMixture mixture({sky, disk});
    double probability = disk->SamplingWeight() / (sky->SamplingWeight() + disk->SamplingWeight());
    vec3f up(0, 1, 0);
    expect_true(std::abs(sky->SamplingWeight() / (4 * M_PI) - 1) < .001);
    expect_true(std::abs(mixture.Pdf(point3f(0), up, 0) -
        ((1 - probability) * sky->Pdf(point3f(0), up, 0) + probability * disk->Pdf(point3f(0), up, 0))) < 1e-6);
    expect_true(mixture.Radiance(point3f(0), up, 0)[0] == 2);
  }
}
#endif
