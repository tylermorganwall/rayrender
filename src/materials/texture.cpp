#include "../materials/texture.h"
#include <algorithm>
#include <cmath>

static constexpr Float rescale = static_cast<Float>(1)/static_cast<Float>(255);

namespace {

struct bilinear_texture_sample {
  int x0;
  int x1;
  int y0;
  int y1;
  Float x_weight;
  Float y_weight;
};

Float wrap_texture_coordinate(Float coordinate, Float repeat, Float offset = 0) {
  // Translate in texture space, after scaling; never pre-wrap the source UV.
  Float mapped = coordinate * repeat + offset;
  Float wrapped = std::fmod(mapped, static_cast<Float>(1));
  if(wrapped < 0) wrapped += 1;
  // Preserve the renderer's existing upper-endpoint convention.
  if(wrapped == 0 && mapped > 0) wrapped = 1;
  return wrapped;
}

bilinear_texture_sample get_bilinear_texture_sample(Float u, Float v,
                                                     Float repeatu, Float repeatv,
                                                     int nx, int ny, Float offsetu = 0, Float offsetv = 0) {
  u = wrap_texture_coordinate(u, repeatu, offsetu);
  v = wrap_texture_coordinate(v, repeatv, offsetv);
  Float x = u * static_cast<Float>(nx - 1);
  Float y = (1 - v) * static_cast<Float>(ny - 1);
  int x0 = static_cast<int>(std::floor(x));
  int y0 = static_cast<int>(std::floor(y));
  int x1 = std::min(x0 + 1, nx - 1);
  int y1 = std::min(y0 + 1, ny - 1);
  return {x0, x1, y0, y1, x - x0, y - y0};
}

Float bilinear_interpolate(Float value00, Float value10,
                           Float value01, Float value11,
                           Float x_weight, Float y_weight) {
  Float top = value00 * (1 - x_weight) + value10 * x_weight;
  Float bottom = value01 * (1 - x_weight) + value11 * x_weight;
  return top * (1 - y_weight) + bottom * y_weight;
}

} // namespace

point3f triangle_texture::value(Float u, Float v, const point3f& p) const {
    return(u * a + v * b + (1 - u - v) * c);
}

point3f image_texture_float::value(Float u, Float v, const point3f& p) const {
  bilinear_texture_sample sample = get_bilinear_texture_sample(
    u, v, repeatu, repeatv, nx, ny, offsetu, offsetv
  );
  auto channel_value = [&](int x, int y, int channel) {
    return data[channels*x + channels*nx*y + channel];
  };
  auto interpolate_channel = [&](int channel) {
    return intensity * bilinear_interpolate(
      channel_value(sample.x0, sample.y0, channel),
      channel_value(sample.x1, sample.y0, channel),
      channel_value(sample.x0, sample.y1, channel),
      channel_value(sample.x1, sample.y1, channel),
      sample.x_weight, sample.y_weight
    );
  };
  return point3f(
    interpolate_channel(0),
    interpolate_channel(1),
    interpolate_channel(2)
  );
}

point3f image_texture_char::value(Float u, Float v, const point3f& p) const {
  bilinear_texture_sample sample = get_bilinear_texture_sample(
    u, v, repeatu, repeatv, nx, ny, offsetu, offsetv
  );
  auto channel_value = [&](int x, int y, int channel) {
    Float value = static_cast<Float>(
      data[channels*x + channels*nx*y + channel]
    ) * rescale * intensity;
    return value * value;
  };
  auto interpolate_channel = [&](int channel) {
    return bilinear_interpolate(
      channel_value(sample.x0, sample.y0, channel),
      channel_value(sample.x1, sample.y0, channel),
      channel_value(sample.x0, sample.y1, channel),
      channel_value(sample.x1, sample.y1, channel),
      sample.x_weight, sample.y_weight
    );
  };
  return point3f(
    interpolate_channel(0),
    interpolate_channel(1),
    interpolate_channel(2)
  );
}

Float alpha_texture::value(Float u, Float v, const point3f& p) const {
  if (!data) return opacity;
  u = wrap_texture_coordinate(u, 1, offsetu);
  v = wrap_texture_coordinate(v, 1, offsetv);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  return(static_cast<Float>(data[channels*i + channels*nx*j + channels-1]) * rescale);
}


Float bump_texture::raw_value(Float u, Float v, const point3f& p) const {
  u = std::fmod(wrap_texture_coordinate(u, repeatu, offsetu), Float(1));
  v = std::fmod(wrap_texture_coordinate(v, repeatv, offsetv), Float(1));
  int i = u * (nx-1);
  int j = (1-v) * (ny-1);
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > nx-2) i = nx-2;
  if (j > ny-2) j = ny-2;
  return((Float)data[channels*i + channels*nx*j] * rescale);
}

point3f bump_texture::value(Float u, Float v, const point3f& p) const {
  u = std::fmod(wrap_texture_coordinate(u, repeatu, offsetu), Float(1));
  v = std::fmod(wrap_texture_coordinate(v, repeatv, offsetv), Float(1));
  int i = u * (nx-1);
  int j = (1-v) * (ny-1);
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > nx-2) i = nx-2;
  if (j > ny-2) j = ny-2;
  Float bu = Float(data[channels*(i+1) + channels*nx*j] - data[channels*(i-1) + channels*nx*j])/2 * rescale;
  Float bv = Float(data[channels*i + channels*nx*(j+1)] - data[channels*i + channels*nx*(j-1)])/2 * rescale;
  return(point3f(intensity*bu,intensity*bv,0));
}

point2f roughness_texture::raw_value(Float u, Float v) const {
  u = wrap_texture_coordinate(u, 1, offsetu);
  v = wrap_texture_coordinate(v, 1, offsetv);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float x = (Float)data[channels*i + channels*nx*j] * rescale;
  Float y = channels > 1 ? (Float)data[channels*i + channels*nx*j+1] * rescale : x;
  return point2f(x, y);
}

point2f roughness_texture::value(Float u, Float v) const {
  const auto raw = raw_value(u, v);
  Float alphax = RoughnessToAlpha(raw[0]);
  Float alphay = RoughnessToAlpha(raw[1]);
  return(point2f(alphax * alphax, alphay * alphay));
}

Float roughness_texture::RoughnessToAlpha(Float roughness) {
  roughness = std::fmax(roughness, (Float)0.0001550155);
  Float x = std::log(roughness);
  return(1.62142f + 0.819955f * x + 0.1734f * x * x +
         0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
}

#ifdef NOT_CRAN
#include <testthat.h>

context("Image texture interpolation") {
  test_that("[UV offsets translate lookups after scaling without resampling]") {
    Float pixels[8 * 8 * 3];
    unsigned char bytes[8 * 8 * 3];
    for (int i = 0; i < 8 * 8 * 3; ++i) {
      pixels[i] = Float((i * 37) % 251) / 255;
      bytes[i] = (i * 37) % 251;
    }
    const point3f p(0, 0, 0);
    image_texture_float base(pixels, 8, 8, 3);
    image_texture_char base_char(bytes, 8, 8, 3);
    for (Float du : {Float(-.37), Float(0), Float(.23)}) {
      const Float dv = -.19;
      image_texture_float shifted(pixels, 8, 8, 3, 1.5, 2.3, 1, du, dv);
      image_texture_char shifted_char(bytes, 8, 8, 3, 1.5, 2.3, 1, du, dv);
      for (Float u : {Float(-2.13), Float(.17), Float(1.31)}) {
        const Float v = .41;
        const Float s = u * Float(1.5) + du;
        const Float t = v * Float(2.3) + dv;
        const auto expected = base.value(s - std::floor(s), t - std::floor(t), p);
        const auto actual = shifted.value(u, v, p);
        const auto expected_char = base_char.value(s - std::floor(s), t - std::floor(t), p);
        const auto actual_char = shifted_char.value(u, v, p);
        for (int c = 0; c < 3; ++c) {
          expect_true(actual[c] == Approx(expected[c]));
          expect_true(actual_char[c] == Approx(expected_char[c]));
        }
      }
      alpha_texture alpha(bytes, 8, 8, 3, du, dv), alpha_base(bytes, 8, 8, 3);
      roughness_texture rough(bytes, 8, 8, 3, du, dv), rough_base(bytes, 8, 8, 3);
      bump_texture bump(bytes, 8, 8, 3, 1, 1.5, 2.3, du, dv);
      bump_texture bump_base(bytes, 8, 8, 3, 1);
      expect_true(alpha.value(.17, .41, p) == Approx(alpha_base.value(.17 + du, .41 + dv, p)));
      const auto r = rough.raw_value(.17, .41);
      const auto r0 = rough_base.raw_value(.17 + du, .41 + dv);
      expect_true(r[0] == Approx(r0[0]));
      expect_true(r[1] == Approx(r0[1]));
      const Float s = Float(.17) * Float(1.5) + du, t = Float(.41) * Float(2.3) + dv;
      expect_true(bump.raw_value(.17, .41, p) == Approx(bump_base.raw_value(s, t, p)));
      const auto b = bump.value(.17, .41, p), b0 = bump_base.value(s, t, p);
      expect_true(b[0] == Approx(b0[0]));
      expect_true(b[1] == Approx(b0[1]));
    }
  }
  test_that("[floating-point textures interpolate between texels]") {
    Float pixels[] = {
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    };
    image_texture_float texture(pixels, 2, 2, 3);
    point3f center = texture.value(0.5, 0.5, point3f(0, 0, 0));
    expect_true(center.xyz.x == Approx(0.25));
    expect_true(center.xyz.y == Approx(0.25));
    expect_true(center.xyz.z == Approx(0.25));

    point3f quarter = texture.value(0.25, 0.75, point3f(0, 0, 0));
    expect_true(quarter.xyz.x == Approx(0.1875));
    expect_true(quarter.xyz.y == Approx(0.1875));
    expect_true(quarter.xyz.z == Approx(0.0625));
  }

  test_that("[floating-point texture endpoints do not wrap]") {
    Float pixels[] = {
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    };
    image_texture_float texture(pixels, 2, 2, 3);
    point3f top_right = texture.value(1, 1, point3f(0, 0, 0));
    expect_true(top_right.xyz.x == Approx(1));
    expect_true(top_right.xyz.y == Approx(0));
    expect_true(top_right.xyz.z == Approx(0));
  }

  test_that("[byte textures interpolate in linear color space]") {
    unsigned char pixels[] = {
      0, 0, 0,
      255, 0, 0,
      0, 255, 0,
      0, 0, 255
    };
    image_texture_char texture(pixels, 2, 2, 3);
    point3f center = texture.value(0.5, 0.5, point3f(0, 0, 0));
    expect_true(center.xyz.x == Approx(0.25));
    expect_true(center.xyz.y == Approx(0.25));
    expect_true(center.xyz.z == Approx(0.25));

    point3f top_right = texture.value(1, 1, point3f(0, 0, 0));
    expect_true(top_right.xyz.x == Approx(1));
    expect_true(top_right.xyz.y == Approx(0));
    expect_true(top_right.xyz.z == Approx(0));
  }
}
#endif
