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

Float wrap_texture_coordinate(Float coordinate, Float repeat) {
  while(coordinate < 0) coordinate += 1;
  while(coordinate > 1) coordinate -= 1;
  Float repeated = coordinate * repeat;
  Float wrapped = std::fmod(repeated, static_cast<Float>(1));
  if(wrapped < 0) wrapped += 1;
  if(wrapped == 0 && repeated > 0) wrapped = 1;
  return wrapped;
}

bilinear_texture_sample get_bilinear_texture_sample(Float u, Float v,
                                                     Float repeatu, Float repeatv,
                                                     int nx, int ny) {
  u = wrap_texture_coordinate(u, repeatu);
  v = wrap_texture_coordinate(v, repeatv);
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
    u, v, repeatu, repeatv, nx, ny
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
    u, v, repeatu, repeatv, nx, ny
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
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  return(static_cast<Float>(data[channels*i + channels*nx*j + channels-1]) * rescale);
}


Float bump_texture::raw_value(Float u, Float v, const point3f& p) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  u = fmod(u * repeatu,1);
  v = fmod(v * repeatv,1);
  int i = u * (nx-1);
  int j = (1-v) * (ny-1);
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > nx-2) i = nx-2;
  if (j > ny-2) j = ny-2;
  return((Float)data[channels*i + channels*nx*j] * rescale);
}

point3f bump_texture::value(Float u, Float v, const point3f& p) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  u = fmod(u * repeatu,1);
  v = fmod(v * repeatv,1);
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

point2f roughness_texture::value(Float u, Float v) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float alphax = RoughnessToAlpha((Float)data[channels*i + channels*nx*j] * rescale);
  Float alphay = channels > 1 ? RoughnessToAlpha((Float)data[channels*i + channels*nx*j+1] * rescale) : alphax;
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
