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

point3f interpolate_float_texture(const Float *data, int nx, int channels, Float intensity,
                                 const bilinear_texture_sample &sample) {
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

} // namespace

point3f triangle_texture::value(Float u, Float v, const point3f& p) const {
  return(u * a + v * b + (1 - u - v) * c);
}

point3f image_texture_float::value(Float u, Float v, const point3f& p) const {
  const auto sample = get_bilinear_texture_sample(u, v, repeatu, repeatv, nx, ny, offsetu, offsetv);
  return interpolate_float_texture(data, nx, channels, intensity, sample);
}

point3f latlong_image_texture::value(Float u, Float v, const point3f&) const {
  // Each texel spans 1/n of its axis. Interpolate across longitude's seam:
  // u=0 lies halfway between the last and first texel centers. Endpoint-based
  // lookup would halve the support of both edge texels and dim an HDR sun.
  const Float x = (u - std::floor(u)) * nx - .5f;
  const Float y = clamp((1 - v) * ny - .5f, Float(0), Float(ny - 1));
  const int ix = int(std::floor(x)), iy = int(std::floor(y));
  const int x0 = (ix + nx) % nx;
  const bilinear_texture_sample sample{x0, (x0 + 1) % nx, iy, std::min(iy + 1, ny - 1),
                                      x - ix, y - iy};
  return interpolate_float_texture(data, nx, channels, intensity, sample);
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
  u = wrap_texture_coordinate(u, repeatu, offsetu);
  v = wrap_texture_coordinate(v, repeatv, offsetv);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  return(opacity * static_cast<Float>(data[channels*i + channels*nx*j + channels-1]) * rescale);
}

HeightImage::HeightImage(int width, int height, std::vector<Float> pixels) {
  if (width < 1 || height < 1 || pixels.size() != size_t(width) * height)
    throw std::invalid_argument("Invalid height image dimensions");
  // PBRT-style power-of-two pyramid. Upsample non-power-of-two dimensions with
  // a separable radius-two windowed sinc, then box-filter successive levels.
  // Heights are signed data: unlike color resampling, never clamp them to zero.
  auto power2 = [](int n) {
    int v = 1;
    while (v < n)
      v *= 2;
    return v;
  };
  const int out_width = power2(width), out_height = power2(height);
  auto resize_axis = [](const std::vector<Float> &source, int w, int h, int target,
                        bool horizontal) {
    const int old = horizontal ? w : h, nw = horizontal ? target : w, nh = horizontal ? h : target;
    std::vector<Float> result(size_t(nw) * nh);
    for (int j = 0; j < nh; ++j)
      for (int i = 0; i < nw; ++i) {
        const double center = ((horizontal ? i : j) + .5) * old / target;
        const int first = int(std::floor(center - 1.5));
        double sum = 0, weight_sum = 0;
        for (int k = 0; k < 4; ++k) {
          const double d = first + k + .5 - center;
          auto sinc = [](double x) {
            return std::abs(x) < 1e-8 ? 1.0 : std::sin(M_PI * x) / (M_PI * x);
          };
          const double weight = std::abs(d) > 2 ? 0 : sinc(d) * sinc(d / 2);
          const int tap = ((first + k) % old + old) % old;
          sum += weight * source[horizontal ? size_t(j) * w + tap : size_t(tap) * w + i];
          weight_sum += weight;
        }
        result[size_t(j) * nw + i] = Float(sum / weight_sum);
      }
    return result;
  };
  if (out_width != width) {
    pixels = resize_axis(pixels, width, height, out_width, true);
    width = out_width;
  }
  if (out_height != height) {
    pixels = resize_axis(pixels, width, height, out_height, false);
    height = out_height;
  }
  levels.push_back({width, height, std::move(pixels)});
  while (width > 1 || height > 1) {
    const auto &previous = levels.back();
    const int nw = std::max(1, width / 2), nh = std::max(1, height / 2);
    std::vector<Float> next(size_t(nw) * nh);
    for (int j = 0; j < nh; ++j)
      for (int i = 0; i < nw; ++i) {
        const int x = 2 * i, y = 2 * j, x1 = std::min(x + 1, width - 1),
                  y1 = std::min(y + 1, height - 1);
        next[size_t(j) * nw + i] =
            Float(.25) *
            (previous.pixels[size_t(y) * width + x] + previous.pixels[size_t(y) * width + x1] +
             previous.pixels[size_t(y1) * width + x] + previous.pixels[size_t(y1) * width + x1]);
      }
    levels.push_back({nw, nh, std::move(next)});
    width = nw;
    height = nh;
  }
}

Float bump_texture::raw_value(Float u, Float v, const point3f &p, TextureFootprint footprint) const {
  if (height_texture) {
    TextureEvalContext context;
    context.p = context.object_p = p;
    context.u = u; context.v = v; context.footprint = footprint;
    return height_texture->value(context)[0];
  }
  const HeightLevel *level = nullptr;
  int width = nx, height = ny;
  if (image) {
    size_t lod = 0;
    if (footprint.valid) {
      const Float extent =
          2 * std::max({std::abs(repeatu * footprint.dudx), std::abs(repeatu * footprint.dudy),
                        std::abs(repeatv * footprint.dvdx), std::abs(repeatv * footprint.dvdy)});
      const Float l = Float(image->levels.size() - 1) + std::log2(std::max(extent, Float(1e-8)));
      lod = size_t(std::clamp(std::floor(l), Float(0), Float(image->levels.size() - 1)));
    }
    level = &image->levels[lod];
    width = level->width;
    height = level->height;
  }
  // PBRT's bilinear texel-center convention, with repeat addressing applied to
  // each tap. Interpolate across periodic seams rather than clamp at the edge.
  const Float s = u * repeatu + offsetu, t = v * repeatv + offsetv;
  const Float x = (s - std::floor(s)) * width - Float(.5);
  const Float y = (1 - (t - std::floor(t))) * height - Float(.5);
  const int ix = int(std::floor(x)), iy = int(std::floor(y));
  const Float fx = x - ix, fy = y - iy;
  auto texel = [&](int i, int j) {
    i = (i % width + width) % width;
    j = (j % height + height) % height;
    const size_t index = size_t(i) + size_t(width) * j;
    return level ? level->pixels[index] : Float(data[size_t(channels) * index]) * rescale;
  };
  const Float a = (1 - fx) * texel(ix, iy) + fx * texel(ix + 1, iy);
  const Float b = (1 - fx) * texel(ix, iy + 1) + fx * texel(ix + 1, iy + 1);
  return (1 - fy) * a + fy * b;
}

point3f bump_texture::value(Float u, Float v, const point3f &p, TextureFootprint footprint) const {
  // PBRT v4 BumpMap: absent ray differentials, du=dv=.0005 in source UV
  // coordinates. Repeat's chain rule is implicit in evaluating shifted UVs.
  Float du =
      footprint.valid ? Float(.5) * (std::abs(footprint.dudx) + std::abs(footprint.dudy)) : 0;
  Float dv =
      footprint.valid ? Float(.5) * (std::abs(footprint.dvdx) + std::abs(footprint.dvdy)) : 0;
  if (du == 0)
    du = .0005;
  if (dv == 0)
    dv = .0005;
  const Float height = raw_value(u, v, p, footprint);
  const Float bu = (raw_value(u + du, v, p, footprint) - height) / du;
  const Float bv = (raw_value(u, v + dv, p, footprint) - height) / dv;
  // Keep the established image-down convention at the legacy slope API.
  return point3f(intensity * bu, -intensity * bv, 0);
}

normal3f bump_texture::perturb(Float u, Float v, const point3f &p, normal3f n, vec3f &dpdu,
                               vec3f &dpdv, normal3f dndu, normal3f dndv,
                               TextureFootprint footprint, const TextureEvalContext* context) const {
  point3f slope;
  Float height;
  if (height_texture) {
    TextureEvalContext h = context ? *context : TextureEvalContext();
    h.u = u; h.v = v; h.p = p; h.footprint = footprint;
    h.geometric_normal = n;
    Float du = footprint.valid ? Float(.5) * (std::abs(footprint.dudx) + std::abs(footprint.dudy)) : 0;
    Float dv = footprint.valid ? Float(.5) * (std::abs(footprint.dvdx) + std::abs(footprint.dvdy)) : 0;
    if (du == 0) du = .0005;
    if (dv == 0) dv = .0005;
    height = intensity * height_texture->value(h)[0];
    auto shifted_u = h, shifted_v = h;
    shifted_u.u += du; shifted_u.p += du * dpdu;
    shifted_v.v += dv; shifted_v.p += dv * dpdv;
    shifted_u.geometric_normal = unit_vector(n + du * dndu);
    shifted_v.geometric_normal = unit_vector(n + dv * dndv);
    slope = point3f((intensity * height_texture->value(shifted_u)[0] - height) / du,
                   -(intensity * height_texture->value(shifted_v)[0] - height) / dv, 0);
  } else {
    slope = value(u, v, p, footprint);
    height = intensity * raw_value(u, v, p, footprint);
  }
  const vec3f du = dpdu + slope[0] * convert_to_vec3(n) + height * convert_to_vec3(dndu);
  const vec3f dv = dpdv - slope[1] * convert_to_vec3(n) + height * convert_to_vec3(dndv);
  const auto crossed = cross(du, dv);
  if (!(crossed.squared_length() > 0) || !std::isfinite(crossed.squared_length()))
    return n;
  auto bumped = convert_to_normal3(unit_vector(crossed));
  // Preserve orientation for mirrored UV charts and reverseOrientation.
  bumped = Faceforward(bumped, n);
  dpdu = du;
  dpdv = dv;
  return bumped;
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
  test_that("[graph bump heights preserve face context and UV units]") {
    class FaceRamp final : public texture {
    public:
      point3f value(Float, Float, const point3f&) const override { return point3f(0); }
      point3f value(const TextureEvalContext& h) const override {
        return point3f(h.face_index == 7 ? Float(.3) * h.u : 0);
      }
    };
    bump_texture bump(std::make_shared<FaceRamp>(), 1);
    TextureEvalContext context;
    context.face_index = 7;
    context.footprint = {Float(.03), 0, 0, Float(.02), true};
    vec3f du(2, 0, 0), dv(0, 2, 0);
    const auto actual = bump.perturb(.5, .5, point3f(0), normal3f(0, 0, 1), du, dv,
                                     normal3f(0), normal3f(0), context.footprint, &context);
    const auto expected = unit_vector(normal3f(-.15, 0, 1));
    expect_true((actual - expected).length() < 1e-5);
    expect_true(std::abs(du[2] - Float(.3)) < 1e-5);
  }
  test_that("[filtered height derivatives use UV units and signed repeats]") {
    const point3f p(0);
    for (int size : {8, 32})
      for (int axis : {0, 1}) {
        std::vector<Float> pixels(size * size);
        for (int j = 0; j < size; ++j)
          for (int i = 0; i < size; ++i)
            pixels[i + size * j] = (Float(axis == 0 ? i : j) + Float(.5)) / size;
        auto image = std::make_shared<HeightImage>(size, size, pixels);
        for (Float repeat : {Float(-2), Float(.5), Float(3)}) {
          // Keep the lookup and forward difference inside the linear ramp.
          bump_texture bump(image, size, size, 1, .012, repeat, repeat, .5 - .4 * repeat,
                            .5 - .4 * repeat);
          auto slope = bump.value(.4, .4, p);
          expect_true(std::abs(slope[0] - (axis == 0 ? Float(.012) * repeat : 0)) < 1e-5);
          expect_true(std::abs(slope[1] - (axis == 1 ? Float(.012) * repeat : 0)) < 1e-5);
        }
      }
  }
  test_that("[height pyramids preserve signed constants and filter small features]") {
    for (int nx : {1, 2, 5})
      for (int ny : {1, 2, 7}) {
        auto image = std::make_shared<HeightImage>(nx, ny, std::vector<Float>(nx * ny, -.012345));
        bump_texture bump(image, nx, ny, 1, 1);
        for (Float uv : {Float(-1), Float(0), Float(.999), Float(2)}) {
          expect_true(bump.value(uv, uv, point3f(0)).squared_length() < 1e-8);
          expect_true(bump.raw_value(uv, uv, point3f(0)) == Approx(-.012345));
        }
      }
    std::vector<Float> pixels(64);
    for (int j = 0; j < 8; ++j)
      for (int i = 0; i < 8; ++i)
        pixels[i + 8 * j] = (i + j) % 2;
    auto image = std::make_shared<HeightImage>(8, 8, pixels);
    bump_texture bump(image, 8, 8, 1, 1);
    TextureFootprint footprint;
    footprint.valid = true;
    footprint.dudx = footprint.dvdy = 1;
    expect_true(bump.raw_value(.123, .456, point3f(0), footprint) == Approx(.5));
    expect_true(bump.value(.123, .456, point3f(0), footprint).squared_length() == 0);
  }
  test_that("[bump construction retains the normal derivative height term]") {
    auto image = std::make_shared<HeightImage>(1, 1, std::vector<Float>{.1});
    bump_texture bump(image, 1, 1, 1, 1);
    vec3f du(2, 0, 0), dv(0, 3, 0);
    auto n = bump.perturb(.5, .5, point3f(0), normal3f(0, 0, 1), du, dv, normal3f(1, 0, 0),
                          normal3f(0, 2, 0));
    expect_true(du[0] == Approx(2.1));
    expect_true(dv[1] == Approx(3.2));
    expect_true(n[2] == Approx(1));
  }
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
      // Float height subtraction at a .0005 UV step amplifies rounding by
      // about 2000. Bound that numerical error in the chain-rule comparison.
      expect_true(std::abs(b[0]-Float(1.5)*b0[0]) < 5e-4);
      expect_true(std::abs(b[1]-Float(2.3)*b0[1]) < 5e-4);
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
