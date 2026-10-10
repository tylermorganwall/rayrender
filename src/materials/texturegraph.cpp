#include "texturegraph.h"
#include "ptextexture.h"
#include "../hitables/hitable.h"
#include <algorithm>
#include <cmath>
#include <cstdint>

point3f texture::value(const hit_record& hit) const {
  // Existing/constant textures do not need to construct a geometric context.
  return value(hit.u, hit.v, hit.p);
}

TextureEvalContext TextureEvalContext::FromHit(const hit_record& h) {
  TextureEvalContext result;
  result.p = h.p; result.u = h.u; result.v = h.v;
  result.geometric_normal = h.geometric_normal.squared_length() > 0 ? h.geometric_normal : h.normal;
  if (result.geometric_normal.squared_length() > 0)
    result.geometric_normal = unit_vector(result.geometric_normal);
  result.dpdx=h.dpdx; result.dpdy=h.dpdy;
  result.has_derivatives=h.has_differentials;
  result.footprint={h.dudx,h.dvdx,h.dudy,h.dvdy,h.has_differentials};
  result.object_p = h.texture_object_p;
  result.object_normal = h.texture_object_normal;
  result.face_index = h.shape ? h.shape->TextureFaceIndex() : -1;
  return result;
}

TextureNodeDescription TextureNode::Describe() const {
  throw std::runtime_error("Ptex or unsupported texture graph node");
}

namespace {
using Node = std::shared_ptr<const TextureNode>;

struct Mapping {
  std::string space;
  point3f scale, offset;
  Float sine, cosine;
  explicit Mapping(const Rcpp::List& p) {
    space = Rcpp::as<std::string>(p["space"]);
    if (space != "uv" && space != "object" && space != "world") Rcpp::stop("Invalid texture coordinate space.");
    Rcpp::NumericVector s = p["scale"], o = p["offset"];
    if (s.size() != 3 || o.size() != 3) Rcpp::stop("Texture mapping requires three-component scale and offset.");
    for (int c = 0; c < 3; ++c) {
      if (!std::isfinite(s[c]) || !std::isfinite(o[c])) Rcpp::stop("Nonfinite texture mapping.");
      scale[c] = s[c]; offset[c] = o[c];
    }
    double angle = Rcpp::as<double>(p["rotation"]);
    if (!std::isfinite(angle)) Rcpp::stop("Nonfinite texture rotation.");
    sine = std::sin(angle); cosine = std::cos(angle);
  }
void Describe(TextureNodeDescription& out) const {
  out.space = space == "uv" ? 0 : space == "object" ? 1 : 2;
  out.scale = scale; out.offset = offset; out.sine = sine; out.cosine = cosine;
}
point3f Map(const TextureEvalContext& h) const {

    point3f p = space == "uv" ? point3f(h.u, h.v, 0) : space == "world" ? h.p : h.object_p;
    p = p * scale;
    return point3f(cosine*p[0] - sine*p[1], sine*p[0] + cosine*p[1], p[2]) + offset;
  }
};

class Constant final : public TextureNode {
  point3f value;
public:
  Constant(bool scalar, point3f value) : TextureNode(scalar), value(value) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Constant;
    out.value = value;
    return out;
  }
  point3f Evaluate(const TextureEvalContext&) const override { return value; }
};

class PtexNode final : public TextureNode {
  std::shared_ptr<const PtexTextureResource> resource;
  int filter;
  std::string encoding;
  Float gamma;
  bool quantize;
public:
  PtexNode(const Rcpp::List& p, TextureCache& cache, bool scalar)
    : TextureNode(scalar), resource(cache.LookupPtex(Rcpp::as<std::string>(p["filename"]))),
      filter(Rcpp::as<int>(p["filter"])), encoding(Rcpp::as<std::string>(p["encoding"])),
      gamma(Rcpp::as<Float>(p["gamma"])), quantize(Rcpp::as<bool>(p["pbrt_encoding"])) {
    if (filter < PTEX_POINT || filter > PTEX_MITCHELL ||
        (encoding != "linear" && encoding != "srgb" && encoding != "gamma") ||
        !std::isfinite(gamma) || gamma <= 0)
      Rcpp::stop("Invalid Ptex filter or encoding.");
  }
  point3f Evaluate(const TextureEvalContext& context) const override {
    point3f result = resource->Sample(context, filter);
    if (encoding != "linear") {
      for (int c = 0; c < 3; ++c) {
        Float x = result[c];
        // PBRT v4 decodes filtered Ptex through an 8-bit color-encoding table.
        // The public constructor keeps full precision; the importer opts in.
        if (quantize) x = std::floor(std::clamp(x * 255 + Float(.5), Float(0), Float(255))) / 255;
        result[c] = encoding == "gamma" ? std::pow(std::max(Float(0), x), gamma) :
          x <= Float(.04045) ? x / Float(12.92) : std::pow((x + Float(.055)) / Float(1.055), Float(2.4));
      }
    }
    return scalar ? point3f((result[0] + result[1] + result[2]) / 3) : result;
  }
};

class Mix final : public TextureNode {
  Node a, b, weight;
public:
  Mix(Node a, Node b, Node weight) : TextureNode(a->scalar && b->scalar), a(a), b(b), weight(weight) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Mix;
    out.children[0] = a.get();
    out.children[1] = b.get();
    out.children[2] = weight.get();
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    Float w = std::clamp(weight->Evaluate(h)[0], Float(0), Float(1));
    if (w == 0) return a->Evaluate(h);
    if (w == 1) return b->Evaluate(h);
    return (1-w)*a->Evaluate(h) + w*b->Evaluate(h);
  }
};

class Direction final : public TextureNode {
  normal3f direction;
  bool absolute, object;
public:
  Direction(normal3f direction, bool absolute, bool object) : TextureNode(true), direction(unit_vector(direction)), absolute(absolute), object(object) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Direction;
    out.value = point3f(direction[0], direction[1], direction[2]);
    out.flag = absolute;
    out.space = object ? 1 : 2;
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    normal3f n = object ? h.object_normal : h.geometric_normal;
    Float weight = n.squared_length() > 0 ? dot(unit_vector(n), direction) : 0;
    return point3f(absolute ? std::abs(weight) : std::max(Float(0), weight));
  }
};

class ScaledTexture final : public TextureNode {
  Node child, factor;
public:
  ScaledTexture(Node child, Node factor) : TextureNode(child->scalar), child(child), factor(factor) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Scale;
    out.children[0] = child.get();
    out.children[1] = factor.get();
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override { return child->Evaluate(h) * factor->Evaluate(h)[0]; }
};

class Channel final : public TextureNode {
  Node child;
  int channel;
public:
  Channel(Node child, int channel) : TextureNode(true), child(child), channel(channel) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Channel;
    out.children[0] = child.get();
    out.argument = channel;
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    point3f c = child->Evaluate(h);
    return point3f(channel < 3 ? c[channel] : channel == 3 ? (c[0]+c[1]+c[2])/3 : .2126f*c[0]+.7152f*c[1]+.0722f*c[2]);
  }
};

// Import adapter: PBRT remaps roughness only after its texture is evaluated.
class Power final : public TextureNode {
  Node child;
  Float exponent;
public:
  Power(Node child, Float exponent) : TextureNode(true), child(child), exponent(exponent) {}
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    out.operation = TextureNodeDescription::Power;
    out.children[0] = child.get();
    out.value = point3f(exponent);
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    return point3f(std::pow(std::max(Float(0), child->Evaluate(h)[0]), exponent));
  }
};

class Pattern final : public TextureNode {
  Mapping mapping;
  Node a, b;
  std::string op;
  int axis, octaves;
  uint32_t seed;
  static uint32_t Hash(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; return x ^ (x >> 16);
  }
  Float Noise(point3f p) const {
    // Seeded smooth value noise. Reduce before integer conversion to avoid
    // overflow for large world coordinates; the lattice repeats at 2^20 cells.
    uint32_t cell[3]; Float f[3];
    for (int i = 0; i < 3; ++i) {
      double base = std::floor(double(p[i]));
      cell[i] = uint32_t(std::fmod(std::fmod(base, 1048576.) + 1048576., 1048576.));
      Float t = p[i] - base;
      f[i] = t*t*t*(t*(t*6-15)+10);
    }
    Float result = 0;
    for (int z=0; z<2; ++z) for (int y=0; y<2; ++y) for (int x=0; x<2; ++x) {
      // Wrap each corner, not just the cell origin. Otherwise the upper corner
      // of cell -1 hashes as 2^20 while the lower corner of cell 0 hashes as 0,
      // creating a seam at the coordinate planes in every octave.
      uint32_t ix = (cell[0] + x) & 0xfffffU;
      uint32_t iy = (cell[1] + y) & 0xfffffU;
      uint32_t iz = (cell[2] + z) & 0xfffffU;
      uint32_t h = Hash(seed ^ Hash(ix) ^ Hash(Hash(iy)) ^ Hash(Hash(Hash(iz))));
      result += Float(h >> 8) / Float(0xffffffU) * (x ? f[0] : 1-f[0]) * (y ? f[1] : 1-f[1]) * (z ? f[2] : 1-f[2]);
    }
    return result;
  }
public:
  Pattern(const Rcpp::List& p, Node a = nullptr, Node b = nullptr) : TextureNode(!a || (a->scalar && b->scalar)),
      mapping(Rcpp::as<Rcpp::List>(p["coordinates"])), a(a), b(b), op(Rcpp::as<std::string>(p["op"])), axis(0), octaves(1), seed(0) {
    if (op == "gradient") {
      axis = Rcpp::as<int>(p["axis"]);
      if (axis < 0 || axis > 2) Rcpp::stop("Invalid gradient axis.");
    }
    if (op == "noise") {
      octaves = Rcpp::as<int>(p["octaves"]); seed = uint32_t(Rcpp::as<int>(p["seed"]));
      if (octaves < 1 || octaves > 16) Rcpp::stop("Invalid noise octave count.");
    }
  }
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    mapping.Describe(out);
    out.operation = op == "checker"    ? TextureNodeDescription::Checker
                    : op == "gradient" ? TextureNodeDescription::Gradient
                                       : TextureNodeDescription::Noise;
    out.children[0] = a.get();
    out.children[1] = b.get();
    out.argument = op == "noise" ? octaves : axis;
    out.seed = seed;
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    point3f p = mapping.Map(h);
    if (op == "gradient") return point3f(std::clamp(p[axis], Float(0), Float(1)));
    if (op == "checker") {
      double parity = std::fmod(std::floor(double(p[0])) + std::floor(double(p[1])) + std::floor(double(p[2])), 2.);
      return parity == 0 ? a->Evaluate(h) : b->Evaluate(h);
    }
    Float sum = 0, amplitude = 1, total = 0;
    for (int i=0; i<octaves; ++i) { sum += amplitude*Noise(p); total += amplitude; amplitude *= .5f; p *= 2; }
    return point3f(sum/total);
  }
};

class Image final : public TextureNode {
  Mapping mapping;
  int width, height;
  bool repeat;
  std::shared_ptr<const DecodedTextureImage> image;
public:
  Image(const Rcpp::List& p, TextureCache& cache, bool scalar) : TextureNode(scalar),
      mapping(Rcpp::as<Rcpp::List>(p["coordinates"])) {
    std::string file = Rcpp::as<std::string>(p["filename"]);
    std::string encoding = Rcpp::as<std::string>(p["encoding"]), wrap = Rcpp::as<std::string>(p["wrap"]);
    if (encoding != "linear" && encoding != "srgb") Rcpp::stop("Invalid texture encoding.");
    if (wrap != "repeat" && wrap != "clamp") Rcpp::stop("Invalid image wrap mode.");
    repeat = wrap == "repeat";
    image = cache.LookupGraphImage(file, encoding);
    width = image->width; height = image->height;
  }
  TextureNodeDescription Describe() const override {
    TextureNodeDescription out;
    mapping.Describe(out);
    out.operation = TextureNodeDescription::Image;
    out.image = image;
    out.flag = repeat;
    return out;
  }
  point3f Evaluate(const TextureEvalContext& h) const override {
    point3f p = mapping.Map(h);
    Float u = repeat ? p[0]-std::floor(p[0]) : std::clamp(p[0], Float(0), Float(1));
    Float v = repeat ? p[1]-std::floor(p[1]) : std::clamp(p[1], Float(0), Float(1));
    Float x = u*width-.5f, y = (1-v)*height-.5f;
    int ix = int(std::floor(x)), iy = int(std::floor(y));
    auto fetch = [&](int i, int j) {
      i = repeat ? (i%width+width)%width : std::clamp(i, 0, width-1);
      j = repeat ? (j%height+height)%height : std::clamp(j, 0, height-1);
      return image->pixels[size_t(j)*width+i];
    };
    Float tx = x-ix, ty = y-iy;
    point3f value = (1-ty)*((1-tx)*fetch(ix,iy)+tx*fetch(ix+1,iy)) + ty*((1-tx)*fetch(ix,iy+1)+tx*fetch(ix+1,iy+1));
    return scalar ? point3f((value[0]+value[1]+value[2])/3) : value;
  }
};
}

std::shared_ptr<const TextureNode> TextureGraphBuilder::Build(const Rcpp::List& p) {
  auto found = nodes.find(p);
  if (found != nodes.end()) return found->second;
  if (++depth > 128) Rcpp::stop("Texture graph exceeds 128 levels or contains a cycle.");
  std::string op = Rcpp::as<std::string>(p["op"]), type = Rcpp::as<std::string>(p["type"]);
  if (type != "scalar" && type != "color") Rcpp::stop("Invalid texture output type.");
  bool scalar = type == "scalar";
  auto child = [&](const char* field) { return Build(Rcpp::as<Rcpp::List>(p[field])); };
  auto scalar_child = [&](const char* field) {
    auto n = child(field);
    if (!n->scalar) Rcpp::stop("Texture input requires a scalar; use texture_channel().");
    return n;
  };
  Node node;
  if (op == "constant") {
    Rcpp::NumericVector v = p["value"];
    if (v.size() != (scalar ? 1 : 3)) Rcpp::stop("Invalid texture constant size.");
    point3f value;
    for (int c=0; c<3; ++c) {
      double x = v[scalar ? 0 : c];
      if (!std::isfinite(x) || !std::isfinite(Float(x))) Rcpp::stop("Nonfinite texture constant.");
      value[c] = x;
    }
    node = std::make_shared<Constant>(scalar, value);
  } else if (op == "mix") node = std::make_shared<Mix>(child("a"), child("b"), scalar_child("weight"));
  else if (op == "scale") node = std::make_shared<ScaledTexture>(child("child"), scalar_child("factor"));
  else if (op == "power") {
    Float exponent = Rcpp::as<Float>(p["exponent"]);
    if (!std::isfinite(exponent) || exponent <= 0) Rcpp::stop("Invalid texture exponent.");
    node = std::make_shared<Power>(scalar_child("child"), exponent);
  }
  else if (op == "channel") {
    int channel = Rcpp::as<int>(p["channel"]);
    if (channel < 0 || channel > 4) Rcpp::stop("Invalid texture channel.");
    node = std::make_shared<Channel>(child("child"), channel);
  } else if (op == "direction") {
    Rcpp::NumericVector d = p["direction"];
    if (d.size()!=3) Rcpp::stop("Texture direction requires three components.");
    normal3f n(d[0], d[1], d[2]);
    if (!std::isfinite(n.squared_length()) || n.squared_length() == 0) Rcpp::stop("Invalid texture direction.");
    std::string space = Rcpp::as<std::string>(p["space"]);
    if (space != "object" && space != "world") Rcpp::stop("Invalid direction texture space.");
    node = std::make_shared<Direction>(n, Rcpp::as<bool>(p["absolute"]), space == "object");
  } else if (op == "checker") node = std::make_shared<Pattern>(p, child("a"), child("b"));
  else if (op == "gradient" || op == "noise") node = std::make_shared<Pattern>(p);
  else if (op == "image") node = std::make_shared<Image>(p, images, scalar);
  else if (op == "ptex") node = std::make_shared<PtexNode>(p, images, scalar);
  else Rcpp::stop("Unknown texture graph operation: " + op);
  if (node->scalar != scalar) Rcpp::stop("Texture graph output type does not match its children.");
  nodes.emplace(p, node);
  --depth;
  return node;
}

#ifdef NOT_CRAN
#include "microfacetdist.h"
#include <testthat.h>

context("Composable texture graphs") {
  test_that("direction mixing uses geometric orientation and remaps after blending") {
    auto a = std::make_shared<Constant>(true, point3f(.0005f));
    auto b = std::make_shared<Constant>(true, point3f(.005f));
    auto direction = std::make_shared<Direction>(normal3f(0, 3, 0), true, false);
    auto mixed = std::make_shared<Mix>(b, a, direction);
    TextureEvalContext h;
    expect_true(std::abs(mixed->Evaluate(h)[0] - .0005f) < 1e-8f);
    h.geometric_normal = normal3f(0, -1, 0);
    expect_true(std::abs(mixed->Evaluate(h)[0] - .0005f) < 1e-8f);
    h.geometric_normal = normal3f(1, 0, 0);
    expect_true(std::abs(mixed->Evaluate(h)[0] - .005f) < 1e-8f);
    h.geometric_normal = normal3f(std::sqrt(.75f), .5f, 0);
    Power remapped(mixed, .5f);
    expect_true(std::abs(remapped.Evaluate(h)[0] - std::sqrt(.00275f)) < 1e-6f);
    // Blending two already remapped endpoints is a different distribution.
    expect_true(std::abs(remapped.Evaluate(h)[0] - .5f*(std::sqrt(.0005f)+std::sqrt(.005f))) > .005f);
    Direction upward(normal3f(0, 1, 0), false, false);
    h.geometric_normal = normal3f(0, -1, 0);
    expect_true(upward.Evaluate(h)[0] == 0);
    h.object_normal = normal3f(0, 1, 0);
    Direction object(normal3f(0, 1, 0), false, true);
    expect_true(object.Evaluate(h)[0] == 1);
  }

  test_that("native graph builder validates types and preserves shared children") {
    TextureCache images;
    TextureGraphBuilder builder(images);
    auto one = Rcpp::List::create(Rcpp::_["op"]="constant", Rcpp::_["type"]="scalar", Rcpp::_["value"]=Rcpp::NumericVector::create(.25));
    auto first = builder.Build(one);
    expect_true(first == builder.Build(one));
    auto color = std::make_shared<Constant>(false, point3f(1, .5, 0));
    auto weight = std::make_shared<Constant>(true, point3f(.5));
    Mix blended(first, color, weight);
    TextureEvalContext h;
    auto result = blended.Evaluate(h);
    expect_true(std::abs(result[0]-.625f) < 1e-7f);
    expect_true(std::abs(result[1]-.375f) < 1e-7f);
    Channel channel(color, 3);
    expect_true(channel.Evaluate(h)[0] == .5f);
  }

  test_that("coordinate mappings and seeded noise are independent of hit world placement") {
    auto mapping = Rcpp::List::create(Rcpp::_["space"]="object", Rcpp::_["scale"]=Rcpp::NumericVector::create(2, 3, 1),
      Rcpp::_["offset"]=Rcpp::NumericVector::create(.1, .2, 0), Rcpp::_["rotation"]=0.);
    TextureEvalContext h;
    h.object_p = point3f(.2, .3, 0); h.p = point3f(100, 100, 0);
    Mapping coordinates(mapping);
    expect_true(std::abs(coordinates.Map(h)[0]-.5f) < 1e-6f);
    expect_true(std::abs(coordinates.Map(h)[1]-1.1f) < 1e-6f);
    auto descriptor = Rcpp::List::create(Rcpp::_["op"]="noise", Rcpp::_["coordinates"]=mapping,
      Rcpp::_["octaves"]=4, Rcpp::_["seed"]=12);
    Pattern noise(descriptor);
    Float value = noise.Evaluate(h)[0];
    h.p = point3f(-200, 50, 10);
    expect_true(noise.Evaluate(h)[0] == value);
    expect_true((value >= 0 && value <= 1));
    descriptor["seed"] = 13;
    Pattern other(descriptor);
    expect_true(other.Evaluate(h)[0] != value);
  }

  test_that("public graph roughness matches numeric roughness and PBRT alpha bypasses it") {
    hit_record h;
    h.u = h.v = .5; h.p = point3f(0); h.normal = normal3f(0, 1, 0);
    TrowbridgeReitzDistribution numeric(.3f*.3f, .3f*.3f, nullptr, false);
    TrowbridgeReitzDistribution graph(.1, .1, nullptr, false);
    graph.roughness_graph = std::make_shared<Constant>(true, point3f(.3f));
    expect_true(std::abs(numeric.GetAlphas(.5,.5)[0] - graph.Resolve(h)[0]) < 1e-6f);
    graph.graph_is_alpha = true;
    expect_true(graph.Resolve(h)[0] == .3f);
    graph.roughness_graph_v = std::make_shared<Constant>(true, point3f(.1f));
    expect_true(graph.Resolve(h)[1] == .1f);
  }

  test_that("noise remains continuous across positive and negative lattice boundaries") {
    auto mapping = Rcpp::List::create(Rcpp::_["space"]="object",
      Rcpp::_["scale"]=Rcpp::NumericVector::create(1, 1, 1),
      Rcpp::_["offset"]=Rcpp::NumericVector::create(0, 0, 0), Rcpp::_["rotation"]=0.);
    for (int octaves : {1, 5}) {
      auto descriptor = Rcpp::List::create(Rcpp::_["op"]="noise", Rcpp::_["coordinates"]=mapping,
        Rcpp::_["octaves"]=octaves, Rcpp::_["seed"]=14);
      Pattern noise(descriptor);
      for (int axis=0; axis<3; ++axis) for (Float boundary : {-1.f, 0.f, 1.f}) {
        TextureEvalContext left, right;
        left.object_p = right.object_p = point3f(.37f, .63f, .29f);
        left.object_p[axis] = boundary - 1e-6f;
        right.object_p[axis] = boundary + 1e-6f;
        expect_true(std::abs(noise.Evaluate(left)[0] - noise.Evaluate(right)[0]) < 1e-4f);
      }
    }
  }
}
#endif
