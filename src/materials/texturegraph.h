#ifndef RAYRENDER_TEXTUREGRAPH_H
#define RAYRENDER_TEXTUREGRAPH_H

#include "texture.h"
#include "texturecache.h"
#include <unordered_map>
#include <cstdint>

class TextureNode;
// Cold, renderer-independent description of an immutable node. Consumers can
// compile a graph without retaining R objects or changing its CPU evaluator.
struct TextureNodeDescription {
  enum Operation { Constant, Mix, Scale, Power, Channel, Direction, Checker, Gradient, Noise, Image };
  Operation operation = Constant;
  const TextureNode *children[3]{};
  point3f value{0}, scale{1}, offset{0};
  Float sine = 0, cosine = 1;
  uint32_t space = 0, argument = 0, seed = 0;
  bool flag = false;
  std::shared_ptr<const DecodedTextureImage> image;
};

// Immutable nodes have scalar or linear RGB output. Scalar nodes replicate their
// result in all three components; the compiler rejects implicit RGB reductions.
class TextureNode {
public:
  explicit TextureNode(bool scalar) : scalar(scalar) {}
  virtual ~TextureNode() = default;
  virtual point3f Evaluate(const TextureEvalContext&) const = 0;
  virtual TextureNodeDescription Describe() const;
  const bool scalar;
};

class TextureGraphBuilder {
public:
  explicit TextureGraphBuilder(TextureCache& images) : images(images) {}
  std::shared_ptr<const TextureNode> Build(const Rcpp::List& descriptor);
private:
  TextureCache& images;
  std::unordered_map<SEXP, std::shared_ptr<const TextureNode>> nodes;
  int depth = 0;
};

class graph_texture final : public texture {
public:
  const TextureNode& Root() const { return *root; }
  explicit graph_texture(std::shared_ptr<const TextureNode> root) : root(std::move(root)) {}
  point3f value(const hit_record& hit) const override { return value(TextureEvalContext::FromHit(hit)); }
  point3f value(const TextureEvalContext& context) const override { return root->Evaluate(context); }
  point3f value(Float u, Float v, const point3f& p) const override {
    TextureEvalContext context;
    context.u = u; context.v = v; context.p = context.object_p = p;
    return value(context);
  }
private:
  std::shared_ptr<const TextureNode> root;
};

#endif
