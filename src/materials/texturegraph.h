#ifndef RAYRENDER_TEXTUREGRAPH_H
#define RAYRENDER_TEXTUREGRAPH_H

#include "texture.h"
#include "texturecache.h"
#include <unordered_map>

// Immutable nodes have scalar or linear RGB output. Scalar nodes replicate their
// result in all three components; the compiler rejects implicit RGB reductions.
class TextureNode {
public:
  explicit TextureNode(bool scalar) : scalar(scalar) {}
  virtual ~TextureNode() = default;
  virtual point3f Evaluate(const TextureEvalContext&) const = 0;
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
