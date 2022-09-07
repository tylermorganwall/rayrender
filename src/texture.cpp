#include "texture.h"

static constexpr Float rescale = 1.f/255.f;

point3f triangle_texture::value(Float u, Float v, const point3f& p) const {
    return(u * a + v * b + (1 - u - v) * c);
}

point3f image_texture_float::value(Float u, Float v, const point3f& p) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  u = fmod(u * repeatu,1);
  v = fmod(v * repeatv,1);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = data[channels*i + channels*nx*j]   * intensity;
  Float g = data[channels*i + channels*nx*j+1] * intensity;
  Float b = data[channels*i + channels*nx*j+2] * intensity;
  return(point3f(r,g,b));
}

point3f image_texture_char::value(Float u, Float v, const point3f& p) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  u = fmod(u * repeatu,1);
  v = fmod(v * repeatv,1);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = Float(data[channels*i + channels*nx*j])   * intensity * rescale;
  Float g = Float(data[channels*i + channels*nx*j+1]) * intensity * rescale;
  Float b = Float(data[channels*i + channels*nx*j+2]) * intensity * rescale;
  return(point3f(r*r,g*g,b*b));
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
  return((Float)data[channels*i + channels*nx*j+channels-1] * rescale);
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
