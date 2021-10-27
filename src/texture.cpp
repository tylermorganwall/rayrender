#include "texture.h"

point3f image_texture::value(Float u, Float v, const point3f& p) const {
  u = fmod(u * repeatu,1);
  v = fmod(v * repeatv,1);
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = data[channels*i + channels*nx*j] * intensity;
  Float g = data[channels*i + channels*nx*j+1] * intensity;
  Float b = data[channels*i + channels*nx*j+2] * intensity;
  return(point3f(r,g,b));
}


point3f triangle_image_texture::value(Float u, Float v, const point3f& p) const {
  Float uu = ((1 - u - v) * a_u + u * b_u + v * c_u);
  Float vv = ((1 - u - v) * a_v + u * b_v + v * c_v);
  while(uu < 0) uu += 1;
  while(vv < 0) vv += 1;
  while(uu > 1) uu -= 1;
  while(vv > 1) vv -= 1;
  int i = uu * nx;
  int j = (1-vv) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = data[channels*i + channels*nx*j];
  Float g = data[channels*i + channels*nx*j+1];
  Float b = data[channels*i + channels*nx*j+2];
  return(point3f(r,g,b));
}


point3f alpha_texture::value(Float u, Float v, const point3f& p) const {
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = data[channels*i + channels*nx*j];
  return(point3f(r,r,r));
}

Float alpha_texture::channel_value(Float u, Float v, const point3f& p) const {
  Float uu = ((1 - u - v) * u_vec.x() + u * u_vec.y() + v * u_vec.z());
  Float vv = ((1 - u - v) * v_vec.x() + u * v_vec.y() + v * v_vec.z());
  while(uu < 0) uu += 1;
  while(vv < 0) vv += 1;
  while(uu > 1) uu -= 1;
  while(vv > 1) vv -= 1;
  int i = uu * nx;
  int j = (1-vv) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  return(data[channels*i + channels*nx*j + 3]);
}


point3f bump_texture::value(Float u, Float v, const point3f& p) const {
  while(u < 0) u += 1;
  while(v < 0) v += 1;
  while(u > 1) u -= 1;
  while(v > 1) v -= 1;
  int i = u * (nx-1);
  int j = (1-v) * (ny-1);
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > nx-2) i = nx-2;
  if (j > ny-2) j = ny-2;
  Float bu = (data[channels*(i+1) + channels*nx*j] - data[channels*(i-1) + channels*nx*j])/2;
  Float bv = (data[channels*i + channels*nx*(j+1)] - data[channels*i + channels*nx*(j-1)])/2;
  return(point3f(intensity*bu,intensity*bv,0));
}

point3f bump_texture::mesh_value(Float u, Float v, const point3f& p) const {
  Float uu = ((1 - u - v) * u_vec.x() + u * u_vec.y() + v * u_vec.z());
  Float vv = ((1 - u - v) * v_vec.x() + u * v_vec.y() + v * v_vec.z());
  while(uu < 0) uu += 1;
  while(vv < 0) vv += 1;
  while(uu > 1) uu -= 1;
  while(vv > 1) vv -= 1;
  int i = uu * (nx-1);
  int j = (1-vv) * (ny-1);
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > nx-2) i = nx-2;
  if (j > ny-2) j = ny-2;
  Float bu = (data[channels*(i+1) + channels*nx*j] - data[channels*(i-1) + channels*nx*j])/2;
  Float bv = (data[channels*i + channels*nx*(j+1)] - data[channels*i + channels*nx*(j-1)])/2;
  return(point3f(intensity*bu,intensity*bv,0));
}

point2f roughness_texture::value(Float u, Float v) const {
  int i = u * nx;
  int j = (1-v) * ny;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float alphax = RoughnessToAlpha(data[channels*i + channels*nx*j]);
  Float alphay = channels > 1 ? RoughnessToAlpha(data[channels*i + channels*nx*j+1]) : alphax;
  return(point2f(alphax * alphax, alphay * alphay));
}

Float roughness_texture::RoughnessToAlpha(Float roughness) {
  roughness = std::fmax(roughness, (Float)0.0001550155);
  Float x = std::log(roughness);
  return(1.62142f + 0.819955f * x + 0.1734f * x * x +
         0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
}
