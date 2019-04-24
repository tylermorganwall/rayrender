#ifndef TONEMAPH
#define TONEMAPH

float reinhard(float color, float sum) {
  color = color*sum/(1 + sum);
  return(pow(color,1/2.2));
}

float A = 0.15;
float B = 0.50;
float C = 0.10;
float D = 0.20;
float E = 0.02;
float F = 0.30;
float W = 11.2;

float uncharted(float x) {
  return(((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F);
}

float hable(float color) {
  float exposure_bias = 2.0f;
  float curr = uncharted(exposure_bias*color);
  float whiteScale = 1.0f/uncharted(W);
  color = curr*whiteScale;
  return(pow(color,1/2.2));
}

float hbd(float color) {
  float x = color-0.004 > 0 ? color - 0.004 : 0;
  float retcolor = (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
  return(retcolor);
}

#endif
