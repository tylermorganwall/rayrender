#ifndef TONEMAPH
#define TONEMAPH

Float reinhard(Float color, Float sum) {
  color = color*sum/(1 + sum);
  return(std::pow(color,1/2.2));
}

Float A = 0.15;
Float B = 0.50;
Float C = 0.10;
Float D = 0.20;
Float E = 0.02;
Float F = 0.30;
Float W = 11.2;

Float uncharted(Float x) {
  return(((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F);
}

Float hable(Float color) {
  Float exposure_bias = 2.0f;
  Float curr = uncharted(exposure_bias*color);
  Float whiteScale = 1.0f/uncharted(W);
  color = curr*whiteScale;
  return(std::pow(color,1/2.2));
}

Float hbd(Float color) {
  Float x = color-0.004 > 0 ? color - 0.004 : 0;
  Float retcolor = (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
  return(retcolor);
}

#endif
