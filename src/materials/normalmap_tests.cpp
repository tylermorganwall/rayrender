#ifdef NOT_CRAN
#include "normalmap.h"
#include "material.h"
#include "../hitables/triangle.h"
#include "../math/loopsubdiv.h"
#include "../hitables/sphere.h"
#include "../hitables/cylinder.h"
#include "../math/animatedtransform.h"
#include <testthat.h>
#include <random>
#include <iostream>

namespace {
using normalmap::Vector;
Vector direction(double z, double phi) {
  double r = std::sqrt(std::max(0.0,1-z*z));
  return {r*std::cos(phi),r*std::sin(phi),z};
}
hit_record fixture() {
  hit_record h;
  h.p=point3f(0);h.pError=vec3f(0);h.u=h.v=0;
  h.normal=h.geometric_normal=normal3f(0,0,1);
  h.bump_normal=normal3f(.5,0,std::sqrt(.75));h.has_bump=true;
  return h;
}
}
context("Analytic diffuse normal mapping") {
 test_that("diffuse adapter support and zero-sigma continuity") {
  auto tex=std::make_shared<constant_texture>(point3f(1));
  diffuse_material lam(tex); diffuse_material rough(tex,1e-8);
  auto h=fixture();Ray r(point3f(0),-convert_to_vec3(h.bump_normal));
  for(int i=0;i<100;i++) {
    double z=-.99+1.98*i/99;
    auto d=direction(z,.2);vec3f wi(d[0],d[1],d[2]);
    auto l=lam.f(r,h,wi),o=rough.f(r,h,wi);
    expect_true(std::abs(l[0]-o[0])<2e-6);
    if(z<=0) expect_true(o[0]==0);
  }
  expect_true(rough.f(Ray(point3f(0),vec3f(1,0,0)),h,vec3f(0,0,-1))[0]==0);
 }
 test_that("EON child preserves white energy and conserves colored energy at every roughness") {
  const Vector n(0,0,1);
  double largest_error=0;
  for(double sigma:{0.,1e-12,.02,.2,.8,1.5707963267948966,100.})
  for(double cosine:{1e-8,.002,.1,.6,1.}) {
    normalmap::DiffuseChild child(sigma);
    const auto wo=direction(cosine,.3);
    double energies[3]={};
    const double colors[]={0.,.35,1.};
    const int nz=256,np=512;
    // Split at ci=co, where max(ci,co) changes branch. An unsplit uniform
    // rule misses the thin grazing interval and falsely reports energy gain.
    bool valid=true,reciprocal=true;
    for(int half=0;half<2;++half)for(int z=0;z<nz;++z)for(int a=0;a<np;++a) {
      const double lower=half==0?0:cosine, width=half==0?cosine:1-cosine;
      if(width==0)continue;
      const auto wi=direction(lower+width*(z+.5)/nz,2/normalmap::inv_pi*(a+.5)/np);
      for(int channel=0;channel<3;++channel) {
        const double rho=colors[channel];
        const double f=child.eval(wo,wi,n,rho);
        valid = valid && std::isfinite(f) && f>=0;
        reciprocal = reciprocal && std::abs(f-child.eval(wi,wo,n,rho))<1e-10*(1+f);
        energies[channel]+=f*wi[2]*width*2/normalmap::inv_pi/(nz*np);
      }
    }
    expect_true(valid);expect_true(reciprocal);
    for(int channel=0;channel<3;++channel) {
      const double error=std::abs(energies[channel]-child.directional_albedo(cosine,colors[channel]));
      largest_error=std::max(largest_error,error);
      expect_true(error<.0002);
      expect_true(energies[channel]<=colors[channel]+.0002);
    }
    expect_true(std::abs(energies[2]-1)<.0002);
    expect_true(child.eval(wo,{1,0,0},n)==0);
    expect_true(child.eval(wo,{0,0,-1},n)==0);
    expect_true(std::isfinite(child.directional_albedo(0)));
  }
  std::cout << "EON maximum furnace/quadrature error: " << largest_error << "\n";
  // At full roughness, opposing grazing tangents have a vanishing single lobe;
  // compensation must still be nonnegative. Zero roughness is exactly Lambert.
  normalmap::DiffuseChild rough(2),smooth(0);
  expect_true(rough.eval(direction(1e-12,0),direction(1e-12,acos(-1)),n)>=0);
  expect_true(smooth.eval(direction(.1,0),direction(.3,2),n,.35)==.35*normalmap::inv_pi);
 }
 test_that("diffuse adapter matches the exact published colored EON BRDF") {
  const double pi=std::acos(-1.0),rho=.35;
  const auto wi=direction(.4,.2),wo=direction(.7,2.1);
  // Independent evaluation of Eqs.10,12,13,15,19 at an interior point.
  const double r=.5,A=1/(1+(.5-2/(3*pi))*r),B=r*A;
  auto E=[&](double mu) {
    const double s=std::sqrt(1-mu*mu);
    const double G=s*(std::acos(mu)-s*mu)+(2./3)*(s/mu*(1-s*s*s)-s);
    return A+B*G/pi;
  };
  const double average=A+B*(2./3-28/(15*pi));
  const double tangent=dot(wi,wo)-wi[2]*wo[2];
  const double angular=tangent>0?tangent/std::max(wi[2],wo[2]):tangent;
  const double expected=(rho*(A+B*angular)+rho*rho*average/(1-rho*(1-average))*
                         (1-E(wi[2]))*(1-E(wo[2]))/(1-average))/pi;
  auto tex=std::make_shared<constant_texture>(point3f(rho));
  diffuse_material mapped(tex,pi/4);
  auto h=fixture();h.has_bump=false;
  Ray ray(point3f(0),vec3f(-wo[0],-wo[1],-wo[2]));
  vec3f incoming(wi[0],wi[1],wi[2]);
  expect_true(std::abs(mapped.f(ray,h,incoming)[0]-expected*wi[2])<2e-7);
  expect_true(opaque_shadow_material(&mapped)==OpaqueShadowType::Opaque);
 }
 test_that("identity, invalid tilt, smooth-normal and neutral-bump continuity") {
  Vector g(0,0,1),wo=direction(.6,.7),wi=direction(.4,2);
  for(double sigma:{0.,.1,1.,4.}) {
   normalmap::DiffuseChild child(sigma);
   double expected=child.eval(wo,wi,g);
   for(double tilt:{0.,1e-10,1e-8,1e-6,1e-4}) {
    normalmap::Model m(g,{sin(tilt),0,cos(tilt)},wo);
    expect_true(std::abs(m.eval_raw(wi,child)-expected)<.001);
   }
   expect_true(normalmap::Model(g,-g,wo).is_identity());
  }
  Vector n(.6,0,.8),u(2,0,0),v(.5,3,0);
  expect_true(normalmap::perturb(n,u,v,0,0).squared_length()>.999999);
  expect_true((normalmap::perturb(n,u,v,0,0)-n).length()<1e-12);
  for(double b:{1e-4,1e-6,1e-8}) {
   expect_true((normalmap::perturb(n,u,v,b,-b)-n).length()<2*b);
   expect_true((normalmap::perturb(n,u,v,b,-b)-normalmap::perturb(n,2.*u,3.*v,2*b,-3*b)).length()<1e-12);
   expect_true((normalmap::perturb(n,u,-v,b,b)-normalmap::perturb(n,u,v,b,-b)).length()<1e-12);
  }
  expect_true((normalmap::perturb(n,u,u,1,1)-n).length()<1e-12);
  expect_true(!normalmap::Model(g,n,wo).is_identity());
  expect_true(!normalmap::Model(Vector(0),n,wo).is_valid());
 }
 test_that("raw reciprocity and independent hemisphere energy quadrature") {
  Vector g(0,0,1);double maxenergy=0,maxerror=0;
  for(double sigma:{0.,.2,.8,1.5707963267948966})for(double tilt:{0.,.4,1.,1.55})
  for(double cosine:{.002,.1,.6,1.})for(double azimuth:{0.,1.3,3.1}) {
   Vector p(sin(tilt),0,cos(tilt)),wo=direction(cosine,azimuth);
   normalmap::Model m(g,p,wo);normalmap::DiffuseChild child(sigma);
   double energy=0;bool finite=true,support=true;
   const int nz=160,np=320;
   for(int z=0;z<nz;z++)for(int a=0;a<np;a++) {
    auto wi=direction((z+.5)/nz,2/normalmap::inv_pi*(a+.5)/np);
    double f=m.eval_raw(wi,child);
    double reverse=normalmap::Model(g,p,wi).eval_raw(wo,child);
    maxerror=std::max(maxerror,std::abs(f-reverse)/(1+f));
    energy+=f*wi[2]*2/normalmap::inv_pi/(nz*np);
    finite = finite && std::isfinite(f) && f>=0;
    support = support && (f==0 || m.pdf(wi)>0);
   }
   expect_true(finite); expect_true(support);
   maxenergy=std::max(maxenergy,energy);
   expect_true(energy<=1.003);
  }
  expect_true(maxerror<1e-10);

 }
 test_that("summed PDF plus null mass, histogram and independent estimators") {
  std::mt19937 rng(71);std::uniform_real_distribution<double> uniform(0,1);
  for(double tilt:{0.,.7,1.4})for(double azimuth:{0.,3.1}) {
   Vector g(0,0,1),p(sin(tilt),0,cos(tilt)),wo=direction(.3,azimuth);
   normalmap::Model m(g,p,wo);normalmap::DiffuseChild child(.7);
   double bins[8]={},observed[8]={},energy=0,estimate=0;
   int nz=240,np=480;
   for(int z=0;z<nz;z++)for(int a=0;a<np;a++) {
    auto wi=direction((z+.5)/nz,2/normalmap::inv_pi*(a+.5)/np);
    int bin=(wi[0]>0)+2*(wi[1]>0)+4*(wi[2]>.5);
    bins[bin]+=m.pdf(wi)*2/normalmap::inv_pi/(nz*np);
    energy+=m.eval_raw(wi,child)*wi[2]*2/normalmap::inv_pi/(nz*np);
   }
   int count=160000;double nullmass=0;
   for(int k=0;k<count;k++) {
    auto wi=m.sample(uniform(rng),uniform(rng),uniform(rng),uniform(rng));
    if(!normalmap::valid(wi)){nullmass+=1./count;continue;}
    expect_true(m.pdf(wi)>0);
    int bin=(wi[0]>0)+2*(wi[1]>0)+4*(wi[2]>.5);
    observed[bin]+=1./count;
    estimate+=m.eval_raw(wi,child)*wi[2]/m.pdf(wi)/count;
   }
   double mass=nullmass;
   for(int k=0;k<8;k++){expect_true(std::abs(bins[k]-observed[k])<.006);mass+=bins[k];}
   expect_true(std::abs(mass-1)<.006);
   expect_true(std::abs(energy-estimate)<.009);
  }
 }
 test_that("material adapter uses geometric cosine, deterministic PDF, both samplers and opaque visibility") {
  auto tex=std::make_shared<constant_texture>(point3f(2,.5,-1));
  diffuse_material mat(tex,.8);auto h=fixture();
  h.physical_shading_normal=h.bump_normal;
  random_gen rng(19),copy(19),rng2(91);RandomSampler sampler(rng2);
  Ray ray(point3f(0),vec3f(0,0,-2));vec3f wi(.3,.2,1);
  scatter_record a,b;
  expect_true(mat.scatter(ray,h,a,rng));expect_true(mat.scatter(ray,h,b,&sampler));
  expect_true(opaque_shadow_material(&mat)==OpaqueShadowType::Opaque);
  expect_true((mat.get_albedo(h)[0]==1&&mat.get_albedo(h)[2]==0));
  Float density=a.pdf_ptr->value(wi,rng);
  normalmap::Model expected({0,0,1},{.5,0,std::sqrt(.75)},{0,0,1});
  const auto input=normalmap::normalized(Vector(wi[0],wi[1],wi[2]));
  expect_true(std::abs(mat.f(ray,h,wi)[1]-expected.eval_raw(input,normalmap::DiffuseChild(.8),.5)*input[2])<2e-6);
  expect_true(mat.f(ray,h,wi)[2]==0);
  for(int i=0;i<30;i++) {
   expect_true(a.pdf_ptr->value(wi,rng)==b.pdf_ptr->value(wi,&sampler));
   expect_true(rng.unif_rand()==copy.unif_rand());
  }
  double means[2]={},nulls[2]={};
  for(int i=0;i<40000;i++) {
   bool da=false,db=false;auto va=a.pdf_ptr->generate(rng,da);auto vb=b.pdf_ptr->generate(&sampler,db);
   expect_true((da&&db));
   means[0]+=va[2]/40000.;means[1]+=vb[2]/40000.;
   nulls[0]+=(va.squared_length()==0)/40000.;nulls[1]+=(vb.squared_length()==0)/40000.;
  }
  expect_true(std::abs(means[0]-means[1])<.012);
  expect_true(std::abs(nulls[0]-nulls[1])<.012);
  expect_true(a.pdf_ptr->value(wi,rng)==density);
  auto m=normalmap::Model({0,0,1},{.5,0,sqrt(.75)},{0,0,1});
  auto v=normalmap::normalized(Vector(wi[0],wi[1],wi[2]));
  expect_true(std::abs(mat.f(ray,h,wi)[0]-m.eval_raw(v,normalmap::DiffuseChild(.8))*v[2])<2e-6);
 }
test_that("raw mesh arrays count vectors rather than scalar components") {
 Transform id;
 float vertices[] = {0,0,0, 1,0,0, 0,1,0};
 float normals[] = {0,0,1, 0,0,1, 0,0,1};
 float uv[] = {0,0, 1,0, 0,1};
 int indices[] = {0,1,2};
 auto mat = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(1)));
 for (bool attributes : {false, true}) {
   TriangleMesh mesh(vertices, indices, attributes ? normals : nullptr,
                     attributes ? uv : nullptr, 3, 3, nullptr, nullptr,
                     mat, &id, &id, false);
   expect_true(mesh.nVertices == 3);
   expect_true(mesh.nNormals == (attributes ? 3 : 0));
   expect_true(mesh.nTex == (attributes ? 3 : 0));
   expect_true(mesh.nTriangles == 1);
   expect_true(mesh.p[2][1] == 1);
   if (attributes) {
     expect_true(mesh.n[2][2] == 1);
     expect_true(mesh.uv[2][1] == 1);
   }
   mesh.ValidateMesh();
   LoopSubdivide(&mesh, 1, false);
   expect_true(mesh.nVertices == 6);
   expect_true(mesh.nTriangles == 4);
   expect_true(mesh.nNormals == 6);
   mesh.ValidateMesh();
 }
}
 test_that("triangle raw input is reciprocal across hit overloads, ray lengths, bump limits and transforms") {
  Transform id;
  float vertices[]={-2,-2,0, 2,-2,0, 0,2,0};int indices[]={0,1,2};
  float normals[]={.6,0,.8, .6,0,.8, .6,0,.8};float uv[]={0,0,1,0,.5,1};
  auto mat=std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(1)));
  TriangleMesh mesh(vertices,indices,normals,uv,3,3,nullptr,nullptr,mat,&id,&id,false);
  mesh.has_consistent_normals=true;mesh.alpha_v={.3,.3,.3};
  triangle tri(&mesh,indices,indices,indices,0,&id,&id,false);
  random_gen rng(52);RandomSampler sampler(rng);
  hit_record a,b,c,d;
  vec3f wo=unit_vector(vec3f(.2,.1,1)),wi=unit_vector(vec3f(-.3,.2,1));
  expect_true(tri.hit(Ray(point3f(wo[0],wo[1],wo[2]),-wo),0,10,a,rng));
  expect_true(tri.hit(Ray(point3f(wo[0],wo[1],wo[2]),-7*wo),0,10,b,rng));
  expect_true(tri.hit(Ray(point3f(wo[0],wo[1],wo[2]),-wo),0,10,c,&sampler));
  expect_true(tri.hit(Ray(point3f(wi[0],wi[1],wi[2]),-wi),0,10,d,&sampler));
  for(int k=0;k<3;k++) {
   expect_true(std::abs(a.normal[k]-b.normal[k])<2e-6);
   expect_true(std::abs(a.normal[k]-c.normal[k])<2e-6);
   expect_true(std::abs(a.physical_shading_normal[k]-d.physical_shading_normal[k])<2e-6);
   expect_true(std::abs(a.p[k]-b.p[k])<2e-6);
  }
  double forward=mat->f(Ray(point3f(0),-wo),a,wi)[0]/wi[2];
  double reverse=mat->f(Ray(point3f(0),-wi),d,wo)[0]/wo[2];
  expect_true(std::abs(forward-reverse)<2e-6);
  // Texture storage is owned here, not by bump_texture (same cache lifetime as production).
  // Zero displacement is neutral; a nonzero constant height also changes
  // curved shading tangents through PBRT's height * normal-derivative term.
  unsigned char pixels[16];std::fill(pixels,pixels+16,0);
  auto bump=std::make_shared<bump_texture>(pixels,4,4,1,1,1,1);
  mesh.bump_textures[0]=bump;
  expect_true(tri.hit(Ray(point3f(0,0,1),vec3f(0,0,-1)),0,10,b,rng));
  expect_true(tri.hit(Ray(point3f(0,0,1),vec3f(0,0,-1)),0,10,c,&sampler));
  for(int k=0;k<3;k++) {
   expect_true(std::abs(a.physical_shading_normal[k]-b.physical_shading_normal[k])<2e-6);
   expect_true(std::abs(b.physical_shading_normal[k]-c.physical_shading_normal[k])<2e-6);
  }
  for(int y=0;y<4;++y)for(int x=0;x<4;++x)pixels[y*4+x]=30*x+10*y;
  expect_true(tri.hit(Ray(point3f(0,0,1),vec3f(0,0,-1)),0,10,b,rng));
  expect_true(tri.hit(Ray(point3f(0,0,1),vec3f(0,0,-1)),0,10,c,&sampler));
  expect_true((b.physical_shading_normal-a.physical_shading_normal).length()>1e-4);
  for(int k=0;k<3;++k) {
    expect_true(std::abs(b.physical_shading_normal[k]-c.physical_shading_normal[k])<2e-6);
    expect_true(b.bump_normal[k]==b.physical_shading_normal[k]);
  }
  // Vertex colors use barycentric texture coordinates in both overloads.
  mesh.has_vertex_colors=true;
  expect_true(tri.hit(Ray(point3f(.2,.1,1),vec3f(0,0,-1)),0,10,b,rng));
  expect_true(tri.hit(Ray(point3f(.2,.1,1),vec3f(0,0,-1)),0,10,c,&sampler));
  expect_true(std::abs(b.u-c.u)<2e-6);expect_true(std::abs(b.v-c.v)<2e-6);
  for(auto transform:{Scale(2,3,.5),Scale(-2,3,.5),RotateY(35)*Scale(2,1,3)}) {
   const hit_record constant=b;
   auto mapped=transform(b),other=transform(constant);
   expect_true(std::abs(mapped.physical_shading_normal.length()-1)<2e-6);
   expect_true(std::abs(mapped.geometric_normal.length()-1)<2e-6);
   for(int k=0;k<3;k++)expect_true(mapped.physical_shading_normal[k]==other.physical_shading_normal[k]);
  }
 }
 test_that("primitive bump frames match both overloads and survive reflected and animated transforms") {
  Transform identity,scaled=Scale(-1.5,2,.7),inverse=Inverse(scaled);
  auto mat=std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(1)));
  unsigned char pixels[64];for(int y=0;y<8;y++)for(int x=0;x<8;x++)pixels[y*8+x]=20*x+10*y;
  auto bump=std::make_shared<bump_texture>(pixels,8,8,1,2);
  sphere ball(1,mat,nullptr,bump,&scaled,&inverse,false);
  cylinder tube(1,2,0,2*M_PI,true,mat,nullptr,bump,&scaled,&inverse,false);
  random_gen rng(199);RandomSampler sampler(rng);
  for(const hitable* object:{static_cast<hitable*>(&ball),static_cast<hitable*>(&tube)}) {
   for(point3f origin:{point3f(2,.2,2),point3f(0,.2,0)}) {
    vec3f d=origin[0]>0?vec3f(-1,0,-1):vec3f(1,0,1);
    Ray ray(scaled(origin),scaled(d));hit_record a,b;
    expect_true(object->hit(ray,0,10,a,rng));
    expect_true(object->hit(ray,0,10,b,&sampler));
    for(int k=0;k<3;k++) {
     expect_true(std::abs(a.geometric_normal[k]-b.geometric_normal[k])<2e-6);
     expect_true(std::abs(a.physical_shading_normal[k]-b.physical_shading_normal[k])<2e-6);
    }
    expect_true(std::abs(a.physical_shading_normal.length()-1)<2e-6);
    expect_true(std::abs(a.geometric_normal.length()-1)<2e-6);
   }
  }
  auto h=fixture();h.physical_shading_normal=h.bump_normal;
  Transform end=Translate(vec3f(1,2,3))*RotateY(35)*Scale(2,3,.5);
  AnimatedTransform animation(&identity,0,&end,1);
  for(Float time:{0.f,.3f,1.f}) {
   Transform transform;animation.Interpolate(time,&transform);
   auto mapped=transform(h);
   expect_true(std::abs(mapped.geometric_normal.length()-1)<2e-6);
   expect_true(std::abs(mapped.physical_shading_normal.length()-1)<2e-6);
  }
 }
 test_that("finite-area light and BSDF proposals and their MIS mixture match area quadrature") {
  normalmap::Model model({0,0,1},{.8,0,.6},direction(.5,2.3));
  normalmap::DiffuseChild child(1.5707963267948966);
  // A unit-radiance 2x2 emitter at z=2, x in [1,3], y in [-1,1].
  // Integrate in area measure independently of either directional proposal.
  double target=0;
  const int resolution=300;
  for(int i=0;i<resolution;i++) for(int j=0;j<resolution;j++) {
    Vector point(1+2*(i+.5)/resolution,-1+2*(j+.5)/resolution,2);
    auto wi=normalmap::normalized(point);
    target+=model.eval_raw(wi,child)*wi[2]*wi[2]/point.squared_length()*4/(resolution*resolution);
  }
  auto light_pdf=[](Vector wi) {
    if(wi[2]<=0) return 0.;
    Vector p=2/wi[2]*wi;
    return p[0]>=1 && p[0]<=3 && p[1]>=-1 && p[1]<=1 ? p.squared_length()/(4*wi[2]) : 0.;
  };
  std::mt19937 rng(2047);std::uniform_real_distribution<double> u(0,1);
  const int count=200000;
  for(int proposal=0;proposal<3;proposal++) {
    double sum=0,squares=0;
    for(int k=0;k<count;k++) {
      bool light=proposal==1 || (proposal==2 && u(rng)<.5);
      Vector wi=light ? normalmap::normalized(Vector(1+2*u(rng),-1+2*u(rng),2)) :
                        model.sample(u(rng),u(rng),u(rng),u(rng));
      double lpdf=light_pdf(wi),bpdf=model.pdf(wi);
      double density=proposal==0?bpdf:proposal==1?lpdf:.5*(lpdf+bpdf);
      double estimate=density>0 && lpdf>0 ? model.eval_raw(wi,child)*wi[2]/density : 0;
      sum+=estimate; squares+=estimate*estimate;
    }
    double mean=sum/count;
    double error=std::sqrt(std::max(0.,squares/count-mean*mean)/count);
    expect_true(std::abs(mean-target)<5*error+2e-5);
  }
 }

}
#endif
