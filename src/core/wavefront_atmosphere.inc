R"MSL(
  // Each atmosphere path/shadow is owned by one 32-lane SIMD group. All lanes
  // traverse identical geometry and make identical stochastic decisions; only
  // spectral work is partitioned. Reductions must remain in convergent control
  // flow, and the calling kernels commit shared state only from lane zero.
  struct AtmosphereData {
    skymodelr_gpu::Layout model;
    float4 transform[3], origin, gain, settings, options;
    float4 weights[16], transmission[16], solar[16], lunar[16], filter[129];
    uint channels, sky, sun, haze, query_altitude, haze_filter, haze_in_volumes, padding;
  };
  struct AirRay { float4 origin, direction, transmission, radiance; };
  static_assert(sizeof(AtmosphereData)==3376 && sizeof(AirRay)==64, "Prague atmosphere ABI mismatch");
  struct AirSegment { float3 transmission, radiance; };
  float3 air_direction(float3 w,device const AtmosphereData &a) {
    return normalize(float3(dot(a.transform[0].xyz,w),dot(a.transform[1].xyz,w),dot(a.transform[2].xyz,w)));
  }
  float3 air_position(float3 p,device const AtmosphereData &a,bool lighting=false) {
    if(lighting && !a.query_altitude) return float3(0,0,a.origin.w);
    float3 d=p-a.origin.xyz;
    return float3(dot(a.transform[0].xyz,d),dot(a.transform[1].xyz,d),dot(a.transform[2].xyz,d))*a.gain.w+float3(0,0,a.origin.w);
  }
  float3 air_clamp(float3 p) {
    float h=skymodelr_gpu::height(p);
    if(h<0 || h>15000) p+=normalize(p+float3(0,0,6378000))*(clamp(h,0.f,15000.f)-h);
    return p;
  }
  skymodelr_gpu::Parameters air_parameters(float3 p,float3 w,device const AtmosphereData &a) {
    return skymodelr_gpu::parameters(air_clamp(p),w,a.settings.x,a.settings.y,a.settings.z,a.settings.w);
  }
  bool air_occluded(float3 p,float3 w,float distance) {
    p=air_clamp(p);
    float b=dot(p,w)+6378000*w.z;
    if(b>=0) return false;
    float c=dot(p,p)+2*6378000*p.z, disc=b*b-c;
    if(disc<=0) return false;
    float t=c/(-b+sqrt(disc));
    return t<distance && t>=-1e-5f;
  }
  float normalized_air_transmission(float zero,float value,float distance) {
    value=zero>0 ? clamp(value/zero,0.f,1.f) : 0;
    return distance<100 ? pow(value,distance/100) : value;
  }
  struct AirSpectrumQuery { skymodelr_gpu::SkyQuery sky[3]; float3 weights; };
  AirSpectrumQuery air_prepare_sky(float3 position,float3 direction,bool smooth,
                                    device const AtmosphereData &a,device const float *data) {
    AirSpectrumQuery query={};
    auto p=air_parameters(position,direction,a);
    float near=smooth ? clamp((5*pi/180-abs(p.theta-pi/2))/(pi/180),0.f,1.f) : 0;
    float3 up=normalize(position+float3(0,0,6378000));
    float3 vertical=normalize(up-direction*dot(up,direction));
    for(int j=-1;j<=1;++j) {
      query.weights[j+1]=j==0 ? 1-.5f*near : .25f*near;
      if(query.weights[j+1]==0) continue;
      float angle=j*.2f*pi/180;
      auto sample=j==0 ? p : air_parameters(position,normalize(direction*cos(angle)+vertical*sin(angle)),a);
      query.sky[j+1]=skymodelr_gpu::prepare_sky(sample,a.model,data);
    }
    return query;
  }
  float air_prepared_sky(thread const AirSpectrumQuery &query,float wavelength,
                         device const AtmosphereData &a,device const float *data) {
    float result=0;
    for(uint j=0;j<3;++j) if(query.weights[j]>0)
      result+=query.weights[j]*skymodelr_gpu::sky_prepared(query.sky[j],wavelength,a.model,data);
    return result;
  }
  float3 air_radiance(float3 origin,float3 w,device const AtmosphereData &a,device const float *data,
                      bool include_sun = true) {
    float3 position=air_position(origin,a,true),direction=air_direction(w,a),rgb=0;
    auto p=air_parameters(position,direction,a);
    bool sun=include_sun && a.sun && p.gamma<=a.options.x && !air_occluded(position,direction,INFINITY);
    AirSpectrumQuery sky={};if(a.sky) sky=air_prepare_sky(position,direction,true,a,data);
    skymodelr_gpu::TransQuery transmission={};
    if(sun) transmission=skymodelr_gpu::prepare_transmission(p,INFINITY,a.model,data);
    for(uint c=simd_prefix_exclusive_sum(1u);c<a.channels;c+=32) {
      float wavelength=a.weights[c].w;
      float value=a.sky ? air_prepared_sky(sky,wavelength,a,data) : 0;
      if(sun) value+=skymodelr_gpu::linear(data,a.model.sun.offset+uint(wavelength-280),fract(wavelength-280))*
        skymodelr_gpu::transmission_prepared(transmission,wavelength,a.model,data);
      rgb+=a.weights[c].xyz*value;
    }
    return simd_sum(rgb)*a.gain.xyz;
  }
  // Camera coverage uses normalized, nonnegative RGB extinction weights, not
  // the solar/lunar spectra used to illuminate the scene. Infinity needs its
  // own query: constructing an endpoint as origin + infinity * direction can
  // create NaNs, and finite-segment in-scattering is irrelevant to coverage.
  float3 air_infinite_transmission(float3 origin,float3 w,device const AtmosphereData &a,
                                    device const float *data) {
    float3 position=air_position(origin,a,true),direction=air_direction(w,a),rgb=0;
    if(air_occluded(position,direction,INFINITY)) return 0;
    auto p=air_parameters(position,direction,a);
    auto zero=skymodelr_gpu::prepare_transmission(p,1e-6f,a.model,data);
    auto far=skymodelr_gpu::prepare_transmission(p,INFINITY,a.model,data);
    for(uint c=simd_prefix_exclusive_sum(1u);c<a.channels;c+=32) {
      float wavelength=a.weights[c].w;
      rgb+=a.transmission[c].xyz*normalized_air_transmission(
        skymodelr_gpu::transmission_prepared(zero,wavelength,a.model,data),
        skymodelr_gpu::transmission_prepared(far,wavelength,a.model,data),INFINITY);
    }
    return clamp(simd_sum(rgb),0.f,1.f);
  }
  float3 air_celestial(float3 origin,float3 w,uint spectrum,device const AtmosphereData &a,device const float *data) {
    float3 position=air_position(origin,a,true),direction=air_direction(w,a),rgb=0;
    if(air_occluded(position,direction,INFINITY)) return 0;
    auto p=air_parameters(position,direction,a);
    auto query=skymodelr_gpu::prepare_transmission(p,INFINITY,a.model,data);
    for(uint c=simd_prefix_exclusive_sum(1u);c<a.channels;c+=32) {
      float tr=skymodelr_gpu::transmission_prepared(query,a.weights[c].w,a.model,data);
      rgb+=(spectrum==1 ? a.solar[c].xyz : a.lunar[c].xyz)*tr;
    }
    return simd_sum(rgb);
  }
  AirSegment air_segment(float3 origin,float3 w,float distance,float uniform,
                         device const AtmosphereData &a,device const float *data) {
    AirSegment result={float3(1),float3(0)};
    if(!a.haze || distance<=0) return result;
    float d=distance*a.gain.w;
    float3 position=air_position(origin,a),direction=air_direction(w,a);
    auto p=air_parameters(position,direction,a);
    float blend=a.haze_filter ? clamp((p.gamma*180/pi-3)/3,0.f,1.f) : 0;
    blend=blend*blend*(3-2*blend);
    float3 filtered=0,source=0,tr=0,filter_direction=direction;
    bool valid=false;
    // Same 129-node finite-path quadrature as the CPU; choose one direction,
    // renormalizing after rejecting paths through Earth, before wavelength work.
    if(blend>0 && a.sky) {
      float3 up=normalize(position+float3(0,0,6378000));
      float3 vertical=up-direction*dot(up,direction);
      if(length(vertical)<1e-10f) vertical=cross(direction,abs(direction.x)<.9f ? float3(1,0,0) : float3(0,1,0));
      vertical=normalize(vertical);
      float mass=0;
      for(uint j=0;j<129;++j) {
        float3 ray=normalize(direction*a.filter[j].x+vertical*a.filter[j].y);
        if(!air_occluded(position,ray,d)) mass+=a.filter[j].z;
      }
      float choice=uniform*mass;
      for(uint j=0;j<129;++j) {
        float3 ray=normalize(direction*a.filter[j].x+vertical*a.filter[j].y);
        if(air_occluded(position,ray,d)) continue;
        filter_direction=ray;valid=true;choice-=a.filter[j].z;
        if(choice<=0) break;
      }
    }
    auto fp=air_parameters(position,filter_direction,a);
    auto zero=skymodelr_gpu::prepare_transmission(p,1e-6f,a.model,data);
    auto far=skymodelr_gpu::prepare_transmission(p,max(100.f,d),a.model,data);
    skymodelr_gpu::SkyQuery filter_start={},filter_end={};
    skymodelr_gpu::TransQuery filter_zero={},filter_far={};
    AirSpectrumQuery direct_start={},direct_end={};
    if(a.sky && blend<1) {
      direct_start=air_prepare_sky(position,direction,true,a,data);
      direct_end=air_prepare_sky(position+direction*d,direction,true,a,data);
    }
    if(a.sky && blend>0 && valid) {
      filter_start=skymodelr_gpu::prepare_sky(fp,a.model,data);
      filter_end=skymodelr_gpu::prepare_sky(air_parameters(position+filter_direction*d,filter_direction,a),a.model,data);
      filter_zero=skymodelr_gpu::prepare_transmission(fp,1e-6f,a.model,data);
      filter_far=skymodelr_gpu::prepare_transmission(fp,max(100.f,d),a.model,data);
    }
    for(uint c=simd_prefix_exclusive_sum(1u);c<a.channels;c+=32) {
      float lambda=a.weights[c].w;
      float transmission=normalized_air_transmission(skymodelr_gpu::transmission_prepared(zero,lambda,a.model,data),
        skymodelr_gpu::transmission_prepared(far,lambda,a.model,data),d);
      tr+=a.transmission[c].xyz*transmission;
      if(!a.sky) continue;
      if(blend<1) {
        float start=air_prepared_sky(direct_start,lambda,a,data);
        float end=air_prepared_sky(direct_end,lambda,a,data);
        source+=a.weights[c].xyz*max(0.f,start-transmission*end);
      }
      if(blend>0 && valid) {
        float start=skymodelr_gpu::sky_prepared(filter_start,lambda,a.model,data);
        float end=skymodelr_gpu::sky_prepared(filter_end,lambda,a.model,data);
        float t=normalized_air_transmission(skymodelr_gpu::transmission_prepared(filter_zero,lambda,a.model,data),
          skymodelr_gpu::transmission_prepared(filter_far,lambda,a.model,data),d);
        filtered+=a.weights[c].xyz*max(0.f,start-t*end);
      }
    }
    result.transmission=clamp(simd_sum(tr),0.f,1.f);
    result.radiance=mix(simd_sum(source),simd_sum(filtered),blend)*a.gain.xyz;
    return result;
  }
  void air_start(thread AirRay &ray,float3 p,float3 w,bool enabled) {
    if(ray.origin.w==0 || bool(ray.direction.w)!=enabled || any(ray.direction.xyz!=w)) {
      ray.origin=float4(p,1);ray.direction=float4(w,enabled);
      ray.transmission=float4(1,1,1,0);ray.radiance=0;
    }
  }
  AirSegment air_advance(thread AirRay &ray,float3 endpoint,float uniform,
                         device const AtmosphereData &a,device const float *data) {
    AirSegment result={float3(1),float3(0)};
    if(!ray.direction.w) return result;
    float next=dot(endpoint-ray.origin.xyz,ray.direction.xyz);
    if(next<=ray.transmission.w) return result;
    auto value=air_segment(ray.origin.xyz,ray.direction.xyz,next,uniform,a,data);
    for(uint c=0;c<3;++c) {
      result.transmission[c]=ray.transmission[c]>0 ? value.transmission[c]/ray.transmission[c] : 0;
      result.radiance[c]=ray.transmission[c]>0 ? (value.radiance[c]-ray.radiance[c])/ray.transmission[c] : 0;
    }
    ray.transmission=float4(value.transmission,next);ray.radiance=float4(value.radiance,0);
    return result;
  }
  float3 air_remaining(AirRay ray,float3 value,bool sky) {
    if(!ray.direction.w) return value;
    return select(float3(0),(value-(sky ? ray.radiance.xyz : float3(0)))/ray.transmission.xyz,ray.transmission.xyz>0);
  }
  uint air_mix(uint x) { x^=x>>16;x*=0x7feb352du;x^=x>>15;x*=0x846ca68bu;return x^(x>>16); }
  float air_uniform(uint pixel,constant Params &p) {
    return min((float(air_mix(pixel^air_mix(p.sample)^air_mix(p.depth)^p.seed))+.5f)*0x1p-32f,0x1.fffffep-1f);
  }
)MSL"
