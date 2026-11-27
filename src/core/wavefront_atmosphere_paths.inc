R"MSL(
  // Two rays per 64-thread group. Counts stay on-device; the ordinary path and
  // shading queues retain their existing one-thread-per-ray dispatch dimensions.
  kernel void atmosphere_dispatch(device Counts &counts [[buffer(0)]],device uint *args [[buffer(1)]]) {
    args[9]=max(1u,(atomic_load_explicit(&counts.rays,memory_order_relaxed)+1)/2);
    args[10]=args[11]=1;
    args[12]=max(1u,(atomic_load_explicit(&counts.shadows,memory_order_relaxed)+1)/2);
    args[13]=args[14]=1;
  }
  float3 air_light(Light light,float3 origin,float3 direction,device const AtmosphereData &air,
                   device const float *data,device const Texture *textures,device const float4 *texels) {
    if(light.flags&128) return air_radiance(origin,direction,air,data);
    float3 value=environment(light,direction,origin,textures,texels,true);
    return any(value!=0) ? value*air_celestial(origin,direction,uint(light.color.w),air,data) : value;
  }
  // Reuse the shadow queue rather than allocating and compacting another one.
  // Visibility has already written a zero contribution for blocked connections.
  // The queue has at most one shadow per pixel per event, so lane zero can add
  // the result without atomics. Encoder order publishes transmittance before
  // this pass reads it, and queue advancement happens only after it finishes.
  kernel void shadow_lights(SCENE_ARGUMENTS) {
    uint lane = i % 32;
    i /= 32;
    if (i >= atomic_load_explicit(&counts.shadows, memory_order_relaxed)) return;
    Shadow shadow = shadow_queue[i];
    if (all(shadow.contribution.xyz == 0)) return;
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const AtmosphereData &air = *reinterpret_cast<device const AtmosphereData *>(
        bytes + storage->atmosphere_offset);
    device const float *data = reinterpret_cast<device const float *>(
        bytes + storage->atmosphere_data_offset);
    float3 value = shadow.padding[1] < p.environments
        ? air_light(lights[shadow.padding[1]], shadow.origin.xyz, shadow.direction.xyz,
                    air, data, textures, texels)
        : float3(1);
    if (lane == 0)
      pixels[shadow.pixel].radiance.xyz += shadow.contribution.xyz * value;
  }
  kernel void atmosphere_paths(SCENE_ARGUMENTS) {
    uint lane=i%32;i/=32;
    if(i>=atomic_load_explicit(&counts.rays,memory_order_relaxed)) return;
    uint pixel=ray_queue[i];
    device AtmosphereRequest *requests=reinterpret_cast<device AtmosphereRequest *>(diffusion_entries+p.camera.padding[0]);
    auto request=requests[pixel];
    if(!request.event) return;
    device AirRay *rays=reinterpret_cast<device AirRay *>(requests+p.camera.padding[0]);
    device const uchar *bytes=reinterpret_cast<device const uchar *>(storage);
    device const AtmosphereData &air=*reinterpret_cast<device const AtmosphereData *>(bytes+storage->atmosphere_offset);
    device const float *data=reinterpret_cast<device const float *>(bytes+storage->atmosphere_data_offset);
    device const Medium *media=reinterpret_cast<device const Medium *>(bytes+storage->medium_offset);
    device const float *fields=reinterpret_cast<device const float *>(bytes+storage->field_offset);
    auto ray=rays[pixel];
    uint medium=request.body ? boundaries[request.body].medium : 0;
    bool enabled=air.haze && !request.glass && (!request.body ||
      (medium && media[medium].emission.w!=0 && air.haze_in_volumes));
    if(request.event==2) {
      bool hide_background=(p.camera.padding[1]&4) && paths[pixel].depth==0;
      if(hide_background && !air.haze) return;
      air_start(ray,request.origin.xyz,request.direction.xyz,enabled);
      // Invisible boundaries and null collisions are not lighting vertices.
      // With finite haze disabled, querying at the last ray origin would move
      // the observer up to a cloud's exit and change the background spectrum.
      float3 anchor=air.haze ? (ray.direction.w ? ray.origin.xyz : request.origin.xyz)
                             : paths[pixel].previous.xyz;
      float3 sum=0;
      for(uint j=0;j<p.environments;++j) {
        if(hide_background && !(lights[j].flags&128)) continue;
        float3 value=hide_background ? air_radiance(anchor,request.direction.xyz,air,data,false)
                                    : air_light(lights[j],anchor,request.direction.xyz,air,data,textures,texels);
        value=air_remaining(ray,value,lights[j].flags&128);
        float weight=(request.padding&1) ? 1 : mis(request.weight.w,light_pdf(lights[j],request.direction.xyz,texels)/p.lights);
        sum+=weight*value;
      }
      if(lane==0) pixels[pixel].radiance.xyz+=paths[pixel].beta.xyz*sum;
    } else {
      auto segment=air_walk(ray,request.origin.xyz,request.direction.xyz,request.direction.w,
                            enabled,media[medium],air_uniform(pixel,p),air,data,fields);
      float weight=1;
      if(!(request.padding&1)) for(uint j=0;j<p.environments;++j)
        if(lights[j].flags&128) weight=mis(request.weight.w,light_pdf(lights[j],request.direction.xyz,texels)/p.lights);
      float3 previous=request.radiance.xyz;
      if(lane==0) pixels[pixel].radiance.xyz=previous+(pixels[pixel].radiance.xyz-previous)*segment.transmission+
        request.weight.xyz*segment.radiance*weight;
      if(lane==0) paths[pixel].beta.xyz*=segment.transmission;
    }
    if(lane==0) rays[pixel]=ray;
  }
)MSL"
