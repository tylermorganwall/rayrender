
# rayballer

<img src="man/figures/rayballer.gif" ></img>

## Overview

**rayballer** is an open source R package for raytracing spheres. Based
off of Peter Shirley’s “Raytracing in One Weekend” book, this package
provides a tidy R API to the underlying raytracer to build scenes out of
diffuse (lambertian), metallic, and glass (dielectric) spheres. The
scene representation is a tibble of sphere properties, and the user can
use the provided material helper functions to easily build the scene out
of individual spheres.

## Installation

``` r
# To install the latest version from Github:
# install.packages("devtools")
devtools::install_github("tylermorganwall/rayballer")
```

## Usage

``` r
library(rayballer)

#Start with the ground
scene = generate_ground(depth=-0.5)
render_scene(scene)
```

![](man/figures/README_basicspheres-1.png)<!-- -->

``` r
#Add a sphere to the center
scene = add_sphere(scene, lambertian(x=0,y=0,z=0,radius=0.5,color=c(1,0,1)))
render_scene(scene)
```

![](man/figures/README_basicspheres-2.png)<!-- -->

``` r
#Add a metal ball (using hexcode color representation)
scene = add_sphere(scene, metal(x=0,y=0,z=1,radius=0.5,color="#ffffff",fuzz=0))
render_scene(scene)
```

![](man/figures/README_basicspheres-3.png)<!-- -->

``` r
#Add a brushed metal ball 
scene = add_sphere(scene, metal(x=0,y=1,z=0,radius=0.5,color=c(0.3,0.6,1),fuzz=0.25))
render_scene(scene)
```

![](man/figures/README_basicspheres-4.png)<!-- -->

``` r
#Add a dielectric (glass) ball
scene = add_sphere(scene, dielectric(x=0,y=0,z=-1,radius=0.5,refraction=1.6))
render_scene(scene)
```

![](man/figures/README_basicspheres-5.png)<!-- -->

``` r
#Add a grid of glass balls in front
glass_array_list = list() 

yloc = seq(-0.33,+0.33,length.out=3)
zloc = seq(-0.33,0.33,length.out=3)
locations = expand.grid(y=yloc,z=zloc)
for(i in 1:9) {
  glass_array_list[[i]] = dielectric(x=1,y=locations$y[i],z=locations$z[i], radius=0.15)
}
glass_array = do.call(rbind,glass_array_list)

scene = add_sphere(scene, glass_array)
render_scene(scene)
```

![](man/figures/README_basicspheres-6.png)<!-- -->

``` r
#Move the camera
render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15)
```

![](man/figures/README_basicspheres-7.png)<!-- -->

``` r
#Change the background gradient to a night time ambience
render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
                 backgroundhigh = "#282375", backgroundlow = "#7e77ea")
```

![](man/figures/README_basicspheres-8.png)<!-- -->

``` r
#'#Increase the aperture to give more depth of field.
render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
                 aperture = 1)
```

![](man/figures/README_basicspheres-9.png)<!-- -->

Here’s the final scene representation:

``` r
scene
```

    ##    x        y     z  radius       type             properties
    ## 1  0 -1000.50  0.00 1000.00 lambertian          0.8, 1.0, 0.0
    ## 2  0     0.00  0.00    0.50 lambertian                1, 0, 1
    ## 3  0     0.00  1.00    0.50      metal             1, 1, 1, 0
    ## 4  0     1.00  0.00    0.50      metal 0.30, 0.60, 1.00, 0.25
    ## 5  0     0.00 -1.00    0.50 dielectric                    1.6
    ## 6  1    -0.33 -0.33    0.15 dielectric                    1.5
    ## 7  1     0.00 -0.33    0.15 dielectric                    1.5
    ## 8  1     0.33 -0.33    0.15 dielectric                    1.5
    ## 9  1    -0.33  0.00    0.15 dielectric                    1.5
    ## 10 1     0.00  0.00    0.15 dielectric                    1.5
    ## 11 1     0.33  0.00    0.15 dielectric                    1.5
    ## 12 1    -0.33  0.33    0.15 dielectric                    1.5
    ## 13 1     0.00  0.33    0.15 dielectric                    1.5
    ## 14 1     0.33  0.33    0.15 dielectric                    1.5
