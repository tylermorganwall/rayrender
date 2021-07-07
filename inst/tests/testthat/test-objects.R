library(testthat)

test_samples = 16

#Set seed per render
render_scene = function(...) {
  set.seed(1)
  rayrender::render_scene(...)
}

image_sums = list()
counter = 1

previous_sums = as.numeric(readLines(system.file("testdata","test-object_image_sums.txt",package="rayrender")))

#Generate a sphere in the cornell box.
generate_cornell() %>%
  add_object(sphere(x = 555/2, y = 555/2, z = 555/2, radius = 100)) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Basic sphere in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})
counter = counter + 1


#Generate a gold sphere in the cornell box
generate_cornell() %>%
  add_object(sphere(x = 555/2, y = 100, z = 555/2, radius = 100, 
                    material = microfacet(color = "gold"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Metal sphere in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Add motion blur and show the sphere moving
generate_cornell() %>%
  add_object(sphere(x = 555/2, y = 100, z = 555/2, radius = 100,
                    material = microfacet(color = "gold"), velocity = c(50, 0, 0))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Moving sphere in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a cube in the cornell box.
generate_cornell() %>%
  add_object(cube(x = 555/2, y = 100, z = 555/2, 
                  xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Basic cube in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Generate a gold cube in the cornell box
generate_cornell() %>%
  add_object(cube(x = 555/2, y = 100, z = 555/2, 
                  xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
                  material = metal(color = "gold", fuzz = 0.2))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Metal cube in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a rotated dielectric box in the cornell box
generate_cornell() %>%
  add_object(cube(x = 555/2, y = 200, z = 555/2, 
                  xwidth = 200, ywidth = 100, zwidth = 200, angle = c(30, 30, 30),
                  material = dielectric())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40,  
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Glass cube in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a purple rectangle in the cornell box.
generate_cornell() %>%
  add_object(xy_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, ywidth = 200,
                     material = diffuse(color = "purple"))) %>%
  render_scene(lookfrom = c(278, 278, -800), lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Basic xyrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a gold plane in the cornell box

generate_cornell() %>%
  add_object(xy_rect(x = 555/2, y = 100, z = 555/2, 
                     xwidth = 200, ywidth = 200, angle = c(0, 30, 0),
                     material = metal(color = "gold"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Metal xyrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a purple rectangle in the cornell box.

generate_cornell() %>%
  add_object(xz_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, zwidth = 200,
                     material = diffuse(color = "purple"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Basic xzrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a gold plane in the cornell box

generate_cornell() %>%
  add_object(xz_rect(x = 555/2, y = 100, z = 555/2, 
                     xwidth = 200, zwidth = 200, angle = c(0, 30, 0),
                     material = metal(color = "gold"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Metal xzrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a purple rectangle in the cornell box.

generate_cornell() %>%
  add_object(yz_rect(x = 100, y = 100, z = 555/2, ywidth = 200, zwidth = 200,
                     material = diffuse(color = "purple"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Basic yzrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Generate a gold plane in the cornell box

generate_cornell() %>%
  add_object(yz_rect(x = 100, y = 100, z = 555/2, 
                     ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
                     material = metal(color = "gold"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5) %>% sum() ->
  image_sums[[counter]]
test_that("Metal yzrect in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate and render the default Cornell box.
scene = generate_cornell()

render_scene(scene, samples=test_samples,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE) %>% sum() ->
  image_sums[[counter]]
test_that("Cornell box render", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` .

#Make a much smaller light in the center of the room.
scene = generate_cornell(lightwidth=200,lightdepth=200)

render_scene(scene, samples=test_samples,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE) %>% sum() ->
  image_sums[[counter]]
test_that("Cornell box render, small light", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` .

#Place a sphere in the middle of the box.
scene = scene %>%
  add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/4))

render_scene(scene, samples=test_samples,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)  %>% sum() ->
  image_sums[[counter]]
test_that("Cornell box render, sphere", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` .

#Reduce "fireflies" by setting a clamp_value in render_scene()

render_scene(scene, samples=test_samples,aperture=0, fov=40, ambient_light=FALSE, 
             parallel=TRUE,clamp_value=3)  %>% sum() ->
  image_sums[[counter]]
test_that("Cornell box render, clamped output", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` .
# Change the color scheme of the cornell box

new_cornell = generate_cornell(leftcolor="purple", rightcolor="yellow")
render_scene(new_cornell, samples=test_samples,aperture=0, fov=40, ambient_light=FALSE, 
             parallel=TRUE,clamp_value=3)  %>% sum() ->
  image_sums[[counter]]
test_that("Cornell box render, cystom color scheme", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` .

#Generate the ground and add some objects
scene = generate_ground(depth=-0.5,
                        material = diffuse(noise=1,noisecolor="blue",noisephase=10)) %>%
  add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
  add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))

render_scene(scene, samples=test_samples, parallel=TRUE,lookfrom=c(0,2,10))  %>% sum() ->
  image_sums[[counter]]
test_that("Generate ground render", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


# Make the sphere representing the ground larger and make it a checkered surface.
scene = generate_ground(depth=-0.5, spheresize=10000,
                        material = diffuse(checkercolor="grey50")) %>%
  add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
  add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))

render_scene(scene, samples=test_samples,parallel=TRUE,lookfrom=c(0,1,10))  %>% sum() ->
  image_sums[[counter]]
test_that("Generate ground render, larger sphere + checkered", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate the ground and add some objects
scene = generate_studio(depth=-1, material = diffuse(color="white")) %>%
  add_object(obj_model(r_obj(),y=-1,x=0.7,material=glossy(color="darkred"),angle=c(0,-20,0))) %>%
  add_object(sphere(x=-0.7,radius=0.5,material=dielectric())) %>% 
  add_object(sphere(y=3,x=-2,z=20,material=light(intensity=600)))

render_scene(scene, parallel=TRUE,lookfrom=c(0,2,10),fov=20,clamp_value=10,samples=test_samples)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate studio render with objects", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Zooming out to show the full default scene

render_scene(scene, parallel=TRUE,lookfrom=c(0,200,400),clamp_value=10,samples=test_samples)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate studio render, zoomed out", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a segment in the cornell box. 

generate_cornell() %>%
  add_object(segment(start = c(100, 100, 100), end = c(455, 455, 455), radius = 50)) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate segment in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Error in render_scene_rcpp(camera_info = camera_info, scene_info = scene_info): Index out of bounds: [index=6; extent=6].

# Draw a line graph representing a normal distribution, but with metal:
xvals = seq(-3, 3, length.out = 30)
yvals = dnorm(xvals)

scene_list = list()
for(i in 1:(length(xvals) - 1)) {
  scene_list[[i]] = segment(start = c(555/2 + xvals[i] * 80, yvals[i] * 800, 555/2),
                            end = c(555/2 + xvals[i + 1] * 80, yvals[i + 1] * 800, 555/2),
                            radius = 10,
                            material = metal())
}
scene_segments = do.call(rbind,scene_list)

generate_cornell() %>% 
  add_object(scene_segments) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate metal segments in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Draw the outline of a cube:

cube_outline = segment(start = c(100, 100, 100), end = c(100, 100, 455), radius = 10) %>%
  add_object(segment(start = c(100, 100, 100), end = c(100, 455, 100), radius = 10)) %>%
  add_object(segment(start = c(100, 100, 100), end = c(455, 100, 100), radius = 10)) %>%
  add_object(segment(start = c(100, 100, 455), end = c(100, 455, 455), radius = 10)) %>%
  add_object(segment(start = c(100, 100, 455), end = c(455, 100, 455), radius = 10)) %>%
  add_object(segment(start = c(100, 455, 455), end = c(100, 455, 100), radius = 10)) %>%
  add_object(segment(start = c(100, 455, 455), end = c(455, 455, 455), radius = 10)) %>%
  add_object(segment(start = c(455, 455, 100), end = c(455, 100, 100), radius = 10)) %>%
  add_object(segment(start = c(455, 455, 100), end = c(455, 455, 455), radius = 10)) %>%
  add_object(segment(start = c(455, 100, 100), end = c(455, 100, 455), radius = 10)) %>%
  add_object(segment(start = c(455, 100, 455), end = c(455, 455, 455), radius = 10)) %>%
  add_object(segment(start = c(100, 455, 100), end = c(455, 455, 100), radius = 10))


generate_cornell() %>%
  add_object(cube_outline) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate cube made out of segments in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Shrink and rotate the cube

generate_cornell() %>%
  add_object(group_objects(cube_outline, pivot_point = c(555/2, 555/2, 555/2),
                           group_angle = c(45,45,45), group_scale = c(0.5,0.5,0.5))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate scaled/rotated cube made out of segments in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Error in render_scene_rcpp(camera_info = camera_info, scene_info = scene_info): Index out of bounds: [index=6; extent=6].

#Generate a cylinder in the cornell box. Add a cap to both ends.


generate_cornell() %>%
  add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
                      length = 300, radius = 100, material = metal())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate cylinder in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Rotate the cylinder

generate_cornell() %>%
  add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
                      length = 300, radius = 100, angle = c(0, 0, 45),
                      material = diffuse())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate rotated cylinder in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


# Only render a subtended arc of the cylinder,

generate_cornell(lightintensity=3) %>%
  add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
                      length = 300, radius = 100, angle = c(45, 0, 0), phi_min = 0, phi_max = 180,
                      material = diffuse())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate partial cylinder in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a triangle in the Cornell box.

generate_cornell() %>%
  add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
                      material = diffuse(color = "purple"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate triangle in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Pass individual colors to each vertex: 

generate_cornell() %>%
  add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
                      color1 = "green", color2 = "yellow", color3 = "red")) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate rainbow triangle in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate an ellipsoid in a Cornell box

generate_cornell() %>%
  add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
                       a = 100, b = 50, c = 50)) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate ellipsoid in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Change the axes to make it taller rather than wide:

generate_cornell() %>%
  add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
                       a = 100, b = 200, c = 100, material = metal())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate tall ellipsoid in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Rotate it and make it dielectric:

generate_cornell() %>%
  add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
                       a = 100, b = 200, c = 100, angle = c(0, 0, 45),
                       material = dielectric())) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate rotated glass ellipsoid in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a disk in the cornell box.

generate_cornell() %>%
  add_object(disk(x = 555/2, y = 50, z = 555/2, radius = 150, 
                  material = diffuse(color = "orange"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate disk in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Rotate the disk.

generate_cornell() %>%
  add_object(disk(x = 555/2, y = 555/2, z = 555/2, radius = 150, angle = c(45, 0, 0), 
                  material = diffuse(color = "orange"))) %>%
  render_scene(lookfrom = c(278, 278, -800) , lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate rotated disk in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Pass a value for the inner radius.

generate_cornell() %>% 
  add_object(disk(x = 555/2, y = 555/2, z = 555/2, 
                  radius = 150, inner_radius = 75, angle = c(45, 0, 0), 
                  material = diffuse(color = "orange"))) %>%
  render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
               ambient_light = FALSE, samples = test_samples, parallel = TRUE, clamp_value = 5)  %>% sum() ->
  image_sums[[counter]]
test_that("Generate disk with hole in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Load the included example R object file, by calling the r_obj() function. This
#returns the local file path to the `r.txt` obj file. The file extension is "txt" 
#due to package constraints, but the file contents are identical and it does not 
#affect the function.


generate_ground(material = diffuse(checkercolor = "grey50")) %>%
  add_object(obj_model(y = -0.8, filename = r_obj(),
                       material = metal(color = "gold", fuzz = 0.025))) %>%
  add_object(obj_model(x = 1.8, y = -0.8, filename = r_obj(), 
                       material = diffuse(color = "lightblue"))) %>%
  add_object(obj_model(x = -1.8, y = -0.8, filename = r_obj() , 
                       material = dielectric(color = "pink"))) %>%
  add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
                    material = light(intensity = 20))) %>%
  render_scene(parallel = TRUE, samples = test_samples, 
               tonemap = "reinhold", aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10))  %>% sum() ->
  image_sums[[counter]]
test_that("Render OBJ file in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Use scale_obj to make objects bigger--this is more robust than the generic scale argument.

generate_ground(material = diffuse(checkercolor = "grey50")) %>%
  add_object(obj_model(y = -0.8, filename = r_obj(), scale_obj = 2,
                       material = diffuse(noise = TRUE, noiseintensity = 10,noisephase=45))) %>%
  add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
                    material = light(intensity = 10))) %>%
  render_scene(parallel = TRUE, samples = test_samples, ambient = TRUE, 
               backgroundhigh="blue", backgroundlow="red",
               aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10),
               lookat = c(0,1,0))   %>% sum() ->
  image_sums[[counter]]
test_that("Render scaled OBJ file in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Load a mesh3d object (from the Rvcg) and render it:
suppressWarnings(library(Rvcg))
data(humface)

generate_studio() %>% 
  add_object(mesh3d_model(humface,y=-0.3,x=0,z=0,
                          material=glossy(color="dodgerblue4"), scale_mesh = 1/70)) %>%
  add_object(sphere(y=5,x=5,z=5,material=light(intensity=50))) %>% 
  render_scene(samples=test_samples,width=400,height=400,
               lookat = c(0,0.5,1), aperture=0.0)  %>% sum() ->
  image_sums[[counter]]
test_that("Render mesh3d file in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a cone in a studio, pointing upwards:

generate_studio() %>% 
  add_object(cone(start=c(0,-1,0), end=c(0,1,0), radius=1,material=diffuse(color="red"))) %>% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
  render_scene(samples=test_samples,clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render cone in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Change the radius, length, and direction
generate_studio() %>% 
  add_object(cone(start=c(0,0,0), end=c(0,-1,0), radius=0.5,material=diffuse(color="red"))) %>% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
  render_scene(samples=test_samples,clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render custom cone in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Give custom start and end points (and customize the color/texture)
generate_studio() %>% 
  add_object(cone(start=c(-1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="red"))) %>%
  add_object(cone(start=c(1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="green"))) %>%
  add_object(cone(start=c(0,1,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="orange"))) %>% 
  add_object(cone(start=c(-1,-0.5,0), end=c(1,-0.5,0), radius=0.25,
                  material = diffuse(color="red",gradient_color="green"))) %>% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
  render_scene(samples=test_samples,clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render start/end cone in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Specify cone via direction and location, instead of start and end positions
#Length is derived from the magnitude of the direction.
gold_mat = microfacet(roughness=0.1,eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661))
generate_studio() %>% 
  add_object(cone(start = c(-1,0,0), direction = c(-0.5,0.5,0), material = gold_mat)) %>% 
  add_object(cone(start = c(1,0,0), direction = c(0.5,0.5,0), material = gold_mat)) %>% 
  add_object(cone(start = c(0,0,-1), direction = c(0,0.5,-0.5), material = gold_mat)) %>% 
  add_object(cone(start = c(0,0,1), direction = c(0,0.5,0.5), material = gold_mat)) %>% 
  add_object(sphere(y=5,material=light())) %>% 
  add_object(sphere(y=3,x=-3,z=-3,material=light(color="red"))) %>% 
  add_object(sphere(y=3,x=3,z=-3,material=light(color="green"))) %>% 
  render_scene(lookfrom=c(0,4,10), clamp_value=10, samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render directional cone in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Render the position from the base, instead of the center of the cone:
noise_mat = material = glossy(color="purple",noisecolor="blue", noise=5)
generate_studio() %>% 
  add_object(cone(start = c(0,-1,0), from_center = FALSE, radius=1, direction = c(0,2,0), 
                  material = noise_mat)) %>% 
  add_object(cone(start = c(-1.5,-1,0), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
                  material = noise_mat)) %>% 
  add_object(cone(start = c(1.5,-1,0), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
                  material = noise_mat)) %>% 
  add_object(cone(start = c(0,-1,1.5), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
                  material = noise_mat)) %>% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
  render_scene(lookfrom=c(0,4,10), clamp_value=10,fov=25, samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render base-centered cone in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1



#Draw a simple arrow from x = -1 to x = 1

generate_studio() %>% 
  add_object(arrow(start = c(-1,0,0), end = c(1,0,0), material=glossy(color="red"))) %>% 
  add_object(sphere(y=5,material=light(intensity=20))) %>% 
  render_scene(clamp_value=10,  samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render arrow in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Change the proportion of tail to top
generate_studio(depth=-2) %>% 
  add_object(arrow(start = c(-1,-1,0), end = c(1,-1,0), tail_proportion = 0.5,
                   material=glossy(color="red"))) %>% 
  add_object(arrow(start = c(-1,0,0), end = c(1,0,0), tail_proportion = 0.75,
                   material=glossy(color="red"))) %>% 
  add_object(arrow(start = c(-1,1,0), end = c(1,1,0), tail_proportion = 0.9,
                   material=glossy(color="red"))) %>% 
  add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
  render_scene(clamp_value=10, fov=25,  samples=test_samples)%>% sum() ->
  image_sums[[counter]]
test_that("Render custom tail arrow in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Change the radius of the tail/top segments
generate_studio(depth=-1.5) %>% 
  add_object(arrow(start = c(-1,-1,0), end = c(1,-1,0), tail_proportion = 0.75,
                   radius_top = 0.1, radius_tail=0.03,
                   material=glossy(color="red"))) %>% 
  add_object(arrow(start = c(-1,0,0), end = c(1,0,0), tail_proportion = 0.75,
                   radius_top = 0.2, radius_tail=0.1,
                   material=glossy(color="red"))) %>% 
  add_object(arrow(start = c(-1,1,0), end = c(1,1,0), tail_proportion = 0.75,
                   radius_top = 0.3, radius_tail=0.2,
                   material=glossy(color="red"))) %>% 
  add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
  render_scene(clamp_value=10, samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render custom radius arrow in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#We can also specify arrows via a midpoint and direction:
generate_studio(depth=-1) %>% 
  add_object(arrow(start = c(-1,-0.5,0), direction = c(0,0,1),
                   material=glossy(color="green"))) %>% 
  add_object(arrow(start = c(1,-0.5,0), direction = c(0,0,-1),
                   material=glossy(color="red"))) %>% 
  add_object(arrow(start = c(0,-0.5,1), direction = c(1,0,0),
                   material=glossy(color="yellow"))) %>% 
  add_object(arrow(start = c(0,-0.5,-1), direction = c(-1,0,0),
                   material=glossy(color="purple"))) %>% 
  add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
  render_scene(clamp_value=10, samples = test_samples, 
               lookfrom=c(0,5,10), lookat=c(0,-0.5,0), fov=16) %>% sum() ->
  image_sums[[counter]]
test_that("Render custom midpoint arrow in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#Plot a 3D vector field for a gravitational well:

r = 1.5
theta_vals = seq(0,2*pi,length.out = 16)[-16]
phi_vals = seq(0,pi,length.out = 16)[-16][-1]
arrow_list = list()
acounter = 1
for(theta in theta_vals) {
  for(phi in phi_vals) {
    rval = c(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta)) 
    arrow_list[[acounter]] = arrow(rval, direction = -1/2*rval/sqrt(sum(rval*rval))^3,
                                  tail_proportion = 0.66, radius_top=0.03, radius_tail=0.01,
                                  material = diffuse(color="red"))
    acounter = acounter + 1
  }
}
vector_field = do.call(rbind,arrow_list)
sphere(material=diffuse(noise=1,color="blue",noisecolor="darkgreen")) %>% 
  add_object(vector_field) %>% 
  add_object(sphere(y=0,x=10,z=5,material=light(intensity=200))) %>% 
  render_scene(fov=20, ambient=TRUE, samples=test_samples,
               backgroundlow="black",backgroundhigh="white")  %>% sum() ->
  image_sums[[counter]]
test_that("Render gravity arrow in a studio", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Manually create a polygon object, here a star:

angles = seq(0,360,by=36)
xx = rev(c(rep(c(1,0.5),5),1) * sinpi(angles/180))
yy = rev(c(rep(c(1,0.5),5),1) * cospi(angles/180))
star_polygon = data.frame(x=xx,y=yy)


generate_ground(depth=0,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(star_polygon,top=0.5,bottom=0,
                              material=diffuse(color="red",sigma=90))) %>%
  add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
  render_scene(parallel=TRUE,lookfrom = c(0,2,3),samples=test_samples,lookat=c(0,0.5,0),fov=60) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded polygon star", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Now, let's add a hole to the center of the polygon. We'll make the polygon
#hollow by shrinking it, combining it with the normal size polygon,
#and specify with the `holes` argument that everything after `nrow(star_polygon)`
#in the following should be used to draw a hole:

hollow_star = rbind(star_polygon,0.8*star_polygon)


generate_ground(depth=-0.01,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon) + 1,
                              material=diffuse(color="red",sigma=90))) %>%
  add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
  render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=test_samples,lookat=c(0,0,0),fov=30) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded polygon star w/ hole", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


# Render one in the y-x plane as well by changing the `plane` argument,
# as well as offset it slightly.

generate_ground(depth=-0.01,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon),
                              material=diffuse(color="red",sigma=90))) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, z=-1.2, 
                              holes = nrow(star_polygon) + 1, plane = "yx", 
                              material=diffuse(color="green",sigma=90))) %>%
  add_object(sphere(y=4,x=-3,material=light(intensity=30))) %>%
  render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=test_samples,lookat=c(0,0.9,0),fov=40) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded polygon stars, 2 planes", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


# Now add the zy plane:

generate_ground(depth=-0.01,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon) + 1,
                              material=diffuse(color="red",sigma=90))) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, z=-1.2, 
                              holes = nrow(star_polygon) + 1, plane = "yx", 
                              material=diffuse(color="green",sigma=90))) %>%
  add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, x=1.2, 
                              holes = nrow(star_polygon) + 1, plane = "zy", 
                              material=diffuse(color="blue",sigma=90))) %>%
  add_object(sphere(y=4,x=-3,material=light(intensity=30))) %>%
  render_scene(parallel=TRUE,lookfrom = c(-4,2,4),samples=test_samples,lookat=c(0,0.9,0),fov=40) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded polygon stars, 3 planes", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#We can also directly pass in sf polygons:
if(length(find.package("spData",quiet=TRUE)) > 0) {
  us_states = spData::us_states
  texas = us_states[us_states$NAME == "Texas",]
  #Fix no sfc class in us_states geometry data
  class(texas$geometry) = c("list","sfc")
}

#This uses the raw coordinates, unless `center = TRUE`, which centers the bounding box
#of the polygon at the origin.

generate_ground(depth=-0.01,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(texas, center = TRUE,
                              material=diffuse(color="#ff2222",sigma=90))) %>%
  add_object(sphere(y=30,x=-30,radius=10,
                    material=light(color="lightblue",intensity=40))) %>%
  render_scene(parallel=TRUE,lookfrom = c(0,10,-10),samples=test_samples,fov=60) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded sf polygon", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Here we use the raw coordinates, but offset the polygon manually.

generate_ground(depth=-0.01,
                material = diffuse(color="grey50",checkercolor="grey20")) %>%
  add_object(extruded_polygon(us_states, x=-96,z=-40, top=2,
                              material=diffuse(color="#ff2222",sigma=90))) %>%
  add_object(sphere(y=30,x=-100,radius=10,
                    material=light(color="lightblue",intensity=200))) %>%
  add_object(sphere(y=30,x=100,radius=10,
                    material=light(color="orange",intensity=200))) %>%
  render_scene(parallel=TRUE,lookfrom = c(0,120,-120),samples=test_samples,fov=20) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded sf polygon, raw coords", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#We can also set the map the height of each polygon to a column in the sf object,
#scaling it down by the maximum population state.


generate_ground(depth=0,
                material = diffuse(color="grey50",checkercolor="grey20",sigma=90)) %>%
  add_object(extruded_polygon(us_states, x=-96,z=-45, data_column_top = "total_pop_15",
                              scale_data = 1/max(us_states$total_pop_15)*5,
                              material=diffuse(color="#ff2222",sigma=90))) %>%
  add_object(sphere(y=30,x=-100,z=60,radius=10,
                    material=light(color="lightblue",intensity=250))) %>%
  add_object(sphere(y=30,x=100,z=-60,radius=10,
                    material=light(color="orange",intensity=250))) %>%
  render_scene(parallel=TRUE,lookfrom = c(-60,50,-40),lookat=c(0,-5,0),samples=test_samples,fov=30) %>% sum() ->
  image_sums[[counter]]
test_that("Render extruded sf polygon, height mapped", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate the default curve:

generate_studio(depth=-0.1) %>%
  add_object(bezier_curve(material=diffuse(color="red"))) %>%
  add_object(sphere(y=3,z=5,x=2,radius=0.3,
                    material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
  render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
               samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render bezier curve", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Change the control points to change the direction of the curve. Here, we place spheres
#at the control point locations.
generate_studio(depth=-0.1) %>%
  add_object(bezier_curve(material=diffuse(color="red"))) %>%
  add_object(sphere(radius=0.075,material=glossy(color="green"))) %>% 
  add_object(sphere(radius=0.075,x=-1,y=0.33,material=glossy(color="green"))) %>% 
  add_object(sphere(radius=0.075,x=1,y=0.66,material=glossy(color="green"))) %>% 
  add_object(sphere(radius=0.075,y=1,material=glossy(color="green"))) %>% 
  add_object(sphere(y=3,z=5,x=2,radius=0.3,
                    material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
  render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=15,
               samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render custom bezier curve", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#We can make the curve flat (always facing the camera) by setting the type to `flat`
generate_studio(depth=-0.1) %>%
  add_object(bezier_curve(type="flat", material=glossy(color="red"))) %>%
  add_object(sphere(y=3,z=5,x=2,radius=0.3,
                    material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
  render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
               samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render flat bezier curve", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#We can also plot a ribbon, which is further specified by a start and end orientation with
#two surface normals.
generate_studio(depth=-0.1) %>%
  add_object(bezier_curve(type="ribbon", width=0.2,
                          p1 = c(0,0,0), p2 = c(0,0.33,0), p3 = c(0,0.66,0), p4 = c(0.3,1,0),
                          normal_end = c(0,0,1),
                          material=glossy(color="red"))) %>%
  add_object(sphere(y=3,z=5,x=2,radius=0.3,
                    material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
  render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
               samples=test_samples) %>% sum() ->
  image_sums[[counter]]
test_that("Render ribbon bezier curve", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Create a single curve and copy and rotate it around the y-axis to create a wavy fountain effect:
scene_curves = list()
for(i in 1:90) {
  scene_curves[[i]] = bezier_curve(p1 = c(0,0,0),p2 = c(0,5-sinpi(i*16/180),2),
                                   p3 = c(0,5-0.5 * sinpi(i*16/180),4),p4 = c(0,0,6),
                                   angle=c(0,i*4,0), type="cylinder",
                                   width = 0.1, width_end =0.1,material=glossy(color="red"))
}
all_curves = do.call(rbind, scene_curves)
generate_ground(depth=0,material=diffuse(checkercolor="grey20")) %>%
  add_object(all_curves) %>%
  add_object(sphere(y=7,z=0,x=0,material=light(intensity=100))) %>% 
  render_scene(lookfrom = c(12,20,50),samples = test_samples,
               lookat=c(0,1,0), fov=15, clamp_value = 10) %>% sum() ->
  image_sums[[counter]]
test_that("Render many bezier curves", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1




#Generate a wavy line, showing the line goes through the specified points:
wave = list(c(-2,1,0),c(-1,-1,0),c(0,1,0),c(1,-1,0),c(2,1,0))
point_mat = glossy(color="green")
generate_studio(depth=-1.5) %>% 
  add_object(path(points = wave,material=glossy(color="red"))) %>% 
  add_object(sphere(x=-2,y=1,radius=0.1,material=point_mat)) %>% 
  add_object(sphere(x=-1,y=-1,radius=0.1,material=point_mat)) %>% 
  add_object(sphere(x=0,y=1,radius=0.1,material=point_mat)) %>% 
  add_object(sphere(x=1,y=-1,radius=0.1,material=point_mat)) %>% 
  add_object(sphere(x=2,y=1,radius=0.1,material=point_mat)) %>% 
  add_object(sphere(z=5,x=5,y=5,radius=2,material=light(intensity=15))) %>% 
  render_scene(samples=test_samples, clamp_value=10,fov=30) %>% sum() ->
  image_sums[[counter]]
test_that("Render path", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Here we use straight lines by setting `straight = TRUE`:
generate_studio(depth=-1.5) %>% 
  add_object(path(points = wave,straight = TRUE, material=glossy(color="red"))) %>% 
  add_object(sphere(z=5,x=5,y=5,radius=2,material=light(intensity=15))) %>% 
  render_scene(samples=test_samples, clamp_value=10,fov=30) %>% sum() ->
  image_sums[[counter]]
test_that("Render straight path", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#We can also pass a matrix of values, specifying the x/y/z coordinates. Here,
#we'll create a random curve:
set.seed(21)
random_mat = matrix(runif(3*9)*2-1, ncol=3)
generate_studio(depth=-1.5) %>% 
  add_object(path(points=random_mat, material=glossy(color="red"))) %>% 
  add_object(sphere(y=5,radius=1,material=light(intensity=30))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render path, matrix init", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#We can ensure the curve is closed by setting `closed = TRUE`
generate_studio(depth=-1.5) %>% 
  add_object(path(points=random_mat, closed = TRUE, material=glossy(color="red"))) %>% 
  add_object(sphere(y=5,radius=1,material=light(intensity=30))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render closed path", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Finally, let's render a pretzel to show how you can render just a subset of the curve:
pretzel = list(c(-0.8,-0.5,0.1),c(0,-0.2,-0.1),c(0,0.3,0.1),c(-0.5,0.5,0.1), c(-0.6,-0.5,-0.1),
               c(0,-0.8,-0.1),
               c(0.6,-0.5,-0.1),c(0.5,0.5,-0.1), c(0,0.3,-0.1),c(-0,-0.2,0.1), c(0.8,-0.5,0.1))

#Render the full pretzel:
generate_studio() %>% 
  add_object(path(pretzel, width=0.17,  material = glossy(color="#db5b00"))) %>% 
  add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render path pretzel", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Here, we'll render only the first third of the pretzel by setting `u_max = 0.33`
generate_studio() %>% 
  add_object(path(pretzel, width=0.17, u_max=0.33, material = glossy(color="#db5b00"))) %>% 
  add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render path pretzel, partial u beginning", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Here's the last third, by setting `u_min = 0.66`
generate_studio() %>% 
  add_object(path(pretzel, width=0.17, u_min=0.66, material = glossy(color="#db5b00"))) %>% 
  add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render path pretzel, partial u end", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Here's the full pretzel, decomposed into thirds using the u_min and u_max coordinates
generate_studio() %>% 
  add_object(path(pretzel, width=0.17, u_max=0.33, x = -0.8, y =0.6,
                  material = glossy(color="#db5b00"))) %>% 
  add_object(path(pretzel, width=0.17, u_min=0.66, x = 0.8, y =0.6,
                  material = glossy(color="#db5b00"))) %>% 
  add_object(path(pretzel, width=0.17, u_min=0.33, u_max=0.66, x=0,
                  material = glossy(color="#db5b00"))) %>% 
  add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
  render_scene(samples=test_samples, clamp_value=10, lookfrom=c(0,3,10)) %>% sum() ->
  image_sums[[counter]]
test_that("Render path pretzel, partial u all", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1


#Generate a pig in the cornell box.


generate_cornell() %>%
  add_object(pig(x=555/2,z=555/2,y=120,scale=c(80,80,80), angle = c(0,135,0))) %>%
  render_scene(parallel=TRUE, samples=test_samples,clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render pig in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})
counter = counter + 1

#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` fov `40` .#> Error in render_scene_rcpp(camera_info = camera_info, scene_info = scene_info): Index out of bounds: [index=6; extent=6].
# Show the pig staring into a mirror, worried 
generate_cornell() %>%
  add_object(pig(x=555/2-70,z=555/2+50,y=120,scale=c(80,80,80),
                 angle = c(0,-40,0), emotion = "worried")) %>%
  add_object(cube(x=450,z=450,y=250, ywidth=500, xwidth=200,
                  angle = c(0,45,0), material = metal())) %>%
  render_scene(parallel=TRUE, samples=test_samples,clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render pig in cornell box w/ mirror", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` fov `40` .#> Error in render_scene_rcpp(camera_info = camera_info, scene_info = scene_info): Index out of bounds: [index=6; extent=6].
# Render many small pigs facing random directions, with an evil pig overlord
set.seed(1)
lots_of_pigs = list() 
for(i in 1:10) {
  lots_of_pigs[[i]] = pig(x=50 + 450 * runif(1), z = 50 + 450 * runif(1), y=50, 
                          scale = c(30,30,30), angle = c(0,360*runif(1),0), emotion = "worried")
}

many_pigs_scene = do.call(rbind, lots_of_pigs) %>%
  add_object(generate_cornell(lightintensity=30, lightwidth=100)) %>%
  add_object(pig(z=500,x=555/2,y=400, emotion = "angry",
                 scale=c(100,100,100),angle=c(30,90,0), order_rotation=c(2,1,3))) 

render_scene(many_pigs_scene,parallel=TRUE,clamp_value=10, samples=test_samples, 
             sample_method="stratified")%>% sum() ->
  image_sums[[counter]]
test_that("Render many pigs in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})
counter = counter + 1

#> Setting default values for Cornell box: lookfrom `c(278,278,-800)` lookat `c(278,278,0)` fov `40` .#> Error in render_scene_rcpp(camera_info = camera_info, scene_info = scene_info): Index out of bounds: [index=6; extent=6].
#Render spiderpig
generate_studio() %>%  
  add_object(pig(y=-1,angle=c(0,-100,0), scale=1/2,spider=TRUE)) %>% 
  add_object(sphere(y=5,z=5,x=5,material=light(intensity=100))) %>% 
  render_scene(samples=test_samples,lookfrom=c(0,2,10),clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render spiderpig in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Generate a label in the cornell box.

generate_cornell() %>% 
  add_object(text3d(label="Cornell Box", x=555/2,y=555/2,z=555/2,text_height=60,
                    material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render label in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1
#Change the orientation
generate_cornell() %>% 
  add_object(text3d(label="YZ Plane", x=550,y=555/2,z=555/2,text_height=100,
                    orientation = "yz",
                    material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
  add_object(text3d(label="XY Plane", z=550,y=555/2,x=555/2,text_height=100,
                    orientation = "xy",
                    material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
  add_object(text3d(label="XZ Plane", z=555/2,y=5,x=555/2,text_height=100,
                    orientation = "xz",
                    material=diffuse(color="grey10"))) %>% 
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render label (diff planes) in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

#Add an label in front of a sphere
generate_cornell() %>% 
  add_object(text3d(label="Cornell Box", x=555/2,y=555/2,z=555/2,text_height=60,
                    material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
  add_object(text3d(label="Sphere", x=555/2,y=100,z=100,text_height=30,
                    material=diffuse(color="white"), angle=c(0,180,0))) %>% 
  add_object(sphere(y=100,radius=100,z=555/2,x=555/2,
                    material=glossy(color="purple"))) %>% 
  add_object(sphere(y=555,radius=100,z=-1000,x=555/2,
                    material=light(intensity=100,
                                   spotlight_focus=c(555/2,100,100)))) %>%                   
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render label in front of sphere in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

set.seed(1)
#A room full of bees
bee_list = list()
for(i in 1:100) {
  bee_list[[i]] = text3d("B", x=20+runif(1)*525, y=20+runif(1)*525, z=20+runif(1)*525, 
                         text_height = 50, angle=c(0,180,0))
}
bees = do.call(rbind,bee_list)

generate_cornell() %>% 
  add_object(bees) %>%                   
  render_scene(samples=test_samples, clamp_value=10) %>% sum() ->
  image_sums[[counter]]
test_that("Render B labels in cornell box", {expect_equal(previous_sums[[counter]], image_sums[[counter]])})

counter = counter + 1

## Note for contributors:
## If your contributions result in intended changes that cause failures to some of these tests,
## Re-run the tests with the below line un-commented. Place the resulting file in the inst/testdata/
## folder, and include a note in the pull request that documents the change and why the change in 
## behavior is intended.

# writeLines(as.character(image_sums),"test-object_image_sums.txt")

