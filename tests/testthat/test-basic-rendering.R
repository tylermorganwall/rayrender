save_test_png = function(code, path) {
  code
  path
}

compare_image = function(path1, path2) {
  image1 = png::readPNG(path1)
  image2 = png::readPNG(path2)
  return(identical(image1, image2))
}

run_tests = function(func, argument_grid, plot_prefix="", ...) {
  stopifnot(inherits(argument_grid,"data.frame"))
  
  for(i in seq_len(nrow(argument_grid))){
    set.seed(123)
    args = unlist(argument_grid[i,], recursive = FALSE)
    test_filename = paste0(sprintf("%s-%s", names(args), args),collapse="_")
    test_filename = sprintf("%s_test%i_%s.png",
                            plot_prefix , i, test_filename)
    path = tempfile(fileext = ".png")
    
    args = append(args, ...)
    args = append(args, list(filename = path))
    save_test_png(do.call(func, args = args), path) |> 
      suppressMessages() |> 
      suppressWarnings() |> 
      expect_snapshot_file(name = test_filename, compare = compare_image)
  }
}

test_that("plot_3d plots basic options", {
  skip_on_cran()
  
  scene1 = generate_ground(depth=-0.5, material = diffuse(color="white", checkercolor="darkgreen")) |> 
    add_object(sphere(x=0,y=0,z=0,radius=0.5,material = diffuse(color=c(1,0,1)))) |> 
    add_object(cube(x=1.1,y=0,z=0,material = diffuse(noise=3))) |> 
    add_object(cylinder(x=-1.1,y=0,z=0,radius=0.5,material = metal(color="gold",fuzz=0.02))) |> 
    add_object(ellipsoid(x=-1.1,y=1.5, a=0.5,b=0.75,c=0.5,
                       material=microfacet(roughness=0.1, 
                                           eta=c(0.216,0.42833,1.3184), 
                                           kappa=c(3.239,2.4599,1.8661)))) |> 
    add_object(cone(start=c(0,0.7,0),end=c(0,2.1,0),radius=0.5,material = glossy(
      color = "green", 
      gradient_point_start = c(0,0.7,0), gradient_point_end = c(0,2.1,0),
      gradient_color = "red"))) |> 
    add_object(pig(y=0.5,x=1,scale = 0.3,angle=c(0,-70,0))) |> 
    add_object(sphere(x=10,y=5,z=10,radius=2,material = light(intensity = 40))) |> 
    add_object(sphere(x=-0.4,y=2,z=-3,radius=0.02,material = light(intensity = 10000))) 
  old_cores = options("cores")
  on.exit(options("cores" = old_cores))
  num_cores = parallel::detectCores()
  if(num_cores > 2) {
    options("cores" = num_cores - 1)
  }
  
  render_args_basic = expand.grid(width         = list(50,100),
                                  height        = list(50,100),
                                  clamp_value   = list(Inf, 10),
                                  parallel      = list(TRUE, FALSE))
  
  run_tests("render_scene", render_args_basic, plot_prefix = "render_basic", 
            list(scene = scene1,
                 samples = 16,
                 ortho_dimensions = c(3,3),
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_args_camera_info = expand.grid(lookfrom = list(c(0,1,10),c(-10,1,10)),
                                        lookat   = list(c(0,0,0), c(0,0.75,0)),
                                        camera_up = list(c(0,1,0), c(1,1,0)),
                                        focal_distance = list(NULL, 5),
                                        fov           = list(0,40,360))
  
  run_tests("render_scene", render_args_camera_info, plot_prefix = "render_camera_info", 
            list(scene = scene1,
                 samples = 16,
                 ortho_dimensions = c(3,3),
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_args_post_processing_info = expand.grid(bloom       = list(FALSE, 10),
                                                 clamp_value = list(Inf, 10),
                                                 tonemap     = list("gamma","reinhold",
                                                                    "uncharted","hbd"))
  
  run_tests("render_scene", render_args_post_processing_info, plot_prefix = "render_post", 
            list(scene = scene1,
                 samples = 16,
                 ortho_dimensions = c(3,3),
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_args_debug = expand.grid(debug_channel= list("none","depth","normals","uv",
                                                      "variance","dpdu","dpdv","color"),
                                                 return_raw_array = list(TRUE, FALSE))
  
  run_tests("render_scene", render_args_debug, plot_prefix = "render_debug", 
            list(scene = scene1,
                 samples = 16,
                 ortho_dimensions = c(3,3),
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_sampling_args = expand.grid(samples = list(1,16),
                                     sample_method = list("sobol","sobol_blue","random", "stratified"),
                                     min_variance  = list(0, 1e-5),
                                     min_adaptive_size = list(1,11),
                                     ambient = list(TRUE, FALSE))
  
  run_tests("render_scene", render_sampling_args, plot_prefix = "render_sampling", 
            list(scene = scene1,
                 width = 200,
                 height = 200,
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_termination_args = expand.grid(max_depth     = list(NA, 10, 100),
                                        roulette_active_depth = list(100,10,3))
  
  run_tests("render_scene", render_termination_args, plot_prefix = "render_termination", 
            list(scene = scene1,
                 samples = 16,
                 width = 200,
                 height = 200,
                 preview  = FALSE, 
                 progress = FALSE))
  
  render_realistic_camera_args = expand.grid(camera_description_file = list("telephoto", "50mm","wide","fisheye"),
                                             aperture = list(1,10),
                                             camera_scale = list(1,3),
                                             film_size = list(22,40))
  #Only include a few examples of 40mm film
  render_realistic_camera_args = subset(render_realistic_camera_args, 
                                        film_size == 22 | (camera_scale == 1 & aperture == 10))
  
  run_tests("render_scene", render_realistic_camera_args, plot_prefix = "render_realistic_camera", 
            list(scene = scene1,
                 samples = 16,
                 iso = 10000,
                 width = 200,
                 height = 200,
                 preview  = FALSE, 
                 progress = FALSE))
                                     
  
  expect_error(render_scene(scene1,lookat=c(0,0.5,0), width  = 3, bloom=TRUE))
  expect_error(render_scene(scene1,lookat=c(0,0.5,0), height = 3, bloom=TRUE))
})
