# Render the authored Blender models with rayrender. Source after library(rayrender).
# Native node, morph, and skeletal clips share exact-time preview/offline evaluation.

#' @param asset_dir Directory containing bird/bird.glb from the Blender generator.
#' @return Shared asset, two independently timed instances, and their lit scene.
#' @keywords internal
blender_birds_scene = function(asset_dir) {
  imported = animation_read_gltf(file.path(asset_dir, "bird", "bird.glb"))
  # The current importer uses explicit rayrender materials. Reassign by authored
  # part name so the cream breast, dark feathers, and eyes remain distinguishable.
  nodes = lapply(imported$nodes, function(node) {
    if (is.null(node$scene)) {
      return(node)
    }
    color = "#167C9B"
    if (grepl("Primary|Secondary|Tail", node$name)) {
      color = "#154A66"
    }
    if (grepl("Covert", node$name)) {
      color = "#52B4C4"
    }
    if (grepl("Cream|Throat|Eye ring|Eye glint", node$name)) {
      color = "#F3DDB4"
    }
    if (grepl("^Eye(\\.[0-9]+)?$", node$name)) {
      color = "#101820"
    }
    if (grepl("beak|foot", node$name)) {
      color = "#D48228"
    }
    node$scene = set_scene_material(node$scene, diffuse(color))
    node
  })
  asset = animation_model_asset(
    nodes,
    imported$clips,
    name = "Blender kingfisher"
  )
  first = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, 1, 2, 3, 4),
      position = rbind(
        c(-1.7, 1.7, 0),
        c(-1.9, 2.1, -.4),
        c(-1.6, 2.2, .2),
        c(-1.4, 1.8, .5),
        c(-1.7, 1.7, 0)
      ),
      angle = rbind(
        c(0, -20, -8),
        c(-5, 10, 10),
        c(5, 25, 5),
        c(8, -5, -12),
        c(0, -20, -8)
      ),
      position_interpolation = "cubic",
      rotation_interpolation = "smooth"
    ),
    animation_clip_player("clip-0", loop = "repeat")
  )
  second = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, .8, 1.9, 3.1, 4),
      position = rbind(
        c(1.7, 2.3, -.5),
        c(1.9, 1.8, .3),
        c(1.5, 1.7, -.3),
        c(1.4, 2.2, -.8),
        c(1.7, 2.3, -.5)
      ),
      angle = rbind(
        c(0, 20, 8),
        c(8, -10, -10),
        c(-5, -25, -5),
        c(-8, 5, 12),
        c(0, 20, 8)
      ),
      position_interpolation = "cubic",
      rotation_interpolation = "smooth"
    ),
    animation_clip_player("clip-0", offset = .65, speed = 1.5, loop = "repeat")
  )
  scene = add_object(first, second) |>
    add_object(rayrender_cycle_stage()) |>
    add_camera(camera(
      lookfrom = c(4, 3.3, 7),
      lookat = c(0, 1.5, 0),
      fov = 44,
      aperture = 0
    ))
  list(asset = asset, first = first, second = second, scene = scene)
}

#' @return A shared ground plane and studio area lights for the rendered examples.
#' @keywords internal
rayrender_cycle_stage = function() {
  xz_rect(
    y = -.015,
    xwidth = 200,
    zwidth = 200,
    material = diffuse("#718591", checkercolor = "#647783", checkerperiod = 1)
  ) |>
    add_object(sphere(
      x = -6,
      y = 12,
      z = 8,
      radius = 2,
      material = light("#FFF2DE", intensity = 32)
    )) |>
    add_object(sphere(
      x = 6,
      y = 10,
      z = -6,
      radius = 1.8,
      material = light("#BCDFFF", intensity = 24)
    ))
}

#' @param asset_dir Directory containing jellyfish/jellyfish.glb.
#' @return Shared morph asset, independently phased instances, and their scene.
#' @keywords internal
blender_jellyfish_scene = function(asset_dir) {
  imported = animation_read_gltf(file.path(
    asset_dir,
    'jellyfish',
    'jellyfish.glb'
  ))
  nodes = lapply(imported$nodes, function(node) {
    if (is.null(node$scene)) {
      return(node)
    }
    color = '#9578CC'
    if (grepl('rim|canal', node$name)) {
      color = '#53C9D9'
    }
    if (grepl('Oral', node$name)) {
      color = '#DC8499'
    }
    if (grepl('Tentacle', node$name)) {
      color = '#BEA6DF'
    }
    node$scene = set_scene_material(node$scene, glossy(color, gloss = .55))
    node
  })
  asset = animation_model_asset(
    nodes,
    imported$clips,
    name = 'Blender jellyfish'
  )
  first = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, 1, 2, 3, 4),
      position = rbind(
        c(-1.1, 2.2, .1),
        c(-1.4, 2.4, -.2),
        c(-1.1, 2.5, -.4),
        c(-.9, 2.3, -.2),
        c(-1.1, 2.2, .1)
      ),
      angle = rbind(
        c(0, 0, -5),
        c(5, 20, 5),
        c(0, 40, 5),
        c(-5, 15, -5),
        c(0, 0, -5)
      ),
      position_interpolation = 'cubic',
      rotation_interpolation = 'smooth'
    ),
    animation_clip_player('clip-0', loop = 'repeat')
  )
  second = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, .7, 1.7, 2.9, 4),
      position = rbind(
        c(1.2, 2.7, -.6),
        c(1.4, 2.4, -.2),
        c(1.1, 2.3, .2),
        c(.9, 2.5, -.3),
        c(1.2, 2.7, -.6)
      ),
      angle = rbind(
        c(5, 15, 5),
        c(-5, 0, -5),
        c(0, -30, -5),
        c(5, 0, 5),
        c(5, 15, 5)
      ),
      position_interpolation = 'cubic',
      rotation_interpolation = 'smooth'
    ),
    animation_clip_player('clip-0', offset = .65, speed = 1.5, loop = 'repeat')
  )
  scene = add_object(first, second) |>
    add_object(rayrender_cycle_stage()) |>
    add_camera(camera(
      lookfrom = c(5, 3.7, 8),
      lookat = c(0, 1.6, 0),
      fov = 40,
      aperture = 0
    ))
  list(asset = asset, first = first, second = second, scene = scene)
}

#' @param asset_dir Directory containing walk/walk.glb.
#' @return Shared skin asset, separately phased walkers on two paths, and scene.
#' @keywords internal
blender_walk_scene = function(asset_dir) {
  asset = animation_read_gltf(
    file.path(asset_dir, 'walk', 'walk.glb'),
    recalculate_normals = TRUE,
    material = glossy('#3689AC', gloss = .55)
  )
  # The authored 0.75-unit steps take one second. Advance placement at that
  # speed while the internal clip loops, preserving heel contact.
  first = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, 1, 2, 3, 4),
      position = rbind(
        c(-.9, 0, -1.5),
        c(-.75, 0, -.75),
        c(-.7, 0, 0),
        c(-.85, 0, .75),
        c(-.9, 0, 1.5)
      ),
      angle = rbind(
        c(0, 10, 0),
        c(0, 6, 0),
        c(0, -4, 0),
        c(0, -10, 0),
        c(0, 0, 0)
      ),
      position_interpolation = 'cubic',
      rotation_interpolation = 'smooth'
    ),
    animation_clip_player('clip-0', loop = 'repeat')
  )
  second = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, .8, 1.8, 3, 4),
      position = rbind(
        c(.9, 0, -1.8),
        c(1.1, 0, -1.2),
        c(1, 0, -.45),
        c(1.05, 0, .45),
        c(1.2, 0, 1.2)
      ),
      angle = rbind(
        c(0, 10, 0),
        c(0, 0, 0),
        c(0, 0, 0),
        c(0, 8, 0),
        c(0, 10, 0)
      ),
      position_interpolation = 'cubic',
      rotation_interpolation = 'smooth'
    ),
    animation_clip_player('clip-0', offset = .7, loop = 'repeat')
  )
  scene = add_object(first, second) |>
    add_object(rayrender_cycle_stage()) |>
    add_camera(camera(
      lookfrom = c(3.5, 2.5, 5.3),
      lookat = c(0, 1, .3),
      fov = 38,
      aperture = 0
    ))
  list(asset = asset, first = first, second = second, scene = scene)
}

#' @param asset_dir Directory containing the three authored GLB subdirectories.
#' @param output Directory for PNG frames, MP4 movies, and timing reports.
#' @param fps Default `12`. Output samples per scene second.
#' @param duration Default `4`. Duration in scene seconds; times exclude the endpoint.
#' @param width Default `640`. Frame width in pixels.
#' @param height Default `480`. Frame height in pixels.
#' @param samples Default `64`. Path samples per pixel.
#' @param integrator_type Default `"metal"`. Rayrender backend, or nee for CPU.
#' @return Invisibly the example render reports. Requires ffmpeg for MP4 encoding.
#' @keywords internal
render_blender_model_examples = function(
  asset_dir,
  output,
  fps = 12,
  duration = 4,
  width = 640,
  height = 480,
  samples = 64,
  integrator_type = 'metal'
) {
  ffmpeg = Sys.which('ffmpeg')
  if (!nzchar(ffmpeg)) {
    stop('Install ffmpeg to encode the example movies.')
  }
  if (
    !is.finite(fps) ||
      fps <= 0 ||
      !is.finite(duration) ||
      duration <= 0 ||
      fps * duration != trunc(fps * duration)
  ) {
    stop('fps * duration must be a positive integer.')
  }
  times = (seq_len(fps * duration) - 1) / fps
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  asset_dir = path.expand(asset_dir)
  output = normalizePath(output, mustWork = TRUE)
  reports = list()
  builders = list(
    bird = blender_birds_scene,
    jellyfish = blender_jellyfish_scene,
    walk = blender_walk_scene
  )
  for (name in names(builders)) {
    directory = file.path(output, name)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    compile = system.time(example <- builders[[name]](asset_dir))
    # Check every requested time cheaply before beginning the full-size movie.
    # This catches device precision/fallback issues without discarding a long run.
    if (integrator_type == 'metal') {
      message('Checking all requested Metal times for ', name)
      probe = render_animation(
        example$scene,
        times = times,
        fps = fps,
        width = 8,
        height = 8,
        samples = 1,
        max_depth = 2,
        preview = FALSE,
        progress = FALSE,
        plot_scene = FALSE,
        denoise = FALSE,
        bloom = FALSE,
        ambient_light = FALSE,
        integrator_type = integrator_type
      )
      probe_backend = lapply(probe, function(frame) attr(frame, 'wavefront'))
      saveRDS(probe_backend, file.path(directory, 'preflight-backend.rds'))
      if (
        !all(vapply(
          probe_backend,
          function(report) isTRUE(report$used),
          logical(1)
        ))
      ) {
        stop(
          'Metal preflight failed; inspect preflight-backend.rds and fallback warnings.'
        )
      }
      rm(probe)
    }
    message('Rendering ', name, ': ', length(times), ' frames')
    # Each frame asks for an absolute scene time. No Blender process, mesh cache,
    # or sequential pose integration participates in these renders.
    elapsed = system.time(
      frames <- render_animation(
        example$scene,
        times = times,
        fps = fps,
        width = width,
        height = height,
        samples = samples,
        preview = FALSE,
        progress = FALSE,
        plot_scene = FALSE,
        denoise = TRUE,
        bloom = FALSE,
        tonemap = 'hbd',
        ambient_light = FALSE,
        backgroundhigh = '#324453',
        backgroundlow = '#9AB3C3',
        integrator_type = integrator_type
      )
    )
    backend = lapply(frames, function(frame) attr(frame, 'wavefront'))
    for (i in seq_along(frames)) {
      rayimage::ray_write_image(
        frames[[i]],
        filename = file.path(directory, sprintf('frame-%04d.png', i))
      )
    }
    saveRDS(backend, file.path(directory, 'frame-backend.rds'))
    if (
      integrator_type == 'metal' &&
        !all(vapply(backend, function(report) isTRUE(report$used), logical(1)))
    ) {
      stop(
        'Metal did not render every requested frame; inspect fallback warnings before continuing.'
      )
    }
    status = system2(
      ffmpeg,
      c(
        '-v',
        'error',
        '-y',
        '-framerate',
        fps,
        '-start_number',
        '1',
        '-i',
        shQuote(file.path(directory, 'frame-%04d.png')),
        '-frames:v',
        length(times),
        '-c:v',
        'libx264',
        '-crf',
        '18',
        '-pix_fmt',
        'yuv420p',
        '-movflags',
        '+faststart',
        shQuote(file.path(output, paste0(name, '.mp4')))
      )
    )
    if (status != 0L) {
      stop('ffmpeg failed for ', name)
    }
    reports[[name]] = list(
      times = times,
      clips = animation_list_clips(example$asset),
      construction_seconds = unname(compile['elapsed']),
      render_seconds = unname(elapsed['elapsed']),
      frames = length(frames),
      width = width,
      height = height,
      samples = samples,
      backend = backend
    )
    rm(frames, example)
    gc(FALSE)
    message(
      'Rendered ',
      name,
      ' through rayrender: ',
      length(times),
      ' exact-time frames.'
    )
  }
  inputs = unlist(lapply(names(builders), function(name) {
    c('-i', shQuote(file.path(output, paste0(name, '.mp4'))))
  }))
  status = system2(
    ffmpeg,
    c(
      '-v',
      'error',
      '-y',
      inputs,
      '-filter_complex',
      shQuote('[0:v][1:v][2:v]hstack=inputs=3[v]'),
      '-map',
      shQuote('[v]'),
      '-c:v',
      'libx264',
      '-crf',
      '18',
      '-pix_fmt',
      'yuv420p',
      '-movflags',
      '+faststart',
      shQuote(file.path(output, 'rayrender-models.mp4'))
    )
  )
  if (status != 0L) {
    stop('ffmpeg failed for the combined movie.')
  }
  saveRDS(reports, file.path(output, 'render-report.rds'))
  jsonlite::write_json(
    reports,
    file.path(output, 'render-report.json'),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  invisible(reports)
}
