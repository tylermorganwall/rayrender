#' Location and Time Sky Light
#' @md
#'
#' @description
#' Creates a location-and-time sky light. Add it with [add_infinite_light()].
#' By default, `skymodelr::generate_sky_latlong()` generates a cached EXR.
#' Set `atmosphere = TRUE` to evaluate the Prague model throughout the scene,
#' including finite-distance atmospheric attenuation and in-scattering by default.
#'
#' @param lat Latitude in degrees, between -90 and 90.
#' @param long Longitude in degrees, between -180 and 180.
#' @param datetime A single `POSIXct` date and time. Specify its time zone when
#' constructing it with `as.POSIXct()`.
#' @param altitude Default `0`. Reference altitude in meters above sea level.
#' Prague supports 0--15000 m. With native atmosphere this is the altitude at
#' `atmosphere_origin`; with an image sky it selects one observer for the image.
#' @param visibility Default `50`. Prague meteorological visibility in kilometers,
#' from 20 to 131.8. Smaller values produce stronger haze.
#' @param albedo Default `0.5`. Uniform ground reflectance for the sky model,
#' between 0 and 1. Local surface materials are specified separately.
#' @param resolution Default `if (isTRUE(atmosphere)) 64 else 2048`. Image height
#' in pixels for a cached sky (width is twice the height). With native atmosphere,
#' the height of directional importance-sampling tables, from 16 to 2048;
#' this does not limit the rendered sky's detail.
#' @param hosek Default `!isTRUE(atmosphere)`. Use the Hosek model for an image
#' sky. Set `FALSE` to use Prague. Must be `FALSE` with native atmosphere.
#' @param render_mode Default `"all"`. Render sky and Sun (`"all"`), sky without
#' the solar disk (`"atmosphere"`), or the solar disk alone (`"sun"`).
#' @param turbidity Default `3`. Hosek turbidity, from 1.7 to 10. Image skies only.
#' @param wide_spectrum Default `FALSE`. Use Prague's 55-channel sea-level
#' coefficients for an image sky. Must be `FALSE` with native atmosphere.
#' @param below_horizon Default `TRUE`. Include atmospheric radiance below the
#' horizon. Must be `TRUE` with native atmosphere for downward queries.
#' @param prague_rgb_correction Default `TRUE`. Apply skymodelr's Prague RGB
#' tint correction. Applies when `hosek = FALSE`.
#' @param prague_rgb_correction_strength Default `1`. Strength of the Prague RGB
#' tint correction. Must be finite and nonnegative: 0 disables correction,
#' and 1 applies the full calibrated correction.
#' @param prague_rgb_correction_gain Default `"auto"`. Calibrated Prague RGB
#' gains, or a numeric vector of three finite, positive linear RGB multipliers.
#' @param stars Default `FALSE`. Composite stars into an image sky.
#' Must be `FALSE` with native atmosphere.
#' @param star_width Default `1`. Stellar point-spread size, passed to
#' `skymodelr::generate_stars()`. Image skies only.
#' @param stars_exposure Default `0`. Artistic exposure adjustment for stars,
#' in stops. Image skies only.
#' @param planets Default `FALSE`. Composite bright planets into an image sky.
#' Must be `FALSE` with native atmosphere.
#' @param moon Default `FALSE`. Composite a Moon image into the cached sky.
#' Must be `FALSE` with native atmosphere; use [moon_light()] for a separate disk.
#' @param moon_atmosphere Default `FALSE`. Include atmospheric scattering of
#' moonlight in an image sky.
#' @param moon_hosek Default `TRUE`. Use Hosek for the image sky's moonlight
#' scattering. Set `FALSE` to use Prague.
#' @param exr_adopted_white Default `"D60"`. Adopted white for cached EXR metadata:
#' `"D60"`, `"D65"`, or numeric XYZ with Y = 1. Does not change image pixels.
#' @param exr_metadata Default `TRUE`. Attach skymodelr color metadata to the
#' cached EXR. Image skies only.
#' @param number_cores Default `1`. CPU threads used to generate a cached sky.
#' Native atmospheric rendering uses the renderer's thread settings.
#' @param verbose Default `FALSE`. Print sky-generation progress information.
#' @param ... Additional named image-sky arguments forwarded by
#' `skymodelr::generate_sky_latlong()` to its star, planet, and Moon generators.
#' Location, datetime, and the cached filename are managed by this light.
#' Native atmosphere rejects unsupported arguments.
#' @param intensity Default `1`. Nonnegative multiplier for this light's radiance.
#' @param rotation Default `0`. Additional rotation in degrees around the world
#' Y axis, using the same convention as [infinite_light()].
#' @param name Default `"sky"`. Unique light name within the scene.
#' @param atmosphere Default `FALSE`. Use the native Prague model, with
#' `attenuation` and `query_altitude` controlling finite haze and observer
#' positions. Requires `integrator_type = "nee"` and the full-altitude dataset.
#' @param attenuation Default `TRUE`. With `atmosphere = TRUE`, include
#' finite-distance attenuation and in-scattering between scene interactions.
#' Set `FALSE` to evaluate native Prague lighting without this finite-distance
#' haze. Sun/sky radiance and celestial disk filtering still include the
#' atmosphere between the query location and space.
#' @param haze_in_volumes Default `TRUE`. With native atmospheric attenuation,
#' integrate clear-air haze inside attached volume materials as well as outside.
#' Set `FALSE` to pause finite haze while a ray is inside any volume boundary,
#' resuming at its exit. This applies to camera paths, shadow connections, and
#' background opacity. The entire enclosed volume is excluded, including empty
#' cells; the volume's own scattering, absorption, and emission remain active.
#' @param deferred_haze Default `TRUE`. With native atmospheric attenuation,
#' defer haze queries between interactions and estimate changes in transport
#' weights using one reservoir sample. This reduces model queries while preserving
#' the estimator's mean, with additional Monte Carlo noise. Extinction is applied
#' at each completed span; emitting medium events also end a span. Works with
#' either `haze_in_volumes` setting. Set `FALSE` to evaluate every haze interval.
#' Signed corrections are averaged before display processing. Keep
#' `render_scene(clamp_value = Inf)` to avoid clipping that estimator.
#' @param haze_correction_probability Default `if (isTRUE(deferred_haze)) 0.5 else 1`.
#' Probability of evaluating the full horizon-smoothing correction at each
#' deferred haze endpoint. Values below
#' `1` require `deferred_haze = TRUE`. The deferred default, `0.5`, reduces sky
#' queries while preserving the expected radiance, with possible additional
#' noise. Must be a finite number
#' greater than zero and at most one. Set `1` to always evaluate the full
#' correction. Disabling `deferred_haze` defaults this probability to `1`.
#' Transmission and direct lighting stay exact.
#' Keep `clamp_value = Inf` so signed corrections are averaged without clipping.
#' @param cache_spectra Default `TRUE`. With `atmosphere = TRUE`, reuse exactly
#' matching Prague sky spectra in a small cache per rendering thread. This changes
#' neither the model nor individual samples. Set `FALSE` to disable this cache.
#' @param transmission_table Default `TRUE`. With `atmosphere = TRUE`, precompute
#' transmission reconstruction at Prague's existing grid points, preserving its
#' interpolation and individual sample results. Adds about 128 MiB per model at
#' 50 km visibility, shared across threads. Tables exceeding
#' `transmission_table_max_mb`, allocation failures, and queries outside the
#' cached visibility slices use the original compressed evaluator. Set `FALSE`
#' to retain that evaluator for all queries.
#' @param transmission_table_max_mb Default `512`. Maximum additional memory per
#' model for the transmission table, in MiB (1024^2 bytes), shared across threads.
#' Accepts nonnegative numbers, including fractional values. Set `0` to disable
#' table allocation or `Inf` to remove the cap. If the complete table does not
#' fit, use the original exact evaluator. Applies when `transmission_table = TRUE`.
#' @param query_altitude Default `TRUE`. With `atmosphere = TRUE`, query lighting
#' at each surface, cloud, or camera position, using `meters_per_unit` and
#' `atmosphere_origin`. Set `FALSE` with `attenuation = FALSE` to evaluate all
#' lighting at the fixed reference point and `altitude`. Finite-distance
#' haze requires position-dependent queries, so `attenuation = TRUE` requires
#' `query_altitude = TRUE`.
#' @param meters_per_unit Default `1`. Physical meters per world-space unit when
#' `atmosphere = TRUE`. Applies to distances along every axis.
#' @param atmosphere_origin Default `c(0, 0, 0)`. World-space location of the
#' geographic reference point, at `altitude` meters above sea level.
#' World +Y is up; horizontal offsets follow the model's spherical Earth.
#'
#' @details skymodelr is a required package dependency. Install any required
#' Prague datasets with `skymodelr::download_sky_data()` before rendering.
#' Atmospheric queries run through skymodelr's registered native API; its
#' implementation and coefficients are shared without calling R from workers.
#'
#' By default this is a static image-based sky: `altitude` selects one
#' observer altitude for the entire environment. With `atmosphere = TRUE` and
#' the default `query_altitude = TRUE`, the
#' sky and Sun change with the location of each surface or cloud interaction.
#' Date and time stay fixed during an animation in both modes.
#' The image uses skymodelr's orientation: north at the image seam, east one
#' quarter across. With zero rotation, north is world +Z and east is world -X.
#' Other infinite lights add to the sky. Adjust `iso` in [render_scene()] for
#' the sky's physical radiance scale; use the same ISO when comparing lights.
#'
#' Atmospheric mode uses the Prague model regardless of the static sky default.
#' Install its data with `skymodelr::download_sky_data(sea_level = FALSE)`.
#' A scene can have one atmospheric sky. Do not combine its haze with a medium
#' modeling the same atmospheric scattering or absorption: that would count the
#' atmosphere twice. Separate cloud volumes can be added normally.
#' The Sun is sampled as a disk independently of the sky sampling resolution.
#' Its visibility accounts for Earth's curvature, allowing high clouds to remain
#' sunlit after ground-level sunset. Atmospheric refraction is not modeled.
#' This also works with `attenuation = FALSE, query_altitude = TRUE`: a cloud
#' base several kilometres above sea level can receive direct sunlight from
#' below its local horizontal plane while Earth hides the Sun from the ground.
#' Ground surfaces still receive diffuse twilight and any indirect cloud light.
#' An image sky or `query_altitude = FALSE` shares one observer's horizon across
#' the scene and cannot reproduce this altitude-dependent lighting.
#'
#' Supported atmospheric sky settings are `altitude` (default 0 m, range
#' 0--15000), `visibility` (default 50 km, range 20--131.8), `albedo` (default
#' 0.5, range 0--1), `render_mode` (`"all"`, `"atmosphere"`, or `"sun"`),
#' and skymodelr's `prague_rgb_correction`, `prague_rgb_correction_strength`,
#' and `prague_rgb_correction_gain`. `resolution` (default 64) controls the
#' height of directional importance-sampling tables, not the rendered sky's
#' detail. `hosek`, `wide_spectrum`, `moon`, `stars`, and `planets` must be
#' `FALSE`; `below_horizon` must be `TRUE`. Solar elevations
#' outside -4.2--90 degrees are rejected. Queries outside the altitude range
#' use the nearest modeled altitude; keep scene interactions within that range.
#'
#' The model precomputes clear-air multiple scattering over a spherical Earth
#' with uniform ground albedo. Local geometry and clouds block direct Sun and
#' sky lighting, but do not cast shadows into this precomputed atmospheric
#' in-scattering field. Haze is disabled inside dielectric solids.
#' Radiance is integrated spectrally and converted to renderer RGB; attenuation
#' of RGB materials uses a broadband approximation. Finite-distance fitted
#' transmittance is normalized at zero distance and interpolated in optical
#' depth over the first 100 m. Ray-anchored cumulative transport prevents those
#' fitting errors from accumulating at each cloud null event or boundary.
#'
#' Additional image lights represent radiance outside the atmosphere and receive
#' atmospheric attenuation. Do not supply an image that already includes the
#' same haze. In atmospheric scenes, [sun_light()] and [moon_light()] automatically
#' request unattenuated textures from skymodelr. The renderer applies spectral
#' atmospheric filtering and Earth occlusion at each interaction, including
#' partial disks and the depressed horizon at altitude. An explicit Sun light
#' replaces this sky's built-in solar disk, preserving the sky and haze.
#' Without an explicit disk altitude, celestial ephemerides use this sky's
#' reference altitude. Light rotations and intensities remain independent;
#' match them when the disk should correspond to this sky's solar illumination.
#' The precomputed haze remains Sun-driven: a Moon disk illuminates surfaces and
#' clouds but does not add moonlit atmospheric in-scattering or a lunar halo.
#' With a transparent background,
#' atmospheric in-scattering remains foreground radiance, with scalar opacity
#' derived from primary-ray atmospheric transmittance. RGB transmission into
#' an arbitrary compositing background remains an approximation.
#'
#' @return A `ray_infinite_light` object containing the sky description.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' if (requireNamespace("skymodelr", quietly = TRUE)) {
#'   scene = sphere(material = diffuse(color = "white")) |>
#'     add_infinite_light(sky_light(
#'       lat = 40.7,
#'       long = -74,
#'       datetime = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York"),
#'       resolution = 256
#'     )) |>
#'     add_camera(camera(aperture = 0, fov=50))
#'   render_scene(scene, integrator_type = "nee", iso = 4)
#' }
#'
#' # Capsule hills and fluffy clouds at several distances. All scene coordinates
#' # below are kilometres; meters_per_unit = 1000 supplies the physical scale.
#' # Install the full-altitude Prague data once before running these examples:
#' # skymodelr::download_sky_data(sea_level = FALSE)
#' # cloud() generates the Perlin volumes and requires ambient.
#' if (requireNamespace("ambient", quietly = TRUE)) {
#'   # Rounded green hills, with their lower capsule ends buried in the ground.
#'   # The rows are roughly 3-7, 17-26, and 60-85 km from the camera.
#'   # Small foreground hills stay crisp while larger distant hills fade.
#'   hills = data.frame(
#'     x = c(-1.1, 3, -6, -2, 3.5, 8, -27, -17, -6, 7, 22, 35, -45),
#'     z = c(-5.6, -2.2, 9, 15, 11, 17, 55, 63, 69, 58, 66, 60, 62),
#'     radius = c(0.30, 0.72, 1.7, 1.4, 2, 2.2, 6, 5, 5.5, 5, 7, 4, 2),
#'     top = c(0.7, 1.7, 3.8, 3, 4.7, 4.2, 10, 9, 22, 20.5, 31, 30, 34)
#'   )
#'   terrain_mat = diffuse(color = "#469D60")
#'   terrain = xz_rect(xwidth = 160, zwidth = 160, material = terrain_mat)
#'   for (i in seq_len(nrow(hills))) {
#'     h = hills[i, ]
#'     terrain = add_object(
#'       terrain,
#'       csg_object(
#'         csg_capsule(
#'           start = c(h$x, -h$radius, h$z),
#'           end = c(h$x, h$top - h$radius, h$z),
#'           radius = h$radius
#'         ),
#'         material = terrain_mat)
#'      )
#'   }
#'
#'   # A river winds around the capsule footprints and turns out of sight behind
#'   # the distant pair at (-6, 69) and (7, 58). Coordinates and width are in km.
#'   river_bends = data.frame(
#'     x = c(0.3, 0.1, 0.7, 0.8, -1.2, -2.6, -3.9, -4.2, 0.2, 2.5, -1.1, 0.7, 1.4, 1.1, -2, -4.5),
#'     z = c(-12, -7, -4, -1, 3, 7, 12, 17, 23, 32, 43, 53, 62, 69, 76, 79)
#'   )
#'   river_curve = stats::splinefun(river_bends$z, river_bends$x, method = "natural")
#'   river_z = seq(min(river_bends$z), max(river_bends$z), length.out = 600)
#'   river_center = cbind(x = river_curve(river_z), z = river_z)
#'
#'   # Offset perpendicular to the tangent, keeping the river 1 km wide even
#'   # through bends. Reverse the second bank to make one closed polygon.
#'   river_width = 1
#'   river_slope = river_curve(river_z, deriv = 1)
#'   bank_offset = river_width / 2 * cbind(1, -river_slope) / sqrt(1 + river_slope^2)
#'   river_banks = rbind(
#'     river_center + bank_offset,
#'     (river_center - bank_offset)[length(river_z):1, ]
#'   )
#'
#'   # Keep the polygon's world x coordinates and lift its top 1 m above
#'   # ground. A thin extrusion gives the river an upward-facing surface.
#'   terrain = add_object(
#'     terrain,
#'     extruded_polygon(
#'       river_banks,
#'       plane = "xz",
#'       top = 0.001,
#'       bottom = -0.001,
#'       flip_horizontal = TRUE,
#'       material = diffuse("#168BC4")
#'     )
#'   )
#'
#'   # Billowing Perlin volumes sit above each row of hills. Optical depth sets
#'   # the cloud's own scattering; sky_light() separately supplies clear-air haze.
#'   cloud_rows = data.frame(
#'     z = c(3, 21, 63),
#'     base = c(5, 6, 12.5),
#'     width = c(16, 32, 90),
#'     depth = c(10, 16, 24)
#'   )
#'   landscape = terrain
#'   for (i in seq_len(nrow(cloud_rows))) {
#'     cl = cloud_rows[i, ]
#'     landscape = add_object(
#'       landscape,
#'       cloud(
#'         z = cl$z,
#'         y = cl$base + 1.8 / 2,
#'         width = cl$width,
#'         depth = cl$depth,
#'         height = 1.8,
#'         resolution = 64,
#'         coverage = 0.4,
#'         detail = 0.4,
#'         optical_depth = 4,
#'         g = 0.65,
#'         seed = 41 + i
#'       )
#'     )
#'   }
#'
#'   day = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York")
#'   sunset = as.POSIXct("2026-06-21 20:35:00", tz = "America/New_York")
#'
#'   # Reuse camera, seed, and sample count, with manually chosen ISO values.
#'   # Append a white caption strip below the rendered image using rayimage.
#'   render_sky = function(light, scene = landscape, iso = 100, caption = "") {
#'     set.seed(2026)
#'     image = scene |>
#'       add_infinite_light(light) |>
#'       render_scene(
#'         lookfrom = c(0, 0.35, -8),
#'         lookat = c(0, 2.4, 7),
#'         fov = 48,
#'         width = 384,
#'         height = 240,
#'         samples = 32,
#'         integrator_type = "nee",
#'         iso = iso,
#'         tonemap = "raw",
#'         plot_scene = FALSE
#'       )
#'     rayimage::render_stack(list(
#'       image,
#'       rayimage::render_text_image(
#'         caption,
#'         size = 14,
#'         font = "Arial",
#'         width = dim(image)[2],
#'         height = 32,
#'         just = "center",
#'         check_text_width = FALSE,
#'         check_text_height = FALSE
#'       )
#'     ))
#'   }
#'
#'   # Keep ISO fixed within each comparison: 4 by day and 175 at sunset.
#'   day_iso = 4
#'   sunset_iso = 175
#'
#'   # Left to right: cached image sky, native lighting without finite haze,
#'   # native lighting with haze. The image has one altitude for the whole scene.
#'   # All three use Prague; the native modes also query each local altitude.
#'   # A 256-row image can undersample the small Sun disk; the next example
#'   # samples that disk separately instead of relying on the image resolution.
#'   image_sky = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = FALSE,
#'     hosek = FALSE,
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 256
#'   )
#'   modes = list(
#'     image = render_sky(image_sky, iso = day_iso, caption = "Fixed EXR sky"),
#'     native = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         day,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'         attenuation = FALSE,
#'         altitude = 0,
#'         visibility = 50,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = day_iso,
#'       caption = "Atmosphere model, no haze"
#'     ),
#'     haze = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         day,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'         altitude = 0,
#'         visibility = 120,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = day_iso,
#'       caption = "Haze (120 km, p = 0.5)"
#'     )
#'   )
#'   rayimage::plot_image_grid(modes, dim = c(1, 3))
#'
#'   # A detailed, independently sampled Sun can also accompany an image sky.
#'   # Exclude its rasterized disk to avoid adding the Sun twice.
#'   image_with_sun = landscape |>
#'     add_infinite_light(sun_light(
#'       40.7,
#'       -74,
#'       day,
#'       sky_args = list(altitude = 0, visibility = 50, albedo = 0.3)
#'     ))
#'   image_atmosphere = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = FALSE,
#'
#'     hosek = FALSE,
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 256,
#'     render_mode = "atmosphere"
#'   )
#'   rayimage::plot_image_grid(
#'     list(render_sky(
#'       image_atmosphere,
#'       image_with_sun,
#'       iso = day_iso,
#'       caption = "Image sky + separate Sun"
#'     )),
#'     dim = c(1, 1)
#'   )
#'
#'   # The distant capsules lose contrast and approach the sky's color as
#'   # visibility decreases. Visibility is in km, independent of scene units.
#'   visibility = list(
#'     hazy_20_km = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         day,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'
#'         altitude = 0,
#'         visibility = 20,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = day_iso,
#'       caption = "Visibility: 20 km"
#'     ),
#'     normal_50_km = modes$haze,
#'     clear_120_km = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         day,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'
#'         altitude = 0,
#'         visibility = 120,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = day_iso,
#'       caption = "Visibility: 120 km"
#'     )
#'   )
#'   rayimage::plot_image_grid(visibility, dim = c(1, 3))
#'
#'   # At sunset the Sun is about 1.7 degrees below the ground horizon. Clouds
#'   # above 5 km can still see it and receive sunlight on their undersides.
#'   # Turn finite haze off in BOTH first images to isolate query_altitude.
#'   twilight = list(
#'     fixed_altitude = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         sunset,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'         attenuation = FALSE,
#'         query_altitude = FALSE,
#'
#'         altitude = 0,
#'         visibility = 50,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = sunset_iso,
#'       caption = "Fixed altitude, no haze"
#'     ),
#'     cloud_altitude = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         sunset,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'         attenuation = FALSE,
#'         query_altitude = TRUE,
#'
#'         altitude = 0,
#'         visibility = 50,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = sunset_iso,
#'       caption = "Local altitude, no haze"
#'     ),
#'     cloud_altitude_haze = render_sky(
#'       sky_light(
#'         40.7,
#'         -74,
#'         sunset,
#'         atmosphere = TRUE,
#'         meters_per_unit = 1000,
#'
#'         altitude = 0,
#'         visibility = 50,
#'         albedo = 0.3,
#'         resolution = 32
#'       ),
#'       iso = sunset_iso,
#'       caption = "Local altitude + haze"
#'     )
#'   )
#'   rayimage::plot_image_grid(twilight, dim = c(1, 3))
#'
#'   # Both images retain the cloud material. FALSE excludes clear-air haze
#'   # throughout its enclosing boxes, including their empty cells.
#'   rayimage::plot_image_grid(
#'     list(
#'       everywhere = modes$haze,
#'       outside_clouds = render_sky(
#'         sky_light(
#'           40.7,
#'           -74,
#'           day,
#'           atmosphere = TRUE,
#'           meters_per_unit = 1000,
#'           haze_in_volumes = FALSE,
#'
#'           altitude = 0,
#'           visibility = 50,
#'           albedo = 0.3,
#'           resolution = 32
#'         ),
#'         iso = day_iso,
#'         caption = "Haze outside cloud volumes"
#'       )
#'     ),
#'     dim = c(1, 2)
#'   )
#'
#'   # Compare haze settings at equal samples. Deferred haze and the 0.5
#'   # correction probability are the defaults. Their radiance estimates converge
#'   # to eager haze, with different noise before the default denoising step.
#'   # Use more samples for final cloud renders.
#'   rayimage::plot_image_grid(
#'     list(
#'       eager = render_sky(
#'         sky_light(
#'           40.7,
#'           -74,
#'           day,
#'           atmosphere = TRUE,
#'           meters_per_unit = 1000,
#'           deferred_haze = FALSE,
#'
#'           altitude = 0,
#'           visibility = 50,
#'           albedo = 0.3,
#'           resolution = 32
#'         ),
#'         iso = day_iso,
#'         caption = "Eager haze"
#'       ),
#'       deferred_full_correction = render_sky(
#'         sky_light(
#'           40.7,
#'           -74,
#'           day,
#'           atmosphere = TRUE,
#'           meters_per_unit = 1000,
#'           deferred_haze = TRUE,
#'           haze_correction_probability = 1,
#'
#'           altitude = 0,
#'           visibility = 50,
#'           albedo = 0.3,
#'           resolution = 32
#'         ),
#'         iso = day_iso,
#'         caption = "Deferred haze, p = 1"
#'       ),
#'       deferred_half_correction = modes$haze
#'     ),
#'     dim = c(1, 3)
#'   )
#'
#'   # Exact acceleration controls: these change cost and memory, not samples.
#'   # Pass any of these lights to render_sky() using the same seed to compare.
#'   no_cache = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     cache_spectra = FALSE,
#'     transmission_table = FALSE,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'   limited_table = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     transmission_table_max_mb = 64,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'   larger_table = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     transmission_table_max_mb = 1024,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'   unlimited_table = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     transmission_table_max_mb = Inf,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'   no_table_allocation = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     transmission_table_max_mb = 0,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'
#'   # Other controls: rotate the illumination, dim it, and name the light.
#'   rayimage::plot_image_grid(
#'     list(
#'       modes$haze,
#'       render_sky(
#'         sky_light(
#'           40.7,
#'           -74,
#'           day,
#'           atmosphere = TRUE,
#'           meters_per_unit = 1000,
#'           rotation = 35,
#'           intensity = 0.7,
#'           name = "rotated_sky",
#'
#'           altitude = 0,
#'           visibility = 50,
#'           albedo = 0.3,
#'           resolution = 32
#'         ),
#'         iso = day_iso,
#'         caption = "Rotation: 35 deg, intensity: 0.7"
#'       )
#'     ),
#'     dim = c(1, 2)
#'   )
#'   # These two reference frames describe the same physical atmosphere:
#'   # sea level at world y=0, or 1000 m altitude at world y=1 (units are km).
#'   raised_reference = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'     atmosphere_origin = c(0, 1, 0),
#'
#'     altitude = 1000,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32
#'   )
#'
#'   # Select sky only or Sun only with render_mode; retain all for normal use.
#'   sky_only = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32,
#'     render_mode = "atmosphere"
#'   )
#'   sun_only = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32,
#'     render_mode = "sun"
#'   )
#'   # Disable or adjust Prague's RGB color correction.
#'   uncorrected = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32,
#'     prague_rgb_correction = FALSE
#'   )
#'   custom_color = sky_light(
#'     40.7,
#'     -74,
#'     day,
#'     atmosphere = TRUE,
#'     meters_per_unit = 1000,
#'
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3,
#'     resolution = 32,
#'     prague_rgb_correction = TRUE,
#'     prague_rgb_correction_strength = 0.5,
#'     prague_rgb_correction_gain = c(1, 0.95, 0.9)
#'   )
#' }
sky_light = function(
  lat,
  long,
  datetime,
  intensity = 1,
  rotation = 0,
  name = "sky",
  atmosphere = FALSE,
  meters_per_unit = 1,
  atmosphere_origin = c(0, 0, 0),
  attenuation = TRUE,
  query_altitude = TRUE,
  haze_in_volumes = TRUE,
  deferred_haze = TRUE,
  haze_correction_probability = if (isTRUE(deferred_haze)) 0.5 else 1,
  cache_spectra = TRUE,
  transmission_table = TRUE,
  transmission_table_max_mb = 512,
  altitude = 0,
  visibility = 50,
  albedo = 0.5,
  resolution = if (isTRUE(atmosphere)) 64 else 2048,
  hosek = !isTRUE(atmosphere),
  render_mode = "all",
  turbidity = 3,
  wide_spectrum = FALSE,
  below_horizon = TRUE,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto",
  stars = FALSE,
  star_width = 1,
  stars_exposure = 0,
  planets = FALSE,
  moon = FALSE,
  moon_atmosphere = FALSE,
  moon_hosek = TRUE,
  exr_adopted_white = "D60",
  exr_metadata = TRUE,
  number_cores = 1,
  verbose = FALSE,
  ...
) {
  extra = list(...)
  if ("sky_args" %in% names(extra)) {
    stop(
      "Pass sky settings directly to sky_light(), not in sky_args.",
      call. = FALSE
    )
  }

  # These settings are shared by the cached image and native Prague paths.
  sky_args = list(
    altitude = altitude,
    visibility = visibility,
    albedo = albedo,
    resolution = resolution,
    hosek = hosek,
    render_mode = render_mode,
    wide_spectrum = wide_spectrum,
    below_horizon = below_horizon,
    prague_rgb_correction = prague_rgb_correction,
    prague_rgb_correction_strength = prague_rgb_correction_strength,
    prague_rgb_correction_gain = prague_rgb_correction_gain,
    stars = stars,
    planets = planets,
    moon = moon,
    number_cores = number_cores,
    verbose = verbose
  )

  # Native atmosphere rejects explicitly requested image-only controls. Keep
  # their defaults out of the native description so they do not change its API.
  image_args = list(
    turbidity = turbidity,
    star_width = star_width,
    stars_exposure = stars_exposure,
    moon_atmosphere = moon_atmosphere,
    moon_hosek = moon_hosek,
    exr_adopted_white = exr_adopted_white,
    exr_metadata = exr_metadata
  )
  if (isTRUE(atmosphere)) {
    supplied = names(match.call(expand.dots = FALSE))
    image_args = image_args[intersect(names(image_args), supplied)]
  }
  sky_args = c(sky_args, image_args, extra)

  result = structure(
    list(
      type = "sky",
      lat = lat,
      long = long,
      datetime = datetime,
      sky_args = sky_args,
      intensity = intensity,
      rotation = rotation,
      name = name,
      atmosphere = atmosphere,
      meters_per_unit = meters_per_unit,
      atmosphere_origin = atmosphere_origin,
      attenuation = attenuation,
      query_altitude = query_altitude,
      haze_in_volumes = haze_in_volumes,
      deferred_haze = deferred_haze,
      haze_correction_probability = haze_correction_probability,
      cache_spectra = cache_spectra,
      transmission_table = transmission_table,
      transmission_table_max_mb = transmission_table_max_mb
    ),
    class = "ray_infinite_light"
  )
  validate_infinite_light(result)
  result
}

#' @keywords internal
validate_sky_light = function(light) {
  # Descriptions serialized before the atmosphere option retain the image path.
  if ("atmosphere" %in% names(light)) {
    if (
      !is.logical(light$atmosphere) ||
        length(light$atmosphere) != 1L ||
        is.na(light$atmosphere)
    ) {
      stop("atmosphere must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (field in c(
    "attenuation",
    "query_altitude",
    "haze_in_volumes",
    "deferred_haze",
    "cache_spectra",
    "transmission_table"
  )) {
    if (
      field %in%
        names(light) &&
        (!is.logical(light[[field]]) ||
          length(light[[field]]) != 1L ||
          is.na(light[[field]]))
    ) {
      stop(field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if ("haze_correction_probability" %in% names(light)) {
    probability = light$haze_correction_probability
    if (
      !is.numeric(probability) ||
        length(probability) != 1L ||
        !is.finite(probability) ||
        probability <= 0 ||
        probability > 1
    ) {
      stop(
        "haze_correction_probability must be a finite number greater than zero and at most one.",
        call. = FALSE
      )
    }
  }
  if ("transmission_table_max_mb" %in% names(light)) {
    limit = light$transmission_table_max_mb
    if (
      !is.numeric(limit) || length(limit) != 1L || is.na(limit) || limit < 0
    ) {
      stop(
        "transmission_table_max_mb must be a nonnegative number or Inf.",
        call. = FALSE
      )
    }
  }
  for (field in c("lat", "long")) {
    value = light[[field]]
    limit = if (field == "lat") 90 else 180
    if (
      !is.numeric(value) ||
        length(value) != 1 ||
        !is.finite(value) ||
        abs(value) > limit
    ) {
      stop(
        "Sky light ",
        field,
        " must be a finite number between -",
        limit,
        " and ",
        limit,
        ".",
        call. = FALSE
      )
    }
  }
  if (
    !inherits(light$datetime, "POSIXct") ||
      length(light$datetime) != 1 ||
      !is.finite(as.numeric(light$datetime))
  ) {
    stop(
      "Sky light datetime must be one finite POSIXct date and time.",
      call. = FALSE
    )
  }
  args = light$sky_args
  if (
    !is.list(args) ||
      (length(args) &&
        (is.null(names(args)) ||
          anyNA(names(args)) ||
          any(!nzchar(names(args))) ||
          anyDuplicated(names(args))))
  ) {
    stop(
      if (light$type == "sky") {
        "Sky arguments must be uniquely named."
      } else {
        "sky_args must be a uniquely named list."
      },
      call. = FALSE
    )
  }
  reserved = intersect(
    names(args),
    c("lat", "lon", "long", "datetime", "filename", "allow_download")
  )
  if (length(reserved)) {
    stop(
      if (light$type == "sky") {
        "Sky arguments cannot override: "
      } else {
        "sky_args cannot override: "
      },
      paste(reserved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (isTRUE(light$atmosphere)) {
    validate_prague_sky_light(light)
  }
  invisible(TRUE)
}

#' @keywords internal
prepare_infinite_light = function(light) {
  validate_infinite_light(light)
  if (isTRUE(light$atmosphere)) {
    return(prepare_prague_sky_light(light))
  }
  if (light$type %in% c("image", "disk")) {
    return(light)
  }
  if (light$type %in% c("sun", "moon")) {
    return(prepare_celestial_light(light))
  }
  if (!requireNamespace("skymodelr", quietly = TRUE)) {
    stop(
      "sky_light() requires skymodelr. Install it with install.packages('skymodelr').",
      call. = FALSE
    )
  }
  args = c(
    list(lat = light$lat, lon = light$long, datetime = light$datetime),
    light$sky_args
  )
  supported = names(formals(skymodelr::generate_sky_latlong))
  unknown = setdiff(names(args), supported)
  if (length(unknown) && !"..." %in% supported) {
    stop(
      "Unsupported sky arguments: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if ("allow_download" %in% supported) {
    args$allow_download = FALSE
  }
  cache = file.path(tempdir(), "rayrender-skies")
  dir.create(cache, showWarnings = FALSE)
  key_file = tempfile(tmpdir = cache)
  on.exit(unlink(key_file), add = TRUE)
  saveRDS(
    list(
      args = args,
      version = as.character(utils::packageVersion("skymodelr"))
    ),
    key_file,
    version = 2
  )
  filename = file.path(cache, paste0(unname(tools::md5sum(key_file)), ".exr"))
  if (!file.exists(filename)) {
    success = FALSE
    on.exit(if (!success) unlink(filename), add = TRUE)
    args$filename = filename
    do.call(skymodelr::generate_sky_latlong, args)
    success = TRUE
  }
  infinite_light(
    filename,
    intensity = light$intensity,
    rotation = light$rotation,
    name = light$name
  )
}
