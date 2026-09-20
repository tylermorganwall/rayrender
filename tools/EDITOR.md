# Native editor workspace

Requires rayimgui 0.0.12 (ABI 1.8). Run `source("tools/editor-demo.R")` from this
rayrender checkout after installing both packages.

- **Left:** render progress, sample count and FPS, exposure, denoising,
  fast-preview and camera controls, plus
  sky controls for scenes with a sky light.
- **Center:** the rendered viewport and selected-object transform gizmo.
- **Right, upper panel:** expandable groups, nested groups, objects and instance
  placements. Meshes remain single selectable objects.
- **Right, lower panel:** selected-object position, rotation and scale, followed
  by the material slot selector and that material's editable properties.

**Escape** closes the editor from any panel, including while editing a text
field, dragging a transform, or playing an animation. Rendering workers finish
safely before the window and its resources are released.

**Undo / Redo** buttons on the left share a session history with
**Control/Command-Z** and **Control/Command-Shift-Z**. Shortcuts work from any
editor panel. An active text field keeps its own text undo; finish text entry
before using scene history. Each drag is one history step. History includes
input drafts and applied object/material/texture edits, camera and render settings,
sky settings, and saved keyframes. Navigation-only selection and animation playback
do not fill the history. Undo during playback stops it before restoring the edit.
A new edit after undo clears redo. The latest 100 edits are kept in memory.
Failed scene/sky restoration leaves the current scene and history position intact
and reports the error beside Undo/Redo.

The **bottom Animation panel** saves camera keyframes with **Save keyframe** (K).
Use **Previous/Next** (Shift-comma/Shift-period), click a snapshot, or **Delete** (/)
to navigate/edit the path. The selected snapshot has a highlighted border. Thumbnails
capture a completed preview at the saved camera, without the selection outline or
gizmo, and stay fixed as rendering continues. If that camera has not rendered yet,
the thumbnail says pending until it does. These are camera keyframes, not snapshots
of object/material edits. They are kept for the current editor session.

**Play/Stop** (M) previews the path and restores the starting camera when finished
or stopped. Save at least two views. **Open path / Closed loop** (Shift-L) controls
whether the last keyframe connects back to the first; it does not repeat playback
forever. **Camera motion blur** (B) and **Shutter (frame fraction)** control camera blur
during movement: 0 freezes motion and 1 exposes the full interval between frames.
Path editing is disabled during playback. Existing
`keyframe_motion_args` settings still determine interpolation and frame count.
The snapshot row scrolls horizontally as more keyframes are saved. The native editor
supports up to 120 saved keyframes; deleting one frees its thumbnail texture.

The default scene is an outdoor sculpture display: a large `r_obj()` logo on
a stepped plinth, three instanced R letters, and metal, glass and glossy surface
examples. Nested groups and shared instances exercise the scene hierarchy.
A 55-degree field of view and low camera angle keep the sky and horizon visible.
The gold metal sphere uses `editor-assets/watermask.00000.png` as a roughness
texture, mapping dark/light pixels to roughness 0.03/0.45. The demo resolves the
asset relative to its script, so it can be sourced from another working directory.

The demo starts with a **Hosek** sky. The left **Sky model** selector also
offers **Prague**. Sun elevation uses a quarter-circle handle (horizon to
zenith); azimuth uses a clockwise circular dial (north 0, east 90). The
compact controls sit side by side, each with an editable numeric value
inside. Negative elevations can be typed for a Sun below the horizon; the
elevation handle then stays at the horizon endpoint. Sun changes apply when
the handle or numeric input is released. **Location and date/time** is
collapsed by default; expand it to edit latitude, longitude and UTC
date/time, then use **Apply location/time**. Applying a date restores the
astronomical Sun direction. Model changes preserve the committed
location/time and any manual Sun direction.

The editor's **Prague** selection uses native atmospheric lighting and
exposes **Haze** and **Altitude queries**, including when the scene starts
with `sky_light_image()`. These settings survive switches to Hosek and back.
Haze requires altitude queries: enabling haze enables both; disabling
altitude queries disables haze. Hosek uses image-based lighting. The
interactive sky editor uses NEE transport so an atmosphere can be enabled
during the session; ordinary image-sky renders retain their chosen
integrator. Prague requires the installed full-altitude skymodelr sky data.
Load failures appear in the panel and preserve the current sky. Existing
explicit sky choices are respected; the demo and `sky_light_image()` default
to Hosek.

The FPS counter beside the sample count measures completed preview frames,
including denoising time, averaged over at least half a second. GUI redraws and
selection-overlay updates do not count as rendered frames.

The **Denoise** checkbox follows the initial `denoise` argument and controls both
normal/fast preview and the final render. It is disabled when denoising support is
unavailable. Toggling it preserves accumulated samples and the selection mask.

Panels start docked. Drag splitters to resize them, or drag tabs to rearrange them.
Click a tree label to select; its arrow only expands/collapses. **Shift-click**
the viewport to select the outer object, group, or instance placement. Another
Shift-click inside the selection descends one hierarchy level toward the hit
object. Further Shift-clicks continue through nested groups and instances; at
a leaf, selection stays there. Clicking a sibling with Shift selects that sibling
at the current depth. Shift-clicking empty space clears selection.

Plain left-click (or Alt-left-click) sets camera target and focus; right-click
sets the target. Holding Shift gives selection priority over idle gizmo handles;
an already-active gizmo drag keeps ownership until release. Movement shortcuts
continue to respect text editing, focus, and gizmo ownership.

Selecting a child in the tree permits editing that child independently of its
group. Instance placements also expose their source objects and nested groups
in the tree. Child transforms and materials affect only that placement's copy;
untouched copies retain shared geometry. Selecting a group or an instance set
transforms all its descendants.
Selecting an **Instances** parent exposes one panel per source material slot;
Apply material changes that slot across every placement. Selecting an individual
placement still edits only that placement. Mixed properties are marked in the
inspector, and unchanged properties retain their individual values.
Material slots retain their type. The inspector exposes their shading inputs:

| Material | Surface controls |
| --- | --- |
| Diffuse / Oren-Nayar | Color and sigma in degrees, including increasing sigma from zero |
| Metal | Color, fuzz, eta RGB and kappa RGB |
| Dielectric | Tint, refraction, RGB absorption and overlap priority |
| Microfacet reflection | Color, anisotropic roughness, distribution, eta and kappa |
| Microfacet transmission | Color, anisotropic roughness, distribution and refraction |
| Glossy | Color, anisotropic gloss, distribution and specular reflectance |
| Light / spotlight | Color, intensity, visibility; spotlight direction, width and falloff |
| Hair | Absorption, color or pigment mode; refraction, longitudinal/azimuthal roughness and scale angle |

Color-bearing materials have a **Texture mode** selector for solid color,
checkers, noise, UV/world gradients, or an image. The original imported texture
can also be retained and tinted. Pattern controls appear for the selected mode.
Choose **Image** and enter a **Color texture file** path to assign a texture;
**Image repeat U/V** controls tiling. File paths may be absolute, relative to the
R working directory, or begin with `~`.

Enable **Use alpha map**, **Use bump map**, or **Use roughness map** and enter the
corresponding file path. Alpha and bump controls appear on geometry that supports
them; roughness maps belong to microfacet/glossy materials. Bump controls include
intensity and tiling; roughness maps include range and flip controls. A blank path
uses the original/embedded map when one exists. Disabling a map preserves its
settings. Grayscale/RGB alpha images describe opacity; RGBA images use their alpha
channel. Bump images must be at least 3 by 3 pixels.

Click **Apply material** to commit the edited inputs together. Invalid values or
unreadable files leave the live scene unchanged. Replacement textures own their
pixel storage, and editing a child/placement does not change another instance.
Mesh opacity caches are refreshed when alpha maps change, including for shadows.
The exported `scene_edits` values include texture paths as strings. Surface type,
volume construction and importance-sampling configuration remain defined by the R
scene; the inspector edits the selected material and its applicable texture maps.

Selection is shown only by a fixed dark outer outline with a light inner band.
The interior is transparent, preserving the object's rendered colors. Both bands
come from the binary mask using all eight neighboring pixels, including diagonals;
render noise cannot change their coverage or color.

Foreground objects occlude the mask. Coverage and its outline are cached until
the camera, scene, selection, or mask dimensions change. Restarting accumulation,
toggling fast preview, and changing lighting/exposure reuse the cached outline.
Selecting a group includes its descendants. ImGui draws the outline as a separate
layer; snapshots, saved images, and materials keep their original colors.

The coverage mask uses centered lens/midpoint shutter rays, like picking, and is
capped at 768 pixels on its longest edge. It indicates geometric selection rather
than depth-of-field, motion-blurred or refracted silhouettes.
A new selection discards unapplied edits. Numeric transforms/material edits use
Apply; releasing a gizmo commits its transform. Scene rebuilds happen after render
workers finish their current sample, and errors preserve the previous live scene.
Returned `scene_edits` records place descendant overrides under `children` on
their enclosing instance entry. Each child row/instance index addresses that
instance's source scene, and its transform delta uses those source coordinates.
The inspector and gizmo continue to show world coordinates.

Downstream editor logic lives here. Rayimgui supplies only generic docking slots,
widgets, viewport input, and gizmos through the public C protocol.

## Export the current scene

The left panel has an **Export file** field and **Export R code** button. The
button writes the committed object/material edits, camera, sky, environment
rotation, denoising and calibrated viewport exposure. Inspector fields that still
need **Apply** are excluded. Exports use the requested full render size and sample
budget; fast-preview resolution, selection outlines and gizmos are display aids.

Each export contains a runnable `.R` script and a companion `-scene.rds` file
holding the original geometry, materials and instance sources. The script lists
cumulative object edits and sky settings as editable R values, applies them with
`apply_scene_edits()`, then calls `render_scene()` with the current camera and
appearance settings. It needs no variables from the original R workspace.
Reopening and exporting an edited scene applies each change once.

Keep these files together. Temporary textures and OBJ/MTL dependencies are copied
to a companion `-assets` directory; permanent model and texture files retain
absolute paths. Moving the export to another machine also requires those external
files and any Prague sky dataset. Existing exports are never overwritten: repeated
clicks produce numbered filenames. Write errors appear in the panel and leave the
renderer running.

The exported render freezes the viewport exposure and omits final-only bloom.
It recreates the scene and view; sampling noise can differ from the progressive
preview. The returned image retains `scene_edits` and additionally provides
`editor_state` with committed camera, sky and appearance values.
