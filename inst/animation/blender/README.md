# Blender animation assets

These scripts build actual three-dimensional models and editable cycles, rather
than using the tiny triangle fixtures as visual demonstrations. Blender 5.2 was
used. Generated files are kept outside the package because the walk uses the
user-supplied `rayman.blend`; that source file is never overwritten.

## Assets

- **bird**: a rounded kingfisher-inspired bird with a beak, eyes, breast, tail,
  layered flight feathers, and separate shoulder/wrist motion. Rigid node clip.
- **jellyfish**: a thick scalloped bell, radial canals, four oral arms, and sixteen
  trailing tentacles. Bell contraction and delayed traveling waves use morph
  targets; all surfaces have thickness and explicit smooth normals.
- **walk**: the male mesh and Rigify rig from `rayman.blend`, with a new in-place
  walk. The cycle has stance/swing phases, evaluated sole contact, alternating
  arm swings, and torso/hip counter-motion. Knees bend in a forward plane with
  one active pole solver per leg, ankles track beneath the hips, and the pelvis
  rises through passing to the up pose. The longer step places the leading heel
  about 0.40 model units ahead of the hip, with raised toes that roll down during
  loading. The shoulders stay upright over the pelvis, with a level gaze.
  Vertical bob is reduced to 40 percent of the earlier cycle, while shoulder
  swing has .50-radian amplitude and a -.25-radian forward bias, giving more
  forward reach and less rearward extension. Hip and chest pose controls soften the lumbar bend while
  preserving the pelvis trajectory and gaze; the waist joint angle remains
  below 16 degrees through the cycle. Contact/down/passing/up markers identify
  the key poses.
  Stray thigh weights on the abdomen are redistributed to its existing spine
  and pelvis weights, with a smooth transition into the groin and flank. This
  removes the navel-side divot during leg swing without moving rest vertices.
  A separate right-hip weight correction places flexion along the diagonal
  pelvic line, with a broad pelvis/thigh transition and a smooth blend through
  the central groin. It removes the horizontal upper-thigh pinch while keeping
  the left outer hip, rest vertices, and skeleton animation unchanged.
  Calf weights use a smooth thigh/shin transition to prevent thigh influence
  from pulling lower-leg vertices into folds. Quad diagonals are selected over
  the cycle before export, avoiding the worst thin-triangle inversions while
  retaining the authored smooth rest normals in the GLB. The rayrender walk
  example uses `read_gltf(..., recalculate_normals = TRUE)` to derive smooth
  shading from the deformed triangles at every requested time, including the
  toes. This replaces authored normals explicitly; other imports keep their
  default normal behavior. All remaining influences export.
  The generated copy uses ordinary linear skinning with single-segment bones
  (glTF does not represent Blender B-Bone interpolation). The original source
  retains its segmented B-Bones.
  The source's saved rig lacks IK/FK influence drivers, so its existing arm FK
  constraints are explicitly enabled in the generated copy. Exported joints are
  flattened to preserve the rig's inheritance rules; the editable control
  hierarchy remains intact.

All clips are two seconds at 24 fps. Frame 49 equals frame 1; playback uses
frames 1 through 48, without duplicating the seam. Each asset has an independent
outer placement object. Camera, lighting, and ground belong to the Blender
preview scene and are excluded from the GLB.

Each output subdirectory contains an editable `.blend`, a `.glb`, four pose
renders, and a PNG frame sequence when movie rendering is requested. The movies
show three repetitions of the same cycle. The character is a smooth mannequin
using the supplied geometry; no facial features or clothing are invented.

## Rebuild

Run from the rayrender checkout. The output directory can be changed. Existing
generated outputs at that destination will be replaced; use a fresh directory
if you have edited them in Blender.

```sh
/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/create_cycles.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles" \
  --rayman "$HOME/Desktop/3dmodels/rayman.blend"

/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/validate_cycles.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"

/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/render_cycles.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"
```

To check the exported surfaces after reimport, assemble the review document,
and encode the completed frame sequences, run:

```sh
/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/check_roundtrip.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"
/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/assemble_showcase.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"
python3 inst/animation/blender/encode_cycles.py \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"
```

Open `animation-cycles.blend` and choose a cycle from the scene selector. Press
Space over the viewport to play. `animation-cycles.mp4` shows all three previews
side by side; the individual subdirectories also contain their own movies.

For a front/side review of contact, loading, passing, and recovery poses:

```sh
/Applications/Blender.app/Contents/MacOS/Blender --background --python-exit-code 1 \
  --python inst/animation/blender/review_walk.py -- \
  --output "$HOME/Desktop/3dmodels/rayrender-cycles"
```

The eight images are saved under `walk/review/`.

The create/render scripts accept `--only bird`, `--only jellyfish`, or
`--only walk`. `create_cycles.py --no-render` skips the four pose renders.
Movie rendering uses Blender Cycles with Metal devices when available on macOS,
and otherwise the saved scene's CPU renderer. These are Blender previews, not
rayrender performance or rendering validation.

Encode a repeating movie from each frame sequence with:

```sh
for asset in bird jellyfish walk; do
  directory="$HOME/Desktop/3dmodels/rayrender-cycles/$asset"
  ffmpeg -y -framerate 24 -start_number 1 -i "$directory/frames/frame-%04d.png" \
    -frames:v 48 -c:v libx264 -crf 18 -pix_fmt yuv420p "$directory/cycle.mp4"
  ffmpeg -y -stream_loop 2 -i "$directory/cycle.mp4" -c copy \
    -movflags +faststart "$directory/$asset.mp4"
done
```

## Validation and current rayrender boundary

`validate_cycles.py` checks saved-file evaluation at multiple poses and reverse
seek orders, exact seam closure, finite exported positions, unit normals,
triangle indices, clip duration, and the expected skin/morph data. Results are
written to `validation.json`. It also checks sole contact on every stance frame
and verifies hand swing, forward knee alignment, narrow ankle lanes, toe
direction, actual heel/toe contact, extended steps, upright torso, reduced bob, and
pelvis height ordering, and the corrected waist bend. A geometric dihedral check catches recurrence
of the horizontal pelvic fold independently of shading normals. `check_roundtrip.py` compares deformed vertex surfaces
at five poses (nine for the walk) after GLB reimport, independently of vertex indexing, and writes
`roundtrip-validation.json`. Pose images and movies also require visual review.

All three GLBs can now be imported by `read_gltf()`: the bird uses rigid nodes,
the jellyfish uses morph channels, and the walker uses skeletal animation.
See [the native rayrender examples](../README.md) to render them on independent
placement paths or scrub them through the preview timeline. Deforming models
require explicit scene times and zero-width exposure. The importer replaces
glTF materials with explicit rayrender materials; the example assigns part
colors and rayrender studio lighting. Blender previews and rayrender movies
are saved separately.

The exporter uses the official Blender
[glTF export API](https://docs.blender.org/api/5.2/bpy.ops.export_scene.html).
Bird/jellyfish geometry is generated by this script. Rights to the supplied
character remain with its original owner; its binary data is not bundled here.
