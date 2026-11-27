These tiny fixtures are original rayrender test data, under the package license.
No external model or network resource is required.

rigid.gltf, rigid-external.gltf plus rigid.bin, and rigid.glb describe the same
two-part tetrahedral model. It uses interleaved positions/normals/normalized UVs,
two sparse position replacements, three named clips, asynchronous channels,
STEP/linear/cubic interpolation, duplicated node/clip labels, and a second scene
with a reflected matrix node. Cubic quaternion signs are deliberately authored
to differ from shortest-arc SLERP.

The invalid-* and network fixtures must fail. skin.gltf is intentionally invalid:
its skin has no paired JOINTS/WEIGHTS attributes. morph.gltf is a small static
morph example.

The deform-skin, deform-morph, and deform-combined glTF/GLB pairs are analytic
triangles with a nonidentity mesh node, two joints, two weight sets, default
morph weights, asynchronous channels, and additional cubic/STEP weight clips.
The deform-bend pair changes surface normals and triangle area. Tests use
independent numerical positions and explicit static mesh render references.
Regenerate these original fixtures without Blender or downloads:

    python3 inst/animation/fixtures/create-deformation-fixtures.py

Regenerate in a checkout with:

    python3 dev/animation/generate-gltf-fixtures.py
