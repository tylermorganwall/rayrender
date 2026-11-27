"""Render front and side views for checking knee alignment and stride poses."""
import argparse
from pathlib import Path
import sys
import bpy
from mathutils import Vector

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
directory=args.output/'walk'
bpy.ops.wm.open_mainfile(filepath=str(directory/'walk.blend'))
scene=bpy.context.scene
scene.cycles.samples=16
scene.render.resolution_x=480
scene.render.resolution_y=640
scene.camera.data.ortho_scale=2.50
(directory/'review').mkdir(exist_ok=True)
for view,position in [('front',(0,-6,1.05)),('side',(6,0,1.05))]:
    scene.camera.location=position
    scene.camera.rotation_euler=(Vector((0,0,1.05))-scene.camera.location).to_track_quat('-Z','Y').to_euler()
    for frame in [1,7,13,19]:
        scene.frame_set(frame)
        scene.render.filepath=str(directory/'review'/('%s-%02d.png'%(view,frame)))
        bpy.ops.render.render(write_still=True)
