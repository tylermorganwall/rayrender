"""Collect the three editable scenes into one Blender review document."""
import argparse
from pathlib import Path
import sys
import bpy

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
bpy.ops.wm.read_factory_settings(use_empty=True)
initial=bpy.context.scene
scenes=[]
for index,name in enumerate(['bird','jellyfish','walk'],1):
    path=args.output/name/(name+'.blend')
    with bpy.data.libraries.load(str(path),link=False) as (source,target):
        target.scenes=list(source.scenes)
    scene=target.scenes[0]
    scene.name='%02d %s cycle'%(index,name.title())
    scenes.append(scene)
bpy.context.window.scene=scenes[0]
bpy.data.scenes.remove(initial)
for scene in scenes:scene.frame_set(1)
for screen in bpy.data.screens:
    for area in screen.areas:
        if area.type=='VIEW_3D':
            area.spaces.active.region_3d.view_perspective='CAMERA'
            area.spaces.active.overlay.show_overlays=False
            area.spaces.active.shading.type='MATERIAL'
bpy.ops.wm.save_as_mainfile(filepath=str(args.output/'animation-cycles.blend'))
print('SHOWCASE_COMPLETE',flush=True)
