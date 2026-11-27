"""Render the authored Blender cycles as PNG sequences for movie encoding."""
import argparse
from pathlib import Path
import sys
import bpy

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--only',choices=['bird','jellyfish','walk'])
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
for name in ['bird','jellyfish','walk']:
    if args.only and args.only!=name:continue
    directory=args.output/name
    bpy.ops.wm.open_mainfile(filepath=str(directory/(name+'.blend')))
    scene=bpy.context.scene
    preferences=bpy.context.preferences.addons['cycles'].preferences
    if sys.platform=='darwin':
        preferences.compute_device_type='METAL'
        preferences.get_devices()
        devices=[device for device in preferences.devices if device.type=='METAL']
        if devices:
            for device in preferences.devices:device.use=device.type=='METAL'
            scene.cycles.device='GPU'
            print('RENDER_DEVICE',[(device.name,device.type) for device in devices],flush=True)
    scene.cycles.samples=24
    scene.render.resolution_x=640
    scene.render.resolution_y=640
    (directory/'frames').mkdir(exist_ok=True)
    scene.render.filepath=str(directory/'frames'/'frame-')
    scene.frame_start=1;scene.frame_end=48
    bpy.ops.render.render(animation=True)
