"""Compare exported GLB surface positions against the saved Blender source."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

# Import the adjacent shared evaluator when run through Blender's --python.
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from validate_cycles import evaluated_vertices

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
results={}
for name,root_name in [('bird','Bird_Placement'),('jellyfish','Jellyfish_Placement'),('walk','Rayman_Placement')]:
    directory=args.output/name
    bpy.ops.wm.open_mainfile(filepath=str(directory/(name+'.blend')))
    source={}
    frames=list(range(1,50,6)) if name=='walk' else [1,13,25,37,49]
    for frame in frames:
        bpy.context.scene.frame_set(frame)
        source[frame]=evaluated_vertices(bpy.data.objects[root_name])
    bpy.ops.wm.read_factory_settings(use_empty=True)
    bpy.context.scene.render.fps=24
    bpy.ops.import_scene.gltf(filepath=str(directory/(name+'.glb')))
    maximum=0
    # glTF splits vertices at normal/UV seams. Compare the surface's vertex
    # positions independent of index ordering, using bounded memory chunks.
    for frame in frames:
        bpy.context.scene.frame_set(frame)
        actual=evaluated_vertices(bpy.data.objects[root_name])
        frame_error=0
        for key,expected in source[frame].items():
            points=actual[key]
            for offset in range(0,len(expected),64):
                chunk=expected[offset:offset+64]
                distances=np.sum((chunk[:,None,:]-points[None,:,:])**2,axis=2)
                error=float(np.sqrt(distances.min(axis=1).max()))
                if error>frame_error:
                    frame_error=error
                    point=chunk[np.argmax(distances.min(axis=1))]
                maximum=max(maximum,error)
        if name=='walk':print('WALK_SURFACE',frame,frame_error,point.tolist(),flush=True)
    results[name]={'maximum_surface_vertex_error':maximum,'frames_checked':frames}
    print('ROUNDTRIP',name,maximum,flush=True)
    assert maximum<.0001,(name,maximum)
(args.output/'roundtrip-validation.json').write_text(json.dumps(results,indent=2)+'\n')
