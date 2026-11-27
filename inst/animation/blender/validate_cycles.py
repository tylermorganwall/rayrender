"""Check the saved Blender assets and their exported GLB geometry/animation."""
import argparse
import json
import math
from pathlib import Path
import struct
import sys

import bpy
import numpy as np


def evaluated_vertices(root):
    graph=bpy.context.evaluated_depsgraph_get()
    result={}
    for obj in root.children_recursive:
        if obj.type!='MESH':continue
        evaluated=obj.evaluated_get(graph)
        mesh=evaluated.to_mesh()
        coords=np.empty(len(mesh.vertices)*3,dtype=np.float64)
        mesh.vertices.foreach_get('co',coords)
        coords=coords.reshape(-1,3)
        matrix=np.array(evaluated.matrix_world)
        result[obj.name]=coords @ matrix[:3,:3].T+matrix[:3,3]
        evaluated.to_mesh_clear()
    return result


def glb_report(path):
    data=path.read_bytes()
    magic,version,length=struct.unpack_from('<III',data)
    assert magic==0x46546c67 and version==2 and length==len(data)
    json_length,json_type=struct.unpack_from('<II',data,12)
    assert json_type==0x4e4f534a
    doc=json.loads(data[20:20+json_length])
    offset=20+json_length
    bin_length,bin_type=struct.unpack_from('<II',data,offset)
    assert bin_type==0x004e4942
    binary=data[offset+8:offset+8+bin_length]
    def accessor(index):
        a=doc['accessors'][index];view=doc['bufferViews'][a['bufferView']]
        dtype={5126:'<f4',5125:'<u4',5123:'<u2',5121:'u1'}[a['componentType']]
        width={'SCALAR':1,'VEC2':2,'VEC3':3,'VEC4':4,'MAT4':16}[a['type']]
        itemsize=np.dtype(dtype).itemsize
        return np.ndarray((a['count'],width),dtype=dtype,buffer=binary,
            offset=view.get('byteOffset',0)+a.get('byteOffset',0),
            strides=(view.get('byteStride',width*itemsize),itemsize)).copy()
    vertices=triangles=targets=0
    normal_error=0
    for mesh in doc['meshes']:
        for primitive in mesh['primitives']:
            assert primitive.get('mode',4)==4
            positions=accessor(primitive['attributes']['POSITION'])
            normals=accessor(primitive['attributes']['NORMAL'])
            assert np.isfinite(positions).all() and np.isfinite(normals).all()
            normal_error=max(normal_error,float(np.max(np.abs(np.linalg.norm(normals,axis=1)-1))))
            indices=accessor(primitive['indices']).ravel()
            assert len(indices)%3==0 and indices.max()<len(positions)
            vertices+=len(positions);triangles+=len(indices)//3
            targets+=len(primitive.get('targets',[]))
    assert normal_error<2e-5
    assert len(doc.get('animations',[]))==1
    duration=0
    for sampler in doc['animations'][0]['samplers']:
        times=accessor(sampler['input']).ravel()
        assert np.all(np.diff(times)>0)
        duration=max(duration,float(times[-1]-times[0]))
    assert abs(duration-2)<1e-6
    return {'meshes':len(doc['meshes']),'vertices':vertices,'triangles':triangles,
        'skins':len(doc.get('skins',[])),'morph_targets':targets,
        'clips':len(doc['animations']),'duration_seconds':duration,
        'maximum_normal_length_error':normal_error}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    options=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    results={}
    for name,root_name in [('bird','Bird_Placement'),('jellyfish','Jellyfish_Placement'),('walk','Rayman_Placement')]:
        directory=options.output/name
        bpy.ops.wm.open_mainfile(filepath=str(directory/(name+'.blend')))
        scene=bpy.context.scene;root=bpy.data.objects[root_name]
        snapshots={}
        for frame in [1,13,25,37,49]:
            scene.frame_set(frame);snapshots[frame]=evaluated_vertices(root)
        seam=max(float(np.max(np.abs(snapshots[1][key]-snapshots[49][key]))) for key in snapshots[1])
        assert seam<1e-5,(name,seam)
        motion=max(float(np.max(np.abs(snapshots[1][key]-snapshots[25][key]))) for key in snapshots[1])
        assert motion>.05,(name,motion)
        seek_error=0
        for frame in [37,1,25,13,49,1]:
            scene.frame_set(frame);actual=evaluated_vertices(root)
            seek_error=max(seek_error,max(float(np.max(np.abs(actual[key]-snapshots[frame][key]))) for key in actual))
        assert seek_error<1e-7,(name,seek_error)
        report=glb_report(directory/(name+'.glb'))
        report.update(loop_seam_error=seam,reverse_seek_error=seek_error,half_cycle_displacement=motion)
        if name=='walk':
            assert report['skins']==1
            report['minimum_height_by_pose']={str(frame):float(min(points[:,2].min() for points in snapshot.values())) for frame,snapshot in snapshots.items()}
            body=bpy.data.objects['body_man']
            # The navel-side patch must follow the torso, never the swinging thigh.
            abdominal_vertices=[v for v in body.data.vertices
                if abs(v.co.x)<.08 and v.co.y<-.08 and .08<v.co.z<.18]
            assert abdominal_vertices
            thigh_groups={g.index for g in body.vertex_groups
                if g.name.startswith('DEF-thigh.')}
            abdominal_thigh_weight=max(sum(g.weight for g in v.groups
                if g.group in thigh_groups) for v in abdominal_vertices)
            assert abdominal_thigh_weight<1e-7,abdominal_thigh_weight
            report['maximum_navel_patch_thigh_weight']=abdominal_thigh_weight
            calf_vertices=[v for v in body.data.vertices if -.9<v.co.z<-.65]
            calf_thigh_weight=max(sum(g.weight for g in v.groups if g.group in thigh_groups)
                                  for v in calf_vertices)
            assert calf_thigh_weight<1e-7,calf_thigh_weight
            report['maximum_lower_calf_thigh_weight']=calf_thigh_weight
            indices={side:[v.index for v in body.data.vertices if v.co.z<-.98 and v.co.x*sign>0]
                     for sign,side in [(1,'L'),(-1,'R')]}
            # Identify heel and toe patches in the undeformed source, then
            # measure their actual skinned surfaces. Control pitch alone does
            # not prove heel-first contact or that the IK chain can reach it.
            foot_regions={}
            for side,foot_indices in indices.items():
                rest_y=np.array([body.data.vertices[i].co.y for i in foot_indices])
                foot_regions[side]={
                    'heel':[i for i,y in zip(foot_indices,rest_y) if y>=np.quantile(rest_y,.8)],
                    'toe':[i for i,y in zip(foot_indices,rest_y) if y<=np.quantile(rest_y,.2)]}
            contacts=[]
            knee_lateral=[]
            knee_forward=[]
            ankle_spacing=[]
            toe_out=[]
            body_lean=[]
            waist_bend=[]
            pelvis_heights={}
            heel_strikes=[]
            rig=bpy.data.objects['rig']
            # Regress the reported horizontal pinch using geometry, independent
            # of smooth shading. The intended flexion crease lies below this band.
            rest=np.array([tuple(v.co) for v in body.data.vertices])
            triangles=np.array([tuple(p.vertices) for p in body.data.polygons])
            edges={}
            for face in body.data.polygons:
                for edge in face.edge_keys:edges.setdefault(edge,[]).append(face.index)
            band=[]
            for edge,faces in edges.items():
                x,y,z=rest[list(edge)].mean(axis=0)
                if len(faces)==2 and -.17<x<-.04 and y<-.045 and -.075<z<-.015:
                    band.append(faces)
            assert band
            band=np.array(band)
            fold_angles=[]
            for frame in range(1,49):
                scene.frame_set(frame)
                points=evaluated_vertices(root)['body_man']
                corners=points[triangles[band]]
                normals=np.cross(corners[:,:,1]-corners[:,:,0],corners[:,:,2]-corners[:,:,0])
                normals/=np.linalg.norm(normals,axis=2)[:,:,None]
                angles=np.degrees(np.arccos(np.clip(np.sum(normals[:,0]*normals[:,1],axis=1),-1,1)))
                fold_angles.append(float(angles.max()))
                shoulders=(rig.pose.bones['ORG-upper_arm.L'].head+rig.pose.bones['ORG-upper_arm.R'].head)/2
                hips=(rig.pose.bones['ORG-thigh.L'].head+rig.pose.bones['ORG-thigh.R'].head)/2
                torso=shoulders-hips
                body_lean.append(math.degrees(math.atan2(-torso.y,torso.z)))
                pelvis_heights[frame]=hips.z
                lower=rig.pose.bones['DEF-spine.001']
                upper=rig.pose.bones['DEF-spine.002']
                waist_bend.append(math.degrees((lower.tail-lower.head).angle(upper.tail-upper.head)))
                ankles=[]
                for side in ['L','R']:
                    hip=rig.pose.bones['ORG-thigh.'+side].head
                    knee=rig.pose.bones['ORG-shin.'+side].head
                    ankle=rig.pose.bones['ORG-foot.'+side].head
                    leg=ankle-hip
                    bend=knee-hip-leg*((knee-hip).dot(leg)/leg.length_squared)
                    knee_lateral.append(abs(bend.x))
                    knee_forward.append(-bend.y)
                    ankles.append(ankle.x)
                    direction=rig.pose.bones['ORG-foot.'+side].tail-ankle
                    toe_out.append(abs(direction.x/direction.y))
                    if (frame,side) in [(1,'L'),(25,'R')]:
                        heel=float(points[foot_regions[side]['heel'],2].min())
                        toe=float(points[foot_regions[side]['toe'],2].min())
                        reach=hip.y-ankle.y
                        knee_angle=math.degrees((knee-hip).angle(ankle-knee))
                        other='R' if side=='L' else 'L'
                        step=rig.pose.bones['ORG-foot.'+other].head.y-ankle.y
                        assert 0<heel<.003 and toe-heel>.06,(frame,heel,toe)
                        assert reach>.35 and step>.70,(frame,reach,step)
                        assert knee_angle<30,(frame,knee_angle)
                        heel_strikes.append({'frame':frame,'side':side,
                            'heel_height':heel,'toe_height':toe,
                            'forward_reach':reach,'step_length':step,
                            'knee_flexion_degrees':knee_angle})
                    if (frame,side) in [(7,'L'),(31,'R')]:
                        toe=float(points[foot_regions[side]['toe'],2].min())
                        assert 0<toe<.003,(frame,toe)
                ankle_spacing.append(ankles[0]-ankles[1])
                for side,offset in [('L',0),('R',.5)]:
                    if ((frame-1)/48+offset)%1<.6:
                        contacts.append(float(points[indices[side],2].min()))
            assert max(fold_angles)<30,fold_angles
            report['maximum_old_pelvic_band_dihedral_degrees']=max(fold_angles)
            assert min(contacts)>0 and max(contacts)<.003,contacts
            report['stance_sole_height_range']=[min(contacts),max(contacts)]
            # Evaluate the actual solved joints, not the authored pole values:
            # the old rig could have a forward target but still bend sideways.
            assert max(knee_lateral)<.008,knee_lateral
            assert min(knee_forward)>-.001 and max(knee_forward)>.10,knee_forward
            assert .14<min(ankle_spacing) and max(ankle_spacing)<.20,ankle_spacing
            assert max(toe_out)<.15,toe_out
            report['maximum_lateral_knee_bend']=max(knee_lateral)
            report['forward_knee_bend_range']=[min(knee_forward),max(knee_forward)]
            report['ankle_spacing_range']=[min(ankle_spacing),max(ankle_spacing)]
            report['maximum_toe_out_ratio']=max(toe_out)
            assert max(abs(angle) for angle in body_lean)<.1,body_lean
            bob=max(pelvis_heights.values())-min(pelvis_heights.values())
            assert .035<bob<.05,bob
            report['pelvis_vertical_excursion']=bob
            for start in [1,25]:
                assert pelvis_heights[start+6]<pelvis_heights[start]<pelvis_heights[start+12]<pelvis_heights[start+18]
            # The old pose had a roughly 32-degree hinge at the stomach.
            # Check the solved joints, including the animated counter-twist.
            assert max(waist_bend)<16,waist_bend
            report['waist_bend_degrees']=[min(waist_bend),max(waist_bend)]
            report['heel_strikes']=heel_strikes
            report['forward_torso_lean_degrees']=[min(body_lean),max(body_lean)]
            rig=bpy.data.objects['rig']
            scene.frame_set(1);first=rig.pose.bones['DEF-hand.L'].head.copy()
            scene.frame_set(25);last=rig.pose.bones['DEF-hand.L'].head.copy()
            report['hand_swing_distance']=(last-first).length
            assert report['hand_swing_distance']>.45
        if name=='jellyfish':assert report['morph_targets']>0
        results[name]=report
    (options.output/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
