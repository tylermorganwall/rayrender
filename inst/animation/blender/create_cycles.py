"""Create editable Blender cycles and glTF assets, without changing the source.

Run with Blender --background --python create_cycles.py -- --output DIRECTORY
--rayman /path/to/rayman.blend. Blender 5.2 was used for validation.
Coordinates are Z up, with the bird and walker facing -Y. Internal cycles stay
in place: downstream placement paths are independent of these actions.
"""
import argparse
import json
import math
from pathlib import Path
import sys

import bmesh
import bpy
from mathutils import Matrix, Quaternion, Vector

TAU = 2 * math.pi
FPS = 24
PERIOD = 48


def material(name, color, roughness=0.4, metallic=0, transmission=0):
    result = bpy.data.materials.new(name)
    result.diffuse_color = (*color, 1)
    result.use_nodes = True
    shader = result.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Base Color'].default_value = (*color, 1)
    shader.inputs['Roughness'].default_value = roughness
    shader.inputs['Metallic'].default_value = metallic
    shader.inputs['Transmission Weight'].default_value = transmission
    return result


def mesh_object(name, vertices, faces, mat, parent=None):
    data = bpy.data.meshes.new(name)
    data.from_pydata(vertices, [], faces)
    data.update()
    # Recalculate coherent outward normals on closed surfaces, rather than
    # relying on mirrored winding or arbitrary per-triangle vertex normals.
    bm = bmesh.new()
    bm.from_mesh(data)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bm.to_mesh(data)
    bm.free()
    obj = bpy.data.objects.new(name, data)
    bpy.context.collection.objects.link(obj)
    obj.parent = parent
    data.materials.append(mat)
    for face in data.polygons:
        face.use_smooth = True
    return obj


def empty(name, position=(0, 0, 0), parent=None):
    obj = bpy.data.objects.new(name, None)
    bpy.context.collection.objects.link(obj)
    obj.location = position
    obj.parent = parent
    return obj


def ellipsoid(name, center, scale, mat, parent=None):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=20)
    obj = bpy.context.object
    obj.name = name
    # Bake positive scale into vertices. Node animation then has a clean rest
    # transform and every smooth normal describes the actual ellipsoid.
    for vertex in obj.data.vertices:
        vertex.co = Vector((vertex.co.x * scale[0], vertex.co.y * scale[1], vertex.co.z * scale[2]))
    obj.location = center
    obj.parent = parent
    obj.data.materials.append(mat)
    for face in obj.data.polygons:
        face.use_smooth = True
    return obj


def feather(name, start, end, width, mat, parent):
    start, end = Vector(start), Vector(end)
    obj = ellipsoid(name, (start + end) / 2, (width, (end - start).length / 2, .035), mat, parent)
    obj.rotation_euler = (end - start).to_track_quat('Y', 'Z').to_euler()
    return obj


def tube(name, points, radii, mat, parent=None, sides=10):
    vertices, faces = [], []
    # Parallel transported frames avoid sudden cross-section flips along curves.
    previous = Vector((1, 0, 0))
    for i, point in enumerate(points):
        tangent = (Vector(points[min(i + 1, len(points) - 1)]) - Vector(points[max(i - 1, 0)])).normalized()
        normal = previous - tangent * previous.dot(tangent)
        if normal.length < 1e-6:
            normal = tangent.cross(Vector((0, 1, 0)))
        normal.normalize()
        binormal = tangent.cross(normal).normalized()
        previous = normal
        for j in range(sides):
            angle = TAU * j / sides
            vertices.append(Vector(point) + radii[i] * (math.cos(angle) * normal + math.sin(angle) * binormal))
    for i in range(len(points) - 1):
        for j in range(sides):
            a = i * sides + j
            b = i * sides + (j + 1) % sides
            faces.append((a, b, b + sides, a + sides))
    faces += [tuple(reversed(range(sides))), tuple((len(points) - 1) * sides + j for j in range(sides))]
    return mesh_object(name, vertices, faces, mat, parent)


def cycle_keys():
    # Each cycle includes an identical endpoint at frame 49. Playback/movie
    # frames 1..48 avoid duplicating that endpoint. Sampled keys are smoothed
    # with clamped handles; the Cycles modifier supports arbitrary seeking.
    for action in bpy.data.actions:
        if not action.name.startswith(('Bird', 'Jellyfish', 'Walk')):
            continue
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:
                            key.interpolation = 'BEZIER'
                            key.handle_left_type = 'AUTO_CLAMPED'
                            key.handle_right_type = 'AUTO_CLAMPED'
                        curve.modifiers.new('CYCLES')


def name_action(obj, name):
    obj.animation_data.action.name = name


def bird():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    root = empty('Bird_Placement')
    body = empty('Bird_Body', parent=root)
    blue = material('Kingfisher - teal plumage', (.025, .27, .34))
    dark = material('Flight feathers - deep ink', (.016, .055, .09))
    light = material('Coverts - blue highlights', (.05, .43, .53))
    cream = material('Warm ivory breast', (.87, .79, .55))
    orange = material('Ochre beak and feet', (.83, .29, .045))
    black = material('Glossy black eyes', (.004, .009, .013), .12)
    ellipsoid('Rounded body', (0, .12, 0), (.33, .60, .35), blue, body)
    ellipsoid('Cream breast', (0, -.19, -.095), (.275, .43, .27), cream, body)
    ellipsoid('Head', (0, -.47, .28), (.265, .29, .27), blue, body)
    ellipsoid('Throat', (0, -.60, .13), (.18, .19, .16), cream, body)
    for side in [-1, 1]:
        ellipsoid('Eye ring', (side * .233, -.57, .33), (.032, .071, .074), cream, body)
        ellipsoid('Eye', (side * .259, -.585, .34), (.026, .041, .043), black, body)
        ellipsoid('Eye glint', (side * .278, -.599, .354), (.006, .011, .011), cream, body)
    bpy.ops.mesh.primitive_cone_add(vertices=32, radius1=.105, radius2=.008, depth=.42)
    beak=bpy.context.object; beak.name='Tapered beak';beak.parent=body
    beak.location=(0,-.87,.255);beak.rotation_euler=(math.pi/2,0,0)
    beak.data.materials.append(orange)
    for face in beak.data.polygons: face.use_smooth=True
    for i in range(7):
        x=(i-3)*.06
        feather('Tail feather %02d'%i,(x,.49,-.02),(x*2.1,1.15-abs(x)*.4,-.11),.075,dark if i%2 else blue,body)
    wings=[]
    for side,label in [(-1,'R'),(1,'L')]:
        shoulder=empty('Shoulder.'+label,(side*.25,.02,.12),body)
        wrist=empty('Wrist.'+label,(side*.66,.04,0),shoulder)
        ellipsoid('Upper wing.'+label,(side*.34,.10,0),(.47,.30,.075),blue,shoulder)
        for i in range(8):
            x=.05+i*.083
            feather('Secondary.%s.%02d'%(label,i),(side*x,.02,0),(side*(x+.10),.51+.08*math.sin(i/8*math.pi),-.018),.076,dark if i%2 else blue,shoulder)
            feather('Covert.%s.%02d'%(label,i),(side*x,-.11,.045),(side*(x+.05),.25,.065),.065,light,shoulder)
        ellipsoid('Hand wing.'+label,(side*.25,.08,0),(.36,.22,.055),blue,wrist)
        for i in range(9):
            x=.08+i*.055
            feather('Primary.%s.%02d'%(label,i),(side*x,-.04,0),(side*(.38+i*.063),.51-i*.032,-.02),.069,dark if i%3 else light,wrist)
        wings.append((side,shoulder,wrist))
        ellipsoid('Tucked foot.'+label,(side*.13,.26,-.31),(.045,.17,.05),orange,body)
    for frame in range(1,PERIOD+2):
        phase=TAU*(frame-1)/PERIOD
        body.location.z=.045*math.sin(phase-.3)
        body.rotation_euler.x=.035*math.sin(phase-.6)
        body.keyframe_insert('location',frame=frame);body.keyframe_insert('rotation_euler',frame=frame)
        for side,shoulder,wrist in wings:
            # Shoulder leads; wrist and feather bank lag through the recovery.
            shoulder.rotation_euler.y=-side*(.13+.70*math.cos(phase))
            shoulder.rotation_euler.z=side*.10*math.sin(phase)
            wrist.rotation_euler.y=-side*.30*math.cos(phase-.65)
            wrist.rotation_euler.z=side*(.10+.16*math.sin(phase-.45))
            shoulder.keyframe_insert('rotation_euler',frame=frame)
            wrist.keyframe_insert('rotation_euler',frame=frame)
    for obj in [body]+[o for _,a,b in wings for o in [a,b]]:name_action(obj,'Bird_Flap_'+obj.name)
    return root, (3.4,-5.6,2.7), (0,.12,.05), 4.4, -.95


def jellyfish():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    root=empty('Jellyfish_Placement')
    swim=empty('Jellyfish_Swim',parent=root)
    bell_mat=material('Bell - lavender glass',(.30,.20,.62),.23,transmission=.25)
    rim_mat=material('Bell edge - luminous cyan',(.07,.62,.76),.25,metallic=.15)
    arm_mat=material('Oral arms - warm coral',(.65,.22,.34),.35)
    tentacle_mat=material('Tentacles - pale lavender',(.57,.44,.80),.3)
    # A closed, thick bell: outer dome down to the rim and back up its inner
    # surface. Separate poles avoid degenerate triangles at the apex.
    n=64; rings=24
    vertices=[(0,0,.73)]; faces=[]
    for inner in [False,True]:
        sequence=range(1,rings+1) if not inner else range(rings,0,-1)
        for j in sequence:
            theta=1.48*j/rings
            for i in range(n):
                angle=TAU*i/n
                radius=(.72 if not inner else .675)*math.sin(theta)
                scallop=.025*math.cos(12*angle)*(j/rings)**6
                z=(.73 if not inner else .685)*math.cos(theta)+scallop
                vertices.append((radius*math.cos(angle),radius*math.sin(angle),z))
    for i in range(n):faces.append((0,1+i,1+(i+1)%n))
    for ring in range(2*rings-1):
        for i in range(n):
            a=1+ring*n+i;b=1+ring*n+(i+1)%n
            faces.append((a,a+n,b+n,b))
    vertices.append((0,0,.685));pole=len(vertices)-1
    for i in range(n):faces.append((pole,1+(2*rings-1)*n+(i+1)%n,1+(2*rings-1)*n+i))
    bell=mesh_object('Thick scalloped bell',vertices,faces,bell_mat,swim)
    bell_parts=[bell]
    points=[(.71*math.cos(TAU*i/128),.71*math.sin(TAU*i/128),.067+.025*math.cos(12*TAU*i/128)) for i in range(129)]
    bell_parts.append(tube('Scalloped rim',points,[.022]*129,rim_mat,swim))
    for i in range(12):
        a=TAU*i/12
        points=[(.728*math.sin(t)*math.cos(a),.728*math.sin(t)*math.sin(a),.738*math.cos(t)) for t in [1.47*j/28+.008 for j in range(29)]]
        bell_parts.append(tube('Radial canal %02d'%i,points,[.008]*29,rim_mat,swim,sides=6))
    pulse_keys=[]
    for obj in bell_parts:
        obj.shape_key_add(name='Basis')
        key=obj.shape_key_add(name='Bell contraction')
        for point in key.data:
            p=point.co
            r=min(1,math.hypot(p.x,p.y)/.72)
            p.x*=1-.24*r;p.y*=1-.24*r;p.z=p.z*1.12+.07*r*r
        pulse_keys.append((obj,key))
    tentacle_keys=[]
    for i in range(16):
        angle=TAU*i/16
        length=1.55+.38*math.sin(i*2.7)**2
        points=[];radii=[]
        for j in range(49):
            t=j/48
            radial=.55+.08*math.sin(7*t+i)*t
            points.append((radial*math.cos(angle)+.09*t*math.sin(9*t+i),radial*math.sin(angle)+.09*t*math.cos(8*t+i),.09-length*t))
            radii.append(.014*(1-.85*t))
        obj=tube('Trailing tentacle %02d'%i,points,radii,tentacle_mat,swim,sides=8)
        obj.shape_key_add(name='Basis')
        keys=[obj.shape_key_add(name=name) for name in ['Wave sine','Wave cosine']]
        for key in keys:key.slider_min=-1
        for index,base in enumerate(obj.data.vertices):
            t=(index//8)/48;lag=5*t+i*.4
            for key,wave in zip(keys,[math.cos(lag),-math.sin(lag)]):
                key.data[index].co.x+=.15*t*t*wave
                key.data[index].co.y+=.09*t*t*wave*math.sin(angle+.7)
        tentacle_keys.append((obj,keys))
    for i in range(4):
        angle=TAU*i/4+.5
        points=[];radii=[]
        for j in range(49):
            t=j/48
            points.append((.20*math.cos(angle)+.09*math.sin(12*t+i)*t,.20*math.sin(angle)+.09*math.cos(12*t+i)*t,-.02-1.05*t))
            radii.append((.065+.035*math.sin(14*t)**2)*(1-.84*t))
        obj=tube('Ruffled oral arm %02d'%i,points,radii,arm_mat,swim,sides=12)
        obj.shape_key_add(name='Basis')
        keys=[obj.shape_key_add(name=name) for name in ['Wave sine','Wave cosine']]
        for key in keys:key.slider_min=-1
        for index,base in enumerate(obj.data.vertices):
            t=(index//12)/48;lag=4*t+i*.5
            for key,wave in zip(keys,[math.cos(lag),-math.sin(lag)]):
                key.data[index].co.x+=.10*t*t*wave
                key.data[index].co.y+=.08*t*t*wave
        tentacle_keys.append((obj,keys))
    # Remap signed waves to weights in [0, 1]. Blender's glTF importer clamps
    # negative animated weights to its imported slider range. An offset basis
    # and doubled deltas represent the same surface without negative weights:
    # (base-ds-dc) + (1+sin)*ds + (1+cos)*dc.
    for obj,keys in tentacle_keys:
        basis=obj.data.shape_keys.key_blocks[0]
        for i,point in enumerate(basis.data):
            original=point.co.copy()
            ds=keys[0].data[i].co-original
            dc=keys[1].data[i].co-original
            point.co=original-ds-dc
            keys[0].data[i].co=original+ds-dc
            keys[1].data[i].co=original-ds+dc
        for key in keys:key.slider_min=0
    for frame in range(1,PERIOD+2):
        phase=TAU*(frame-1)/PERIOD
        # A brief contraction and slower release, followed by trailing motion.
        pulse=((1-math.cos(phase))/2)**2
        swim.location.z=.10*math.sin(phase-.8)
        swim.keyframe_insert('location',frame=frame)
        for obj,key in pulse_keys:
            key.value=pulse;key.keyframe_insert('value',frame=frame)
        for obj,keys in tentacle_keys:
            for key,value in zip(keys,[(1+math.sin(phase))/2,(1+math.cos(phase))/2]):
                key.value=value;key.keyframe_insert('value',frame=frame)
    name_action(swim,'Jellyfish_Swim')
    for obj,key in pulse_keys:name_action(obj.data.shape_keys,'Jellyfish_Pulse_'+obj.name)
    for obj,keys in tentacle_keys:name_action(obj.data.shape_keys,'Jellyfish_Trail_'+obj.name)
    return root,(3,-5,1.9),(0,0,-.55),3.5,-2.15


def gait_curve(phase, keys):
    """Sample cubic Hermite keys (phase, value, derivative per cycle)."""
    for left, right in zip(keys, keys[1:]):
        if phase <= right[0]:
            duration = right[0] - left[0]
            t = (phase - left[0]) / duration
            return ((2*t**3 - 3*t*t + 1)*left[1]
                    + (t**3 - 2*t*t + t)*duration*left[2]
                    + (-2*t**3 + 3*t*t)*right[1]
                    + (t**3 - t*t)*duration*right[2])
    return keys[-1][1]


def walk_foot_pose(phase):
    """Return ankle Y, sole clearance, and X pitch for one in-place stride.

    Facing -Y, negative pitch raises the toes. Contact is at phase 0, and
    liftoff at .6. The passing foot is under the pelvis at .75 rather than
    remaining behind the body. Late swing extends the leg before heel contact.
    """
    if phase < .6:
        y = -.40 + 1.50*phase
        clearance = 0
    else:
        y = gait_curve(phase, [(.6, .50, 1.50), (.75, .04, -4.6),
                               (.875, -.32, -1.8), (1, -.40, 1.50)])
        clearance = .16*math.sin(math.pi*(phase-.6)/.4)**2
    pitch = gait_curve(phase, [(0, -.62, 0), (.10, 0, 0), (.35, 0, 0),
                               (.50, .48, 2), (.60, .75, 0), (.75, .40, -3),
                               (.875, -.25, -3), (.96, -.62, 0), (1, -.62, 0)])
    return y, clearance, pitch


def walking(source):
    bpy.ops.wm.open_mainfile(filepath=str(source))
    scene=bpy.context.scene
    rig=bpy.data.objects['rig'];body=bpy.data.objects['body_man']
    # Keep the source's control rig and start with its weights. Clear action/NLA in
    # this new document only, so no pre-existing pose contaminates a seek.
    for obj in bpy.data.objects:
        if obj.animation_data:
            obj.animation_data.action=None
            for track in list(obj.animation_data.nla_tracks):obj.animation_data.nla_tracks.remove(track)
    for bone in rig.pose.bones:bone.matrix_basis=Matrix.Identity(4)
    # Author this new clip with ordinary linear skinning. glTF has no Blender
    # B-Bone interpolation; leaving segmented B-Bones enabled changes the
    # exported surface. Retain the existing two deform joints per limb; repair
    # abdominal weight leakage below. The original .blend keeps its segmented rig.
    for bone in rig.data.bones:bone.bbone_segments=1
    for side in ['L','R']:
        rig.pose.bones['thigh_parent.'+side]['IK_FK']=0
        rig.pose.bones['thigh_parent.'+side]['IK_Stretch']=0
        # The saved rig has both alternative IK constraints enabled, with
        # poles far outside the legs. Use only the pole-controlled solver;
        # animate its target in the forward knee plane below.
        shin=rig.pose.bones['MCH-shin_ik.'+side]
        # The inherited local-axis locks conflict with a forward pole plane.
        # Let the pole define that plane instead of retaining the splayed hinge.
        shin.lock_ik_y=False
        shin.lock_ik_z=False
        shin.ik_stretch=0
        rig.pose.bones['thigh_ik.'+side].ik_stretch=0
        for constraint in shin.constraints:
            if constraint.type=='IK':
                constraint.influence=1 if constraint.pole_subtarget else 0
                constraint.mute=not bool(constraint.pole_subtarget)
                constraint.use_stretch=False
        rig.pose.bones['upper_arm_parent.'+side]['IK_FK']=1
        # This saved rig has no IK/FK influence drivers. Explicitly enable its
        # existing FK constraints for the arms; changing the custom slider
        # alone leaves the later IK Copy Transforms constraint fully active.
        for part in ['upper_arm','forearm','hand']:
            for constraint in rig.pose.bones['ORG-'+part+'.'+side].constraints:
                if constraint.type=='COPY_TRANSFORMS':
                    if '_ik.' in constraint.subtarget:constraint.influence=0
                    if '_fk.' in constraint.subtarget:constraint.influence=1
    bpy.context.view_layer.update()
    # Use only the supplied male character; Rigify widgets and the alternative
    # mesh remain available in the source file, not in this export scene.
    for obj in list(bpy.data.objects):
        if obj not in [rig,body]:bpy.data.objects.remove(obj,do_unlink=True)
    for obj in [rig,body]:
        obj.hide_set(False);obj.hide_viewport=False;obj.hide_render=False
    body.data.materials.clear()
    body.data.materials.append(material('Rayman - blue ceramic',(.075,.32,.46),.36,.08))
    for face in body.data.polygons:face.use_smooth=True
    bpy.context.view_layer.objects.active=body;body.select_set(True)
    subdiv=body.modifiers.new('Smooth source mesh','SUBSURF');subdiv.levels=2
    bpy.ops.object.modifier_move_up(modifier=subdiv.name)
    bpy.ops.object.modifier_apply(modifier=subdiv.name)
    # Stray thigh weights above the hip pull the navel inward during swing.
    from refine_walk_torso import repair_abdominal_weights, repair_knee_weights, repair_pelvic_weights
    repair_abdominal_weights(body)
    print('PELVIC_WEIGHTS',repair_pelvic_weights(body),flush=True)
    repair_knee_weights(body, rig)
    root=empty('Rayman_Placement');rig.parent=root
    rig.location.z=1.14
    scene.frame_set(1)
    scene.timeline_markers.clear()
    for frame, label in [(1, 'Contact L'), (7, 'Down'), (13, 'Passing'),
                         (19, 'Up'), (25, 'Contact R'), (31, 'Down'),
                         (37, 'Passing'), (43, 'Up'), (49, 'Contact L')]:
        scene.timeline_markers.new(label, frame=frame)
    for frame in range(1,PERIOD+2):
        t=(frame-1)/PERIOD;phase=TAU*t
        torso=rig.pose.bones['torso']
        # Contact, down, passing, up, contact: the high point follows passing.
        height=gait_curve(t % .5, [(0, -.105, -.2), (.125, -.12, 0),
                                 (.25, -.015, .3), (.375, -.008, 0),
                                 (.5, -.105, -.2)])
        # Keep contact height/stride reach, but reduce the vertical excursion
        # to 40 percent of the earlier cycle. Preserve down/passing/up timing.
        height=-.105+.4*(height+.105)
        torso.location=(.018*math.sin(phase),-.035,height)
        # Start from the reference hip trajectory. align_walking_spine brings
        # the shoulders upright afterwards without shifting that trajectory.
        # Conjugate the world-X rotation through the Rigify control's rest basis.
        rest=torso.bone.matrix_local.to_quaternion()
        torso.rotation_mode='QUATERNION'
        torso.rotation_quaternion=rest.inverted() @ Quaternion((1,0,0),.185) @ rest
        torso.keyframe_insert('location',frame=frame)
        torso.keyframe_insert('rotation_quaternion',frame=frame)
        head=rig.pose.bones['head']
        rest=head.bone.matrix_local.to_quaternion()
        head.rotation_mode='QUATERNION'
        head.rotation_quaternion=rest.inverted() @ Quaternion((1,0,0),-.10) @ rest
        head.keyframe_insert('rotation_quaternion',frame=frame)
        for sign,side,offset in [(1,'L',0),(-1,'R',.5)]:
            u=(t+offset)%1
            foot=rig.pose.bones['foot_ik.'+side]
            y, z, pitch = walk_foot_pose(u)
            # Rest ankle Y is .1023, so author the path in armature coordinates
            # instead of unintentionally shifting the entire stride rearward.
            foot.location=(-sign*.085, y-foot.bone.head_local.y, z)
            # Keep toe direction consistent while the foot rolls around X.
            # A constant yaw would splay the toes during the high heel lift.
            rest_foot=rig.data.bones['ORG-foot.'+side]
            direction=Quaternion((1,0,0),pitch) @ (rest_foot.tail_local-rest_foot.head_local)
            yaw=sign*.08-math.atan2(direction.x,-direction.y)
            foot.rotation_mode='XYZ';foot.rotation_euler=(pitch,0,yaw)
            foot.keyframe_insert('location',frame=frame);foot.keyframe_insert('rotation_euler',frame=frame)
            # Lower the source A-pose arms and bias the counter-swing forward.
            # Negative world-X pitch brings a hanging arm toward the front (-Y).
            arm=rig.pose.bones['upper_arm_fk.'+side]
            arm.rotation_mode='QUATERNION'
            rest=arm.bone.matrix_local.to_quaternion()
            world=Quaternion((1,0,0),-.25+.50*math.cos(phase+TAU*offset)) @ Quaternion((0,1,0),sign*.48)
            arm.rotation_quaternion=rest.inverted() @ world @ rest
            arm.keyframe_insert('rotation_quaternion',frame=frame)
            forearm=rig.pose.bones['forearm_fk.'+side]
            rest=forearm.bone.matrix_local.to_quaternion()
            forearm.rotation_mode='QUATERNION'
            forearm.rotation_quaternion=rest.inverted() @ Quaternion((1,0,0),-.18-.10*math.cos(phase+TAU*offset)) @ rest
            forearm.keyframe_insert('rotation_quaternion',frame=frame)
        hips=rig.pose.bones['hips'];hips.rotation_mode='XYZ'
        hips.rotation_euler=(0,.025*math.sin(phase),.045*math.cos(phase))
        hips.keyframe_insert('rotation_euler',frame=frame)
        chest=rig.pose.bones['chest'];chest.rotation_mode='XYZ'
        chest.rotation_euler=(0,0,-.035*math.cos(phase))
        chest.keyframe_insert('rotation_euler',frame=frame)
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        for side in ['L','R']:
            hip=rig.pose.bones['ORG-thigh.'+side].head
            pole=rig.pose.bones['thigh_ik_target.'+side]
            # Set the pole in armature space, then key its converted local
            # location. Its inherited torso transform must not swing it out.
            pose=pole.matrix.copy()
            pose.translation=(hip.x,-1.15,-.50)
            pole.matrix=pose
            pole.keyframe_insert('location',frame=frame)
    # Align the pelvis/rib cage through the existing controls before the final
    # sole fit, retaining the hip trajectory and keeping the torso upright.
    from refine_walk_torso import align_walking_spine
    align_walking_spine(rig)
    # Fit each sole to the ground during stance, using the actual evaluated
    # source mesh, including skinning. During swing preserve the authored lift.
    # Three small corrections account for the source rig's constraint response.
    feet={side:[v.index for v in body.data.vertices if v.co.z<-.98 and v.co.x*sign>0]
          for sign,side in [(1,'L'),(-1,'R')]}
    for frame in range(1,PERIOD+2):
        scene.frame_set(frame)
        for iteration in range(3):
            graph=bpy.context.evaluated_depsgraph_get()
            evaluated=body.evaluated_get(graph)
            mesh=evaluated.to_mesh()
            heights={side:min((evaluated.matrix_world @ mesh.vertices[i].co).z for i in indices)
                     for side,indices in feet.items()}
            evaluated.to_mesh_clear()
            for side,offset in [('L',0),('R',.5)]:
                u=((frame-1)/PERIOD+offset)%1
                _, lift, _ = walk_foot_pose(u)
                foot=rig.pose.bones['foot_ik.'+side]
                foot.location.z+=.002+lift-heights[side]
            bpy.context.view_layer.update()
        for side in ['L','R']:
            rig.pose.bones['foot_ik.'+side].keyframe_insert('location',frame=frame)
    from refine_walk_torso import stabilize_walk_triangles
    print('WALK_TRIANGULATION',stabilize_walk_triangles(body),flush=True)
    name_action(rig,'Walk_Rayman_InPlace')
    return root,(3.1,-4.5,2.3),(0,0,1.03),2.7,0


def stage(camera_position, target, scale, ground):
    scene=bpy.context.scene
    scene.world=bpy.data.worlds.new('Studio world')
    scene.world.use_nodes=True
    scene.world.node_tree.nodes['Background'].inputs[0].default_value=(.12,.16,.21,1)
    scene.world.node_tree.nodes['Background'].inputs[1].default_value=.35
    floor=material('Studio floor',(.055,.075,.10),.65)
    bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,ground))
    obj=bpy.context.object;obj.name='STUDIO_Floor';obj.data.materials.append(floor)
    for name,location,power,size,color in [
        ('Key',(-3,-4,6),700,5,(.80,.91,1)),
        ('Fill',(4,-1,3),420,4,(.61,.80,1)),
        ('Rim',(1,4,5),1000,3,(1,.72,.47))]:
        data=bpy.data.lights.new(name,'AREA');data.energy=power;data.shape='DISK';data.size=size;data.color=color
        light=bpy.data.objects.new('STUDIO_'+name,data);scene.collection.objects.link(light)
        light.location=location;light.rotation_euler=(Vector(target)-light.location).to_track_quat('-Z','Y').to_euler()
    data=bpy.data.cameras.new('Preview camera');camera=bpy.data.objects.new('STUDIO_Camera',data)
    scene.collection.objects.link(camera);camera.location=camera_position
    camera.rotation_euler=(Vector(target)-camera.location).to_track_quat('-Z','Y').to_euler()
    data.type='ORTHO';data.ortho_scale=scale;scene.camera=camera
    scene.render.engine='CYCLES';scene.cycles.samples=32;scene.cycles.use_denoising=True
    scene.render.resolution_x=720;scene.render.resolution_y=720;scene.render.resolution_percentage=100
    scene.render.fps=FPS;scene.frame_start=1;scene.frame_end=PERIOD
    scene.render.image_settings.file_format='PNG'
    scene.view_settings.view_transform='AgX'
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type=='VIEW_3D':
                area.spaces.active.region_3d.view_perspective='CAMERA'
                area.spaces.active.overlay.show_overlays=False
    return scene


def save_asset(name, setup, output, render):
    root,camera,target,scale,ground=setup
    cycle_keys()
    scene=stage(camera,target,scale,ground)
    scene.name={'bird':'Bird_Flap','jellyfish':'Jellyfish_Pulse','walk':'Rayman_Walk'}[name]
    scene.frame_set(1)
    directory=output/name;directory.mkdir(parents=True,exist_ok=True)
    # Export only model descendants, keeping stage/camera outside the asset.
    bpy.ops.object.select_all(action='DESELECT')
    for obj in [root]+list(root.children_recursive):obj.select_set(True)
    bpy.context.view_layer.objects.active=root
    scene.frame_end=PERIOD+1
    bpy.ops.export_scene.gltf(filepath=str(directory/(name+'.glb')),
        export_format='GLB',use_selection=True,export_animations=True,
        export_animation_mode='SCENE',export_anim_scene_split_object=False,
        export_force_sampling=True,export_frame_range=True,
        export_skins=True,export_all_influences=True,export_morph=True,export_morph_normal=True,
        # Flatten only the exported joint hierarchy so Rigify's inheritance
        # rules and stretch transforms are baked into independent joint TRS.
        # The editable Blender rig retains its full control hierarchy.
        export_def_bones=True,export_hierarchy_flatten_bones=True,
        export_normals=True,export_yup=True)
    scene.frame_end=PERIOD
    bpy.ops.wm.save_as_mainfile(filepath=str(directory/(name+'.blend')))
    if render:
        for frame in [1,13,25,37]:
            scene.frame_set(frame);scene.render.filepath=str(directory/('pose-%02d.png'%frame))
            bpy.ops.render.render(write_still=True)
    print('ASSET_COMPLETE',name,flush=True)


def main():
    args=sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else []
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rayman',type=Path,required=True)
    parser.add_argument('--only',choices=['bird','jellyfish','walk'])
    parser.add_argument('--no-render',action='store_true')
    options=parser.parse_args(args)
    sys.dont_write_bytecode=True
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    options.output.mkdir(parents=True,exist_ok=True)
    for name,create in [('bird',bird),('jellyfish',jellyfish),('walk',lambda:walking(options.rayman))]:
        if options.only and options.only!=name:continue
        save_asset(name,create(),options.output,not options.no_render)


if __name__=='__main__':main()
