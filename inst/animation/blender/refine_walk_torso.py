"""Refine the walking pose and skin weights on the supplied Rigify character."""
import math

import bpy
from mathutils import Matrix, Quaternion


def align_walking_spine(rig, first_frame=1, last_frame=49):
    """Reduce lumbar hinging while retaining the pelvis path and standing upright.

    The source's hips/chest controls distribute their rotation across two
    spine joints each. Counter-tilting the pelvis 25 degrees and the chest
    12 degrees softens the abrupt waist bend. Reorient the whole torso to
    zero forward lean, then restore each original hip-center position.
    Feet and knee poles remain in their authored armature-space positions.
    All original samples are captured before editing, so evaluation order
    cannot feed a corrected pose into the next frame's target.
    """
    scene = bpy.context.scene
    torso = rig.pose.bones['torso']
    hips = rig.pose.bones['hips']
    chest = rig.pose.bones['chest']
    head = rig.pose.bones['head']
    samples = []
    for frame in range(first_frame, last_frame + 1):
        scene.frame_set(frame)
        hip_center = (rig.pose.bones['ORG-thigh.L'].head +
                      rig.pose.bones['ORG-thigh.R'].head) / 2
        samples.append((frame, hip_center.copy(), torso.location.copy(),
                        torso.rotation_quaternion.copy(),
                        hips.rotation_euler.to_quaternion(),
                        chest.rotation_euler.to_quaternion(),
                        head.matrix.to_quaternion()))
    for frame, hip_target, location, rotation, hip_rotation, chest_rotation, head_rotation in samples:
        scene.frame_set(frame)
        torso.location = location
        torso.rotation_quaternion = rotation
        for bone, original, degrees in [(hips, hip_rotation, -25),
                                        (chest, chest_rotation, 12)]:
            rest = bone.bone.matrix_local.to_quaternion()
            bone.rotation_mode = 'QUATERNION'
            correction = Quaternion((1, 0, 0), math.radians(degrees))
            bone.rotation_quaternion = (rest.inverted() @ correction @ rest) @ original
        bpy.context.view_layer.update()
        hip_center = (rig.pose.bones['ORG-thigh.L'].head +
                      rig.pose.bones['ORG-thigh.R'].head) / 2
        shoulders = (rig.pose.bones['ORG-upper_arm.L'].head +
                     rig.pose.bones['ORG-upper_arm.R'].head) / 2
        direction = shoulders - hip_center
        lean = math.atan2(-direction.y, direction.z)
        rest = torso.bone.matrix_local.to_quaternion()
        correction = Quaternion((1, 0, 0), -lean)
        torso.rotation_quaternion = (rest.inverted() @ correction @ rest) @ rotation
        bpy.context.view_layer.update()
        hip_center = (rig.pose.bones['ORG-thigh.L'].head +
                      rig.pose.bones['ORG-thigh.R'].head) / 2
        # Pose matrices are armature-space; assigning the translated control
        # matrix lets Blender convert back through the control's parent basis.
        pose = torso.matrix.copy()
        pose.translation += hip_target - hip_center
        torso.matrix = pose
        torso.keyframe_insert('location', frame=frame)
        torso.keyframe_insert('rotation_quaternion', frame=frame)
        hips.keyframe_insert('rotation_quaternion', frame=frame)
        chest.keyframe_insert('rotation_quaternion', frame=frame)
        bpy.context.view_layer.update()
        # Keep the original gaze level as the chest changes underneath it.
        head.matrix = Matrix.LocRotScale(head.matrix.translation, head_rotation,
                                         head.matrix.to_scale())
        head.keyframe_insert('rotation_quaternion', frame=frame)
        bpy.context.view_layer.update()
    scene.frame_set(first_frame)


def repair_abdominal_weights(body):
    """Remove stray thigh influence above the hip on the authored male mesh.

    The source's automatic weights attach a patch beside the navel to the
    right thigh. Use rest-space coordinates (Z up, front -Y, model units) to
    blend that patch back into its existing spine/pelvis influences. Smooth
    transition bands retain the groin and outer hip weights. Total vertex
    weight, rest coordinates, topology, and every bone pose are preserved.
    """
    groups = {group.index: group for group in body.vertex_groups}
    changed = 0
    for vertex in body.data.vertices:
        x, y, z = vertex.co
        # Full correction above z=.06, inside the abdomen and in front of the
        # spine. Fade over the groin, flank, and rear boundaries, not in time.
        bands = [(z + .07) / .13, (.22 - abs(x)) / .06, (.09 - y) / .07]
        mask = 1.0
        for value in bands:
            t = max(0.0, min(1.0, value))
            mask *= t * t * (3 - 2 * t)
        if mask == 0:
            continue
        weights = {g.group: g.weight for g in vertex.groups}
        thighs = [index for index in weights if groups[index].name.startswith('DEF-thigh.')]
        torso = [index for index in weights if groups[index].name.startswith(('DEF-spine', 'DEF-pelvis.'))]
        removed = sum(weights[index] * mask for index in thighs)
        retained = sum(weights[index] for index in torso)
        if removed < 1e-8 or retained < 1e-8:
            continue
        for index in thighs:
            weight = weights[index] * (1 - mask)
            if weight == 0:
                groups[index].remove([vertex.index])
            else:
                groups[index].add([vertex.index], weight, 'REPLACE')
        for index in torso:
            weight = weights[index] + removed * weights[index] / retained
            groups[index].add([vertex.index], weight, 'REPLACE')
        changed += 1
    return changed


def repair_knee_weights(body, rig):
    """Remove thigh influence leaking down the calf in the source weights.

    Redistribute each leg's existing thigh/shin total through a smooth knee
    transition in rest space. The old asymmetric weights still assigned over
    half of some calf vertices to the thigh, folding triangles during flexion.
    Keep the existing split between each limb's two deform joints, and leave
    foot, opposite-leg, and torso influences unchanged.
    """
    for side in ('L', 'R'):
        thigh = [body.vertex_groups['DEF-thigh.' + side],
                 body.vertex_groups['DEF-thigh.' + side + '.001']]
        shin = [body.vertex_groups['DEF-shin.' + side],
                body.vertex_groups['DEF-shin.' + side + '.001']]
        knee_z = rig.data.bones['DEF-shin.' + side].head_local.z
        for vertex in body.data.vertices:
            if vertex.co.z >= knee_z + .20:
                continue
            weights = {g.group: g.weight for g in vertex.groups}
            thigh_weights = [weights.get(g.index, 0) for g in thigh]
            shin_weights = [weights.get(g.index, 0) for g in shin]
            thigh_total, shin_total = sum(thigh_weights), sum(shin_weights)
            total = thigh_total + shin_total
            if total < 1e-8:
                continue
            # An eight-centimeter band on each side of the knee. Fade the
            # correction into the unchanged upper-thigh weights above it.
            t = max(0.0, min(1.0, (vertex.co.z - knee_z + .08) / .16))
            thigh_fraction = t * t * (3 - 2 * t)
            t = max(0.0, min(1.0, (knee_z + .20 - vertex.co.z) / .10))
            blend = t * t * (3 - 2 * t)
            target_thigh = thigh_total + blend * (total * thigh_fraction - thigh_total)
            for groups, previous, subtotal, target, fallback in [
                (thigh, thigh_weights, thigh_total, target_thigh, 1),
                (shin, shin_weights, shin_total, total - target_thigh, 0)]:
                for index, group in enumerate(groups):
                    share = previous[index] / subtotal if subtotal > 1e-8 else float(index == fallback)
                    weight = target * share
                    if weight <= 0:
                        if group.index in weights:
                            group.remove([vertex.index])
                    else:
                        group.add([vertex.index], weight, 'REPLACE')


def stabilize_walk_triangles(body, first_frame=1, last_frame=49):
    """Improve fixed quad diagonals using poses throughout this clip.

    glTF stores triangles, while Blender can retessellate deformed quads at
    every pose. Thin source quads around the knee can invert one exported
    triangle when a rest-pose diagonal is used. Score both diagonals by their
    worst triangle/skinned-normal agreement over the authored poses, including
    the reported fractional pose. Commit one fixed topology for native skinning.
    """
    import bmesh
    import numpy as np

    scene = bpy.context.scene
    quads = np.array([tuple(p.vertices) for p in body.data.polygons])
    if quads.shape[1] != 4:
        raise ValueError('Walk diagonal selection expects the subdivided quad mesh.')
    worst = np.ones((len(quads), 2))
    rig = body.find_armature()
    joints = [bone for bone in rig.pose.bones if body.vertex_groups.get(bone.name)]
    columns = {body.vertex_groups[bone.name].index: i for i, bone in enumerate(joints)}
    weights = np.zeros((len(body.data.vertices), len(joints)))
    rest_normals = np.array([tuple(v.normal) for v in body.data.vertices])
    for vertex in body.data.vertices:
        for group in vertex.groups:
            if group.group in columns:
                weights[vertex.index, columns[group.group]] = group.weight
    weights /= weights.sum(axis=1)[:, None]
    object_to_rig = rig.matrix_world.inverted() @ body.matrix_world
    rig_to_object = object_to_rig.inverted()
    poses = list(range(first_frame, last_frame + 1)) + [15.8]
    for time in poses:
        frame = math.floor(time)
        scene.frame_set(frame, subframe=time-frame)
        evaluated = body.evaluated_get(bpy.context.evaluated_depsgraph_get())
        mesh = evaluated.to_mesh()
        points = np.array([tuple(v.co) for v in mesh.vertices])[quads]
        evaluated.to_mesh_clear()
        # Use the same inverse-transpose blend as native glTF skinning. A
        # deformed quad normal alone can already have flipped on a thin quad.
        palette = np.array([rig_to_object @ bone.matrix @ bone.bone.matrix_local.inverted()
                            @ object_to_rig for bone in joints])[:, :3, :3]
        blended = np.einsum('ij,jkl->ikl', weights, palette)
        normals = np.linalg.solve(blended.transpose(0, 2, 1), rest_normals[:, :, None])[:, :, 0]
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        for option, triangles in enumerate([[(0, 1, 2), (0, 2, 3)],
                                             [(0, 1, 3), (1, 2, 3)]]):
            for a, b, c in triangles:
                reference = normals[quads[:, a]] + normals[quads[:, b]] + normals[quads[:, c]]
                reference /= np.maximum(np.linalg.norm(reference, axis=1)[:, None], 1e-20)
                normal = np.cross(points[:, b]-points[:, a], points[:, c]-points[:, a])
                normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-20)
                worst[:, option] = np.minimum(worst[:, option], np.sum(normal*reference, axis=1))
    # Ties retain the first diagonal for deterministic source generation.
    alternate = worst[:, 1] > worst[:, 0]
    scene.frame_set(first_frame)
    normals = [tuple(v.normal) for v in body.data.vertices]
    mesh = bmesh.new()
    mesh.from_mesh(body.data)
    mesh.faces.ensure_lookup_table()
    first = [mesh.faces[i] for i in range(len(quads)) if not alternate[i]]
    second = [mesh.faces[i] for i in range(len(quads)) if alternate[i]]
    bmesh.ops.triangulate(mesh, faces=first, quad_method='FIXED')
    bmesh.ops.triangulate(mesh, faces=second, quad_method='ALTERNATE')
    mesh.to_mesh(body.data)
    mesh.free()
    body.data.update()
    # Retain the original smooth rest normals instead of shading the new
    # diagonal as an authored crease. Consumers can retain these authored
    # normals or explicitly recalculate smooth normals from the posed surface.
    body.data.normals_split_custom_set_from_vertices(normals)
    return {'alternate_diagonals': int(alternate.sum()),
            'worst_selected_normal_agreement': float(np.max(worst, axis=1).min())}


def repair_pelvic_weights(body):
    """Place the right hip's flexion transition along the pelvic crease.

    The source attaches front abdominal vertices to the raised thigh across a
    nearly horizontal band. Paint a diagonal groin-to-outer-hip transition in
    rest coordinates (Z up, model units), blending the rigid pelvis and proximal
    thigh across a .24-unit band. Smooth spatial boundaries
    retain the left outer hip, back, knee, and abdomen outside this correction.
    Rest vertices and all bone transforms remain unchanged. Normalize the local
    weights so Blender does not clamp a consolidated influence above one.
    """
    groups = {g.name: g for g in body.vertex_groups}
    changed = 0
    for vertex in body.data.vertices:
        x, y, z = vertex.co
        bands = [(.06-x)/.10, (.30+x)/.06, (z+.36)/.10,
                 (.16-z)/.10, (.05-y)/.08]
        blend = 1.0
        for value in bands:
            t = max(0.0, min(1.0, value))
            blend *= t*t*(3-2*t)
        if blend == 0:
            continue
        current = {body.vertex_groups[g.group].name:g.weight for g in vertex.groups}
        current_total = sum(current.values())
        current = {n:w/current_total for n,w in current.items()}
        affected = {n:w for n,w in current.items()
                    if n.startswith(('DEF-pelvis','DEF-spine','DEF-thigh'))}
        total = sum(affected.values())
        if total < 1e-8:
            continue
        # The crease rises toward the outer hip instead of cutting horizontally
        # across the abdomen. Above it, the surface follows the pelvis; below
        # it, the thigh. Cubic smoothstep avoids a discontinuity in skin weights.
        line = -.17 + .65*abs(x)
        t = max(0.0, min(1.0, (line+.12-z)/.24))
        thigh_fraction = t*t*(3-2*t)
        target = {n:w for n,w in current.items() if n not in affected}
        target.update({'DEF-pelvis.R': total*(1-thigh_fraction),
                       'DEF-thigh.R': total*thigh_fraction})
        for name in current.keys() | target.keys():
            weight = (1-blend)*current.get(name,0) + blend*target.get(name,0)
            if weight > 0:
                groups[name].add([vertex.index], weight, 'REPLACE')
            elif name in current:
                groups[name].remove([vertex.index])
        changed += 1
    return {'vertices': changed}
