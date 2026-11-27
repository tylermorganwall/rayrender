"""Generate original tiny analytic glTF/GLB deformation fixtures."""
import base64
import copy
import json
from pathlib import Path
import struct

output = Path(__file__).resolve().parent
blob = bytearray()
views = []
accessors = []


def accessor(values, width, kind='VEC3', component=5126):
    while len(blob) % 4:
        blob.append(0)
    offset = len(blob)
    code = {5126: 'f', 5123: 'H', 5121: 'B'}[component]
    blob.extend(struct.pack('<'+code*len(values), *values))
    views.append({'buffer': 0, 'byteOffset': offset, 'byteLength': len(blob)-offset})
    accessors.append({'bufferView': len(views)-1, 'componentType': component,
                      'count': len(values)//width, 'type': kind})
    return len(accessors)-1


positions = accessor([0,0,0, 1,0,0, 0,1,0],3)
normals = accessor([0,0,1]*3,3)
joints0 = accessor([0,1,0,1]*3,4,'VEC4',5121)
weights0 = accessor([.2,.3,0,0]*3,4,'VEC4')
joints1 = accessor([1,0,1,0]*3,4,'VEC4',5121)
weights1 = accessor([.4,.1,0,0]*3,4,'VEC4')
delta0 = accessor([0,1,0]*3,3)
delta1 = accessor([0,0,1]*3,3)
identity = [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]
bind1 = identity.copy(); bind1[13] = -1
binds = accessor(identity+bind1,16,'MAT4')
time = accessor([0,2],1,'SCALAR')
joint_motion = accessor([3,1,0, 5,1,0],3)
mesh_motion = accessor([3,0,0, 7,0,0],3)
weight_time = accessor([.25,1.25],1,'SCALAR')
weight_values = accessor([0,0, 1,.5],1,'SCALAR')
# Two targets per key: incoming, value, outgoing, each in target order.
cubic_values = accessor([0,0, .2,.3, .4,-.2, 0,0, .8,.1, 0,0],1,'SCALAR')
base = {'asset': {'version':'2.0','generator':'rayrender original analytic deformation fixture'},
        'scene':0,'scenes':[{'nodes':[0,1,2]}],
        'nodes':[{'name':'mesh','mesh':0,'translation':[3,0,0],'skin':0,'weights':[.5,0]},
                 {'name':'joint0','translation':[3,0,0]},
                 {'name':'joint1','translation':[3,1,0]}],
        'meshes':[{'weights':[.25,.5],'primitives':[{'attributes':{
            'POSITION':positions,'NORMAL':normals,'JOINTS_0':joints0,'WEIGHTS_0':weights0,
            'JOINTS_1':joints1,'WEIGHTS_1':weights1},'targets':[{'POSITION':delta0},{'POSITION':delta1}]}]}],
        'skins':[{'joints':[1,2],'inverseBindMatrices':binds}],
        'animations':[{'name':'combined','samplers':[
            {'input':time,'output':joint_motion,'interpolation':'LINEAR'},
            {'input':time,'output':mesh_motion,'interpolation':'LINEAR'},
            {'input':weight_time,'output':weight_values,'interpolation':'LINEAR'}],
            'channels':[{'sampler':0,'target':{'node':2,'path':'translation'}},
                        {'sampler':1,'target':{'node':0,'path':'translation'}},
                        {'sampler':2,'target':{'node':0,'path':'weights'}}]},
          {'name':'cubic morph','samplers':[{'input':time,'output':cubic_values,'interpolation':'CUBICSPLINE'}],
           'channels':[{'sampler':0,'target':{'node':0,'path':'weights'}}]},
          {'name':'step morph','samplers':[{'input':weight_time,'output':weight_values,'interpolation':'STEP'}],
           'channels':[{'sampler':0,'target':{'node':0,'path':'weights'}}]}],
        'bufferViews':views,'accessors':accessors}


def write(name, document):
    document=copy.deepcopy(document)
    document['buffers']=[{'byteLength':len(blob),'uri':'data:application/octet-stream;base64,'+base64.b64encode(blob).decode()}]
    (output/(name+'.gltf')).write_text(json.dumps(document,separators=(',',':'))+'\n')
    del document['buffers'][0]['uri']
    payload=json.dumps(document,separators=(',',':')).encode()
    payload+=b' '*((-len(payload))%4)
    binary=bytes(blob)+b'\0'*((-len(blob))%4)
    data=struct.pack('<III',0x46546c67,2,28+len(payload)+len(binary))
    data+=struct.pack('<II',len(payload),0x4e4f534a)+payload
    data+=struct.pack('<II',len(binary),0x004e4942)+binary
    (output/(name+'.glb')).write_bytes(data)


write('deform-combined',base)
only_skin=copy.deepcopy(base)
del only_skin['meshes'][0]['primitives'][0]['targets']
del only_skin['meshes'][0]['weights']; del only_skin['nodes'][0]['weights']
only_skin['animations']=only_skin['animations'][:1]
only_skin['animations'][0]['channels']=only_skin['animations'][0]['channels'][:2]
only_skin['animations'][0]['samplers']=only_skin['animations'][0]['samplers'][:2]
write('deform-skin',only_skin)
only_morph=copy.deepcopy(base)
del only_morph['skins']; del only_morph['nodes'][0]['skin']
for key in ['JOINTS_0','JOINTS_1','WEIGHTS_0','WEIGHTS_1']:
    del only_morph['meshes'][0]['primitives'][0]['attributes'][key]
write('deform-morph',only_morph)

bend_position=accessor([0,0,0, 0,0,1, 0,0,0],3)
bend_normal=accessor([-1,0,0]*3,3)
bend_weight=accessor([0,1],1,'SCALAR')
bend=copy.deepcopy(only_morph)
bend['accessors']=accessors; bend['bufferViews']=views
bend['nodes'][0]['weights']=[0]; bend['meshes'][0]['weights']=[0]
bend['meshes'][0]['primitives'][0]['targets']=[{'POSITION':bend_position,'NORMAL':bend_normal}]
bend['animations']=[{'name':'bend','samplers':[{'input':time,'output':bend_weight,'interpolation':'LINEAR'}],
                     'channels':[{'sampler':0,'target':{'node':0,'path':'weights'}}]}]
write('deform-bend',bend)
