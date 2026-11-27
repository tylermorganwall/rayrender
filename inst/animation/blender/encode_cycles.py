"""Encode and decode-check three loop movies and a side-by-side preview."""
import argparse
import json
from pathlib import Path
import subprocess

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
movies=[]
for name in ['bird','jellyfish','walk']:
    directory=args.output/name
    frames=list((directory/'frames').glob('frame-*.png'))
    assert len(frames)==48,(name,len(frames))
    cycle=directory/'cycle.mp4'
    movie=directory/(name+'.mp4')
    subprocess.run(['ffmpeg','-v','error','-y','-framerate','24','-start_number','1',
        '-i',str(directory/'frames'/'frame-%04d.png'),'-frames:v','48',
        '-c:v','libx264','-crf','18','-pix_fmt','yuv420p',str(cycle)],check=True)
    subprocess.run(['ffmpeg','-v','error','-y','-stream_loop','2','-i',str(cycle),
        '-c','copy','-movflags','+faststart',str(movie)],check=True)
    movies.append(movie)
combined=args.output/'animation-cycles.mp4'
command=['ffmpeg','-v','error','-y']
for movie in movies:command+=['-i',str(movie)]
command+=['-filter_complex','[0:v][1:v][2:v]hstack=inputs=3[v]','-map','[v]',
    '-c:v','libx264','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(combined)]
subprocess.run(command,check=True)
report={}
for movie in movies+[combined]:
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-select_streams','v:0',
        '-show_entries','stream=width,height,nb_frames,r_frame_rate:format=duration',
        '-of','json',str(movie)]))
    assert int(probe['streams'][0]['nb_frames'])==144
    assert abs(float(probe['format']['duration'])-6)<.02
    subprocess.run(['ffmpeg','-v','error','-i',str(movie),'-f','null','-'],check=True)
    report[movie.name]=probe
(args.output/'movie-validation.json').write_text(json.dumps(report,indent=2)+'\n')
print('All four movies have 144 frames, six-second duration, and decode cleanly.')
