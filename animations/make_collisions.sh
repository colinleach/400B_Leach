ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/collision_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="Local Group simulation, 3 orthogonal views" -metadata year="2020" \
   -metadata author="Colin Leach" collisions.mp4
