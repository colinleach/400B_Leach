ffmpeg -r 10  -start_number 430  -s 1920x1080 -i png_files/remnant_stellar_disp_y_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="M33 disk particles, cylindrical coordinates" -metadata year="2020" \
   -metadata author="Colin Leach" remnant_stellar_disp_y.mp4

ffmpeg -r 10  -start_number 430  -s 1920x1080 -i png_files/remnant_stellar_disp_z_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="M33 disk particles, cylindrical coordinates" -metadata year="2020" \
   -metadata author="Colin Leach" remnant_stellar_disp_z.mp4
