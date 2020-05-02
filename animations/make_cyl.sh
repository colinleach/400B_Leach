#ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/cyl_MW_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="MW disk particles, cylindrical coordinates" -metadata year="2020" \
   -metadata author="Colin Leach" cyl_MW_disk.mp4
#ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/cyl_M31_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="M31 disk particles, cylindrical coordinates" -metadata year="2020" \
   -metadata author="Colin Leach" cyl_M31_disk.mp4
ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/cyl_M33_%03d.png -vcodec libx264 -vf fps=25 -crf 25\
  -pix_fmt yuv420p -metadata title="M33 disk particles, cylindrical coordinates" -metadata year="2020" \
   -metadata author="Colin Leach" cyl_M33_disk.mp4
