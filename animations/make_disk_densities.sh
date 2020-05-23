ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/MW_density_early_%03d.png -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p MW_disk_60.mp4
ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/M31_density_early_%03d.png -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p M31_disk_60.mp4
ffmpeg -r 10  -start_number 0  -s 1920x1080 -i png_files/M33_density_early_%03d.png -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p M33_disk_60.mp4
