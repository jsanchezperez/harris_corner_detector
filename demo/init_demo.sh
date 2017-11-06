#!/bin/sh
set -e
measure=$1 # (0.Harris, 1.Shi-Tomasi, 2.Szeliski)
k_harris=$2
shi_tomasi_threshold=$3
sigma_denosing=$4
sigma_windows=$5
radius=$5
subpixel_precision=$6

k=$k_harris
if [ "$measure" == "1" ];then
        k=$shi_tomasi_threshold
fi

harris_corner_detector input_0.png -o selected_points.png -f harris_corners.txt -k $k -i $sigma_denosing -s $sigma_windows -w $radius -p 0.0001 -m $measure -q $subpixel_precision
