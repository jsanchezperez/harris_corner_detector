gnuplot -e "dir='$1'" graphics.plg

mkdir -p $1

ps2pdf subpixel_accuracy.eps $1/subpixel_accuracy_$1.pdf
ps2pdf subpixel_quadratic.eps $1/subpixel_quadratic_$1.pdf
ps2pdf subpixel_quartic.eps $1/subpixel_quartic_$1.pdf
ps2pdf affine.eps $1/affine_$1.pdf
ps2pdf gaussian_comparison.eps $1/gaussian_comparison_$1.pdf
ps2pdf gaussian_xi_std.eps $1/gaussian_xi_std_$1.pdf
ps2pdf gaussian_xi_fast.eps $1/gaussian_xi_fast_$1.pdf
ps2pdf gradient_comparison.eps $1/gradient_comparison_$1.pdf
ps2pdf rotation.eps $1/rotation_$1.pdf
ps2pdf scale.eps $1/scale_$1.pdf
ps2pdf noise.eps $1/noise_$1.pdf
ps2pdf illumination.eps $1/illumination_$1.pdf

rm *.eps
