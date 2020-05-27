# Set linestyle 1 to blue (#0060ad)
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5

# Set linestyle 2 to red (#ff0000)
set style line 2 \
    linecolor rgb '#ff0000' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 0.5

nis_line = 5.991

plot 'nis_log.dat' with linespoints linestyle 1
replot nis_line with linespoints linestyle 2