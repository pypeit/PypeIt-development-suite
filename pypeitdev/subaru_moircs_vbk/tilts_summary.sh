#!/bin/sh
# Per-slit tilt summary from a run_pypeit log: slit, lines found,
# lines used in the final fit, fit RMS (px).  Dev script (JXP and Claude).
sed 's/\x1b\[[0-9;]*m//g' "$1" | awk '
/Computing tilts for slit\/order/ {slit=$(NF-1); split($NF,a,"[(/)]")}
/Modeling arc line tilts with/ {nl=$(NF-2)}
/Number of usable arc lines for tilts/ {use=$NF}
/RMS \(pixels\)/ {printf "%5s found=%3s used=%6s rms=%.3f\n", slit, nl, use, $NF; nl="-"}
/Did not recover any lines for slit/ {printf "%5s NO LINES\n", slit}'
