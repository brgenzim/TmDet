ARGS="--copyright-prefix spdx-symbol --license 'CC-BY-NC-4.0' --year 2003"
ARGS="$ARGS --copyright 'HUN-REN Research Center of Natural Sciences'"
ARGS="$ARGS --copyright 'Gabor E. Tusnady <tusnady.gabor@ttk.hu'"
ARGS="$ARGS --contributor 'Csongor Gerdan <gerdan.csongor@ttk.hu'"

FILES="src/cli/tmdet.cpp src/CMakeLists.txt"

eval pipx run reuse annotate $ARGS $FILES
