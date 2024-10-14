docker run --user $(id -u):$(id -g) \
    -it --mount type=bind,source=$(pwd),target=/usr/local/src/tmdet \
    tmdet:4.0.0 \
    /bin/bash -c 'cmake -B build && make -C build tmdet; /bin/bash'
