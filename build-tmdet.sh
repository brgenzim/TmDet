docker run --rm --user $(id -u):$(id -g) -it \
    --mount type=bind,source=$(pwd),target=/usr/local/src/tmdet \
    --mount type=bind,source=/zfs,target=/zfs \
    --mount type=bind,source=/home/tusi,target=/home/tusi \
    --mount type=bind,source=/home/B,target=/home/B \
    --mount type=bind,source=/home/data,target=/home/data \
    --mount type=bind,source=/data,target=/data \
    tmdet:4.0.0 \
    /bin/bash -c "$1"
    #/bin/bash -c 'cmake -B build && make -C build tmdet; /bin/bash'
