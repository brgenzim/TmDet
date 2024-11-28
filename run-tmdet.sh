WD=/home/A/csongor/dev/UniTmp/Tmdet.CPP
cd "$WD"
docker run --rm --user $(id -u):$(id -g) -t \
    --mount type=bind,source=$(pwd),target=/usr/local/src/tmdet \
    --mount type=bind,source=/zfs,target=/zfs \
    --mount type=bind,source=/home/tusi,target=/home/tusi \
    --mount type=bind,source=/home/B,target=/home/B \
    --mount type=bind,source=/home/data,target=/home/data \
    --mount type=bind,source=/data,target=/data \
    --mount type=bind,source=/home/A/csongor/dev/UniTmp/unitmp/storage/app/tmdet_job_files,target=/home/A/csongor/dev/UniTmp/unitmp/storage/app/tmdet_job_files \
    tmdet:4.0.0 \
    /usr/local/src/tmdet/build/src/tmdet $@

