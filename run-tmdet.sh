docker run --rm --user $(id -u):$(id -g) \
     -v "$PWD:/work" \
    "brgenzim/tmdet:4.1.0" $@

