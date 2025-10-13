[Get and Install TmDet](#local-install) section describes setup on user's operating system.
By installing the software and configuring the necessary dependencies, the environment will be altered to support the application.

Read [Build and run Docker image from local source directory](#docker-install) section,
if you would like to avoid modifying your local environment.
We recommend running the software in a Docker container in this case.
This ensures the software runs in an isolated environment with minimal changes to your system's core settings.

# Table of Contents

* [Get and Install **TmDet**](#local-install)
* [Build and run Docker image from local source directory](#docker-install)
* [Command line arguments](#cmd-args)
* [Set environment](#variables)
* [Chemical component directory](#ccd)

<a name="local-install"></a>
# Get and Install **TmDet**

## Setup required software environment

> **NOTE**
>
> These commands were tested on Ubuntu 24.04. Other Linux distributions are out of
> scope of this document.
>
> Below commands requires ```sudo``` permission.

1. Install build tools

```
sudo apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gcc g++
```

2. Install dependencies

```
sudo apt-get install -y \
    curl \
    git \
    libcurl4-openssl-dev \
    libeigen3-dev \
    libzip-dev \
    nlohmann-json3-dev

GEMMI_VERSION=0.7.0
cd /tmp
curl -L -O https://github.com/project-gemmi/gemmi/archive/refs/tags/v$GEMMI_VERSION.tar.gz && \
    tar -xzf v$GEMMI_VERSION.tar.gz && \
    rm v$GEMMI_VERSION.tar.gz && \
    cd gemmi-$GEMMI_VERSION && \
    cmake -B build && \
    make -j4 -C build && \
    sudo make -C build install
```

```
PUGIXML_VERSION=1.14
cd /tmp && mkdir -p contrib && cd contrib
curl -L -O https://github.com/zeux/pugixml/archive/refs/tags/v$PUGIXML_VERSION.tar.gz && \
    tar -xzf v$PUGIXML_VERSION.tar.gz && \
    rm v$PUGIXML_VERSION.tar.gz && \
    ln -s pugixml-$PUGIXML_VERSION pugixml
```

## Build commands

1. Clone the TmDet repository:

```
cd /tmp
git clone https://github.com/brgenzim/TmDet.git
```

2. Change to the TmDet folder:

```
cd TmDet
```

3. Compile the TmDet:

```
cmake -B build && make -j4 -C build && sudo make -C build install
```

4. The binary is located in the ```/usr/local/bin``` folder. Enjoy it by typing ```tmdet -h```!

<a name="docker-install"></a>
# Build and run Docker image from local source directory

## Prerequisite

Current user must have ```sudo``` right or user must be member of ```docker``` user group.
In the latter case ```sudo``` can be omitted from the command lines below.

## Commands

1. Clone the TmDet repository:

```
git clone https://github.com/brgenzim/TmDet.git
```

2. Build docker image:

```
cd TmDet && sudo docker compose --env-file .env.example build
```

3. Run tmdet in docker container:

```
sudo bash run-tmdet.sh -pi ./data/1a0s.cif -po ./data/1a0s.tr.cif -x ./data/1a0s.xml
```

> **NOTE**
>
> Input file must be in the current directory where ```run-tmdet.sh``` is executed.

<a name="cmd-args"></a>
# Command line arguments
- Get help:

```
tmdet -h
```

- Set input:

    - by path:
        >-pi /path/to/the/pdb/struct[.cif|.cif.gz|ent|ent.gz]
    - by pdbCode (```PDB_ENT_DIR``` and/or ```PDB_CIF_DIR``` have to be defined in this case, see *Set environments*)
        >-c pdbCode

        - refine by assembly id (using assemblies downloaded from rcsb.org):
        >-a assembly_id (default is 1)

    - in case of NMR structures model number also can be set:
    >-m model_number

    - unselect chain(s):
    >-uc chainId1,chainId2

    - antibodies are unselected automatically, user can prevent this by 'force no delete antibodies' switch:
    >-fa

- Set output:
    - by path (if -c is not used (see input parameters)):
        >-po /path/to/the/pdb/transformed_struct.cif.gz
        >-x /path/to/the/tmdet/output.xml

- Set main operation mode:
    - Search for curved membrane:
        >-cm
    - Search for double membrane:
        >-dm
    - Search membrane plane in inaccurately modeled structure (fragment analysis):
        >-fr

- Other parameters:
    | Short | Long | Type | Description |
    |-------|------|------|-------------|
    | -nc | --no_cache| Bool | Do not use cached data (default: *false*)|
    | -ns | --no_symmetry | Bool | Do not use symmetry axes as membrane normal (default: *false*)|
    | -lq | --lower_qvalue | float | Lower qValue, above it is membrane (default: *30*)|
    | -hq | --higher_qvalue | float | Higher qValue, limit for transmembrane type (default: *36*)|
    | -hq2 | --higher_qvalue2 | float | Higher qValue2, limit for second membrane (default: *55*)|
    | -minht | --minimum_of_half_thickness | float | Minimum value of half thickness (default: *10.0*)|
    | -maxht | --maximum_of_half_thickness | float | Maximum value of half thickness (default: *20.0*)|
    | -maxcht | --maximum_of_curved_half_thickness | float | Maximum value of half thickness for curved membrane detection (default: *14*)|
    | -ihml | --ifh_hydrph_limit | float | Hydrophobicity momentum limit for ifh detection (default: *1.6*)|
    | -ias | --ifh_avg_surface | float | Average free solvent accessible surface limit for ifh detection (default: *40*)|
    | -ian | --ifh_angle | float | Maximum angle between membrane plane and ifh (default: *5*)|
    | -iml | --ifh_min_length | int | Minimum length of ifhs (default: *6*)|
    | -ba | --boost_angle | float | Boost secondary structure element angle in optimization (default: *0.7*)|
    | -bba | --boost_beta_angle | float | Boost beta sheet angle in optimization (default: *0.55*)|
    | -bp | --boost_polarity | float | Boost polarity calculation in optimization (default: *0.55*)|
    | -lmhp | --loop_min_helix_part | float | Minimum of a helix be part as re-entrant loop (default: *0.25*)|
    | -lmd | --loop_min_depth | float | Minimum depth of a re-entrant loop in angstrom (default: *3.0*)|
    | -lmnss | --loop_min_no_sec_str | int | Minimum number of residues in a re-entrant loop that has no secondary structure (default: *1*)|
    | -mums | --max_unannotated_memb_segm | int | Maximum number of residues that can not be annotated (default: *10*)|
    | -sm | --shift_membrane | float | Shift membrane with the given distance (default: *0*)|
    | -mltmh | --min_length_of_tmh | int | Minimum length of transmembrane helix (default: *12*)|
    | -spen | --straigth_penalty | float | Additional value for normalizing straigth (default: *0.0*)|
    | -mcbs | --min_contacts_between_sheets | int | Minimum of contacts between sheets for barrel detection (default: *5*)|
    | -bh | --broken_helix | string | Type of broken helix (loop (L) or transmembrane helix (H)) (default: *L*)|
    | -minbs | --min_number_of_beta_sheets | int | Minimum number of beta sheets in beta barrel (default: *8*)|
    | -bi | --barrel_inside | string | Indicate chains those are within a beta-barrel (but not part of barrel, like chain B in 5iv8)|

<a name="variables"></a>
# Set environment

Default environment file is ```.env``` in the current directory. The path of the environment file can be set by the ```-e``` command line argument or by the ```TMDET_ENV_FILE``` environment variable. The priority of environment files is as follows:

* if ```-e``` is given, this file takes precedence over ```TMDET_ENV_FILE``` variable or ```.env``` in the current directory
* if ```-e``` is omitted but ```TMDET_ENV_FILE``` variable is set, the given file will be processed
* if neither ```-e``` nor ```TMDET_ENV_FILE``` are used, TMDET attempts to read ```.env``` file in the current directory

<a name="ccd"></a>
# Chemical component directory

Chemical component directory is important during parsing CIF files. If it is missing, TmDet
downloads it from WWPDB and install it automatically.
