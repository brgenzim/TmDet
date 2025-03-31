# Get and Install **TmDet**

1. Clone the TmDet repository:
> git clone https://github.com/brgenzim/TmDet

2. Change to the TmDet folder:
> cd TmDet

3. Compile the TmDet:
> cmake -B build; make -C build; make -C build install

4. Enjoy it!


# Command line arguments
- Get help:
>tmdet -h

- Set input:

    - by path:
        >-pi /path/to/the/pdb/struct[.cif|.cif.gz|ent|ent.gz]
    - by pdbCode (PDB_ENT_DIR and/or PDB_CIF_DIR have to be defined in this case, see Set environments)
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

# Set environment

Default environment file is *.env* in the current directory. The peth of the environment file can be set by the *-e* command line arguments.

# Chemical component directory

Chemical component directory is important during parsing CIF files. If it is missing, TmDet
download it from WWPDB and install it automatically.
