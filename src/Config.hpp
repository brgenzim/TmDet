#pragma once

#define DEFAULT_APP_NAME "tmdet"

#define DEFAULT_PDB_DATA_ROOT "/home/data/wwPDB"
#define DEFAULT_PDB_STATUS_DIR "/home/data/wwPDB/status"
#define DEFAULT_PDB_DATA_DIR "/home/data/wwPDB/data"
#define DEFAULT_PDB_ENT_DIR "/home/data/wwPDB/data/structures/divided/pdb"
#define DEFAULT_PDB_CIF_DIR "/home/data/wwPDB/data/assemblies/divided/mmCIF"
#define DEFAULT_PDB_CIF_EXT ".cif.gz"

#define DEFAULT_TMDET_DATA_ROOT "/tmp/TmDet/data"
#define DEFAULT_TMDET_CACHE_ROOT "/tmp/TmDet/cache"
#define DEFAULT_TMDET_LOG_DIR "/tmp/TmDet/log"
#define DEFAULT_TMDET_LOG_FILE "/tmp/TmDet/log/tmdet.log"
#define DEFAULT_TMDET_CC_DIR "/tmp/TmDet/data/ccd"
#define DEFAULT_TMDET_CC_FILE "/tmp/TmDet/data/components.cif.gz"
#define DEFAULT_TMDET_CC_URL "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
#define DEFAULT_TMDET_POLYMER_FILTER_FILE "/tmp/TmDet/data/polymer_filter.txt"

#define DEFAULT_PDBTM_DATA_ROOT "/tmp/PDBTM"

#define DEFAULT_UNITMP_SCHEMA "https://"
#define DEFAULT_UNITMP_DOMAIN "unitmp.org"

#define DEFAULT_TMDET_MIN_NUMBER_OF_RESIDUES_IN_CHAIN "15"
#define DEFAULT_TMDET_BALL_DIST "0.15"
#define DEFAULT_TMDET_SURF_PROBSIZE "1.4"
#define DEFAULT_TMDET_SURF_ZSLICE "0.05"
#define DEFAULT_TMDET_SURF_DIST "0.1"
#define TMDET_TINY 1e-10
#define TMDET_MINIMUM_QVALUE 48
#define TMDET_MEMBRANE_QVALUE 40
#define TMDET_MEMBRANE_MIN_HALFTHICKNESS 6.5
#define TMDET_MEMBRANE_MAX_HALFTHICKNESS 15
#define TMDET_SECSTRVEC_MERGE_DIST 6.0
#define TMDET_SECSTRVEC_MERGE_ANGLE 15.0

#ifndef TMDET_LOG_LEVEL
#define TMDET_LOG_LEVEL "off"
#endif


#include <System/Environment.hpp>
#include <System/Logger.hpp>

extern Tmdet::System::Environment environment;
extern Tmdet::System::Logger logger;

namespace Tmdet {
    static std::vector<std::string> ANTIBODY_NAMES = {
        "ANTIBODY",
        "NANOBODY",
        "MEGABODY",
        "MONOBODY",
        "PROMACROBODY",
        "SYBODY",
        "FAB",
    //    "LIGHT CHAIN",
    //    "HEAVY CHAIN",
    };
}
