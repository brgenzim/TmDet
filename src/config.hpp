#ifndef __TMDET_CONFIG__
#define __TMDET_CONFIG__

#include "System/Environment.hpp"

#define DEFAULT_APP_NAME "tmdet"

#define DEFAULT_PDB_DATA_ROOT "/home/data/wwPDB"
#define DEFAULT_PDB_STATUS_DIR "/home/data/wwPDB/status"
#define DEFAULT_PDB_DATA_DIR "/home/data/wwPDB/data"
#define DEFAULT_PDB_ENT_DIR "/home/data/wwPDB/data/structures/divided/pdb"
#define DEFAULT_PDB_CIF_DIR "/home/data/wwPDB/data/structures/divided/updated_mmcif"

#define DEFAULT_TMDET_DATA_ROOT "/usr/local/share/tmdet/data"
#define DEFAULT_TMDET_CC_DIR "/usr/local/share/tmdet/ccd"
#define DEFAULT_TMDET_CC_FILE "/usr/local/share/tmdet/components.cif.gz"
#define DEFAULT_TMDET_CC_URL "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"

#define DEFAULT_UNITMP_SCHEMA "https://"
#define DEFAULT_UNITMP_DOMAIN "unitmp.org"

#define DEFAULT_TMDET_MIN_MEMBRANE_WIDTH 10.0
#define DEFAULT_TMDET_MAX_MEMBRANE_WIDTH 18.0
#define DEFAULT_TMDET_OPTIM_BALL_DIST 0.15

extern Tmdet::System::Environment environment;
#endif