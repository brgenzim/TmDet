// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace XmlRW
 */
namespace Tmdet::DTOs::XmlRW {
            
    inline const char* XML3_ATTR_CHAINID="CHAINID";
    inline const char* XML3_ATTR_ID="ID";
    inline const char* XML3_ATTR_NEW_CHAINID="NEW_CHAINID";
    inline const char* XML3_ATTR_NUM_TM="NUM_TM";
    inline const char* XML3_ATTR_PDB_BEG="pdb_beg";
    inline const char* XML3_ATTR_PDB_BEGI="pdb_begi";
    inline const char* XML3_ATTR_PDB_END="pdb_end";
    inline const char* XML3_ATTR_PDB_ENDI="pdb_endi";
    inline const char* XML3_ATTR_SEQ_BEG="seq_beg";
    inline const char* XML3_ATTR_SEQ_END="seq_end";
    inline const char* XML3_ATTR_TMP="TMP";
    inline const char* XML3_ATTR_TYPE="TYPE";
    inline const char* XML3_ATTR_type="type";
    inline const char* XML3_ATTR_T="T";
    inline const char* XML3_ATTR_X="X";
    inline const char* XML3_ATTR_Y="Y";
    inline const char* XML3_ATTR_Z="Z";
    
    inline const char* XML3_NODE_APPLY_TO_CHAIN="APPLY_TO_CHAIN";
    inline const char* XML3_NODE_BIOMATRIX="BIOMATRIX";
    inline const char* XML3_NODE_CHAIN="CHAIN";
    inline const char* XML3_NODE_CREATE_DATE="CREATE_DATE";
    inline const char* XML3_NODE_DATE="DATE";
    inline const char* XML3_NODE_DELETE="DELETE";
    inline const char* XML3_NODE_DESCR="DESCR";
    inline const char* XML3_NODE_MATRIX="MATRIX";
    inline const char* XML3_NODE_MEMBRANE="MEMBRANE";
    inline const char* XML3_NODE_MODIFICATION="MODIFICATION";
    inline const char* XML3_NODE_NORMAL="NORMAL";
    inline const char* XML3_NODE_ROOT="pdbtm";
    inline const char* XML3_NODE_PDBKWORD="PDBKWORD";
    inline const char* XML3_NODE_PDBKWRES="PDBKWRES";
    inline const char* XML3_NODE_RAWRES="RAWRES";
    inline const char* XML3_NODE_REGION="REGION";
    inline const char* XML3_NODE_ROWX="ROWX";
    inline const char* XML3_NODE_ROWY="ROWY";
    inline const char* XML3_NODE_ROWZ="ROWZ";
    inline const char* XML3_NODE_SEQ="SEQ";
    inline const char* XML3_NODE_SPRES="SPRES";
    inline const char* XML3_NODE_TMATRIX="TMATRIX";
    inline const char* XML3_NODE_TMRES="TMRES";
    inline const char* XML3_NODE_TMTYPE="TMTYPE";
    inline const char* XML3_NODE_TMDET_VERSION="TMDET_VERSION";
    
}
