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
            
    inline const char* XML_ATTR_AUTH_ID="authId";
    inline const char* XML_ATTR_END_AUTH_ID="endAuthId";
    inline const char* XML_ATTR_END_AUTH_ICODE="endAuthIcode";
    inline const char* XML_ATTR_END_LABEL_ID="endLabelId";
    inline const char* XML_ATTR_DATE="date";
    inline const char* XML_ATTR_HALF_THICKNESS="halfThickness";
    inline const char* XML_ATTR_LABEL_ID="labelId";
    inline const char* XML_ATTR_PDB_CODE="pdbCode";
    inline const char* XML_ATTR_NUM_TM="numTM";
    inline const char* XML_ATTR_SEQ="seq";
    inline const char* XML_ATTR_SIZE="size";
    inline const char* XML_ATTR_SPHERE_RADIUS="sphereRadius";
    inline const char* XML_ATTR_START_AUTH_ID="startAuthId";
    inline const char* XML_ATTR_START_AUTH_ICODE="startAuthIcode";
    inline const char* XML_ATTR_START_LABEL_ID="startLabelId";
    inline const char* XML_ATTR_TRANSMEMBRANE="transmembrane";
    inline const char* XML_ATTR_TYPE="type";
    inline const char* XML_ATTR_X="x";
    inline const char* XML_ATTR_Y="y";
    inline const char* XML_ATTR_Z="z";
    
    inline const char* XML_NODE_ASSEMBLY="assembly";
    inline const char* XML_NODE_COPYRIGHT="copyright";
    inline const char* XML_NODE_CHAIN="chain";
    inline const char* XML_NODE_CHAINS="chains";
    inline const char* XML_NODE_CREATED="created";
    inline const char* XML_NODE_MEMBRANE="membrane";
    inline const char* XML_NODE_MEMBRANES="membranes";
    inline const char* XML_NODE_MODIFICATION="modification";
    inline const char* XML_NODE_MODIFICATIONS="modifications";
    inline const char* XML_NODE_ROOT="pdbtm";
    inline const char* XML_NODE_RAWDATA="rawData";
    inline const char* XML_NODE_REGION="region";
    inline const char* XML_NODE_REGIONS="regions";
    inline const char* XML_NODE_ROTATE="rotate";
    inline const char* XML_NODE_ROWX="rowX";
    inline const char* XML_NODE_ROWY="rowY";
    inline const char* XML_NODE_ROWZ="rowZ";
    inline const char* XML_NODE_SEQENCE="sequence";
    inline const char* XML_NODE_TRANSFORMATION="transformation";
    inline const char* XML_NODE_TRANSLATE="translate";
    inline const char* XML_NODE_QVALUE="qValue";
    inline const char* XML_NODE_TMTYPE="tmType";
    inline const char* XML_NODE_TMDET_VERSION="tmdetVersion";
    
}
