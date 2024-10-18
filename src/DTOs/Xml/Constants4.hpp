#pragma once

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::Xml {
            
    const char* XML_ATTR_AUTH_ASYM_ID="pdb_auth_asym_id";
    const char* XML_ATTR_AUTH_SEQ_ID="pdb_auth_seq_id";
    const char* XML_ATTR_DATE="date";
    const char* XML_ATTR_HALF_THICKNESS="halfThickness";
    const char* XML_ATTR_LABEL_ASYM_ID="pdb_label_asym_id";
    const char* XML_ATTR_LABEL_SEQ_ID="pdb_label_seq_id";
    const char* XML_ATTR_PDB_CODE="pdbCode";
    const char* XML_ATTR_NUM_TM="numTM";
    const char* XML_ATTR_SEQ="seq";
    const char* XML_ATTR_TRANSMEMBRANE="transmembrane";
    const char* XML_ATTR_TYPE="type";
    const char* XML_ATTR_X="x";
    const char* XML_ATTR_Y="y";
    const char* XML_ATTR_Z="z";
    
    const char* XML_NODE_ASSEMBLY="assembly";
    const char* XML_NODE_COPYRIGHT="copyright";
    const char* XML_NODE_CHAIN="chain";
    const char* XML_NODE_CREATED="created";
    const char* XML_NODE_END="end";
    const char* XML_NODE_MEMBRANE="membrane";
    const char* XML_NODE_MODIFICATION="modification";
    const char* XML_NODE_MODIFICATIONS="modifications";
    const char* XML_NODE_ROOT="pdbtm";
    const char* XML_NODE_RAWDATA="rawData";
    const char* XML_NODE_REGION="region";
    const char* XML_NODE_ROTATE="rotate";
    const char* XML_NODE_ROWX="rowX";
    const char* XML_NODE_ROWY="rowY";
    const char* XML_NODE_ROWZ="rowZ";
    const char* XML_NODE_SEQENCE="sequence";
    const char* XML_NODE_START="start";
    const char* XML_NODE_TRANSFORMATION="transformation";
    const char* XML_NODE_TRANSLATE="translate";
    const char* XML_NODE_QVALUE="qValue";
    const char* XML_NODE_TMTYPE="tmType";
    const char* XML_NODE_TMDET_VERSION="tmdetVersion";
    
}
