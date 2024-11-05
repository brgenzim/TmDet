#pragma once

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::XmlRW {
            
    inline const char* XML_ATTR_AUTH_ASYM_ID="pdbAuthAsymId";
    inline const char* XML_ATTR_AUTH_SEQ_ID="pdbAuthSeqId";
    inline const char* XML_ATTR_AUTH_SEQ_ICODE="pdbAuthSeqIcode";
    inline const char* XML_ATTR_DATE="date";
    inline const char* XML_ATTR_HALF_THICKNESS="halfThickness";
    inline const char* XML_ATTR_LABEL_ASYM_ID="pdbLabelAsymId";
    inline const char* XML_ATTR_LABEL_SEQ_ID="pdbLabelSeqId";
    inline const char* XML_ATTR_PDB_CODE="pdbCode";
    inline const char* XML_ATTR_NUM_TM="numTM";
    inline const char* XML_ATTR_SEQ="seq";
    inline const char* XML_ATTR_SIZE="size";
    inline const char* XML_ATTR_SPHERE_RADIUS="sphereRadius";
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
    inline const char* XML_NODE_END="end";
    inline const char* XML_NODE_MEMBRANE="membrane";
    inline const char* XML_NODE_MEMBRANES="membranes";
    inline const char* XML_NODE_MODIFICATION="modification";
    inline const char* XML_NODE_MODIFICATIONS="modifications";
    inline const char* XML_NODE_ROOT="pdbtm";
    inline const char* XML_NODE_RAWDATA="rawData";
    inline const char* XML_NODE_REGION="region";
    inline const char* XML_NODE_ROTATE="rotate";
    inline const char* XML_NODE_ROWX="rowX";
    inline const char* XML_NODE_ROWY="rowY";
    inline const char* XML_NODE_ROWZ="rowZ";
    inline const char* XML_NODE_SEQENCE="sequence";
    inline const char* XML_NODE_START="start";
    inline const char* XML_NODE_TRANSFORMATION="transformation";
    inline const char* XML_NODE_TRANSLATE="translate";
    inline const char* XML_NODE_QVALUE="qValue";
    inline const char* XML_NODE_TMTYPE="tmType";
    inline const char* XML_NODE_TMDET_VERSION="tmdetVersion";
    
}
