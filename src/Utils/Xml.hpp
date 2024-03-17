#ifndef __TMDET_UTILS_XML__
#define __TMDET_UTILS_XML__

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>

using namespace std;
using namespace gemmi;

namespace Tmdet::Utils {

    class Xml {
        private:
            pugi::xml_document _doc;
            pugi::xml_node _root;
            Tmdet::ValueObjects::TMatrix getTMatrix(pugi::xml_node& node);
            void setTMatrix(pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix);

            const string _pdbtm_xml=R"(
<pdbtm xmlns="https://pdbtm.unitmp.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="https://pdbtm.unitmp.org/data/pdbtm.xsd pdbtm.xsd" ID="xxxx" TMP="unk">
  <COPYRIGHT>
     All  information, data  and  files are copyright.  PDBTM database is
     produced  in  the  Institute  of  Molecular Life Sciences, Budapest, 
     Hungary.  There  are  no  restrictions  on  its  use  by  non-profit 
     institutions as long as its  content is  in no way modified and this 
     statement is not removed from entries.  Usage  by and for commercial 
     entities  requires  a  license  agreement (send an email to pdbtm at 
     enzim dot hu).
  </COPYRIGHT>
  <CREATE_DATE>YYYY-MM-DD</CREATE_DATE>
  <RAWRES>
    <TMRES>00.00</TMRES>
    <TMTYPE>xxxxx</TMTYPE>
    <SPRES>xxxx</SPRES>
    <PDBKWRES>xxx</PDBKWRES>
  </RAWRES>
</pdbtm>
)";
            const string _biomatrix_xml=R"(
<MATRIX ID="xxx">
  <APPLY_TO_CHAIN CHAINID="xxx" NEW_CHAINID="xxx"/>
  <TMATRIX>
    <ROWX X="1.00000000" Y="0.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWY X="0.00000000" Y="1.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWZ X="0.00000000" Y="0.00000000" Z="1.00000000" T="0.00000000"/>
  </TMATRIX>
</MATRIX>
)";

            const string _membrane_xml=R"(
<MEMBRANE>
  <NORMAL X="0.0" Y="0.0" Z="15.0"/>
  <TMATRIX>
    <ROWX X="1.00000000" Y="0.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWY X="0.00000000" Y="1.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWZ X="0.00000000" Y="0.00000000" Z="1.00000000" T="0.00000000"/>
  </TMATRIX>
</MEMBRANE>
)";
            const string _chain_xml=R"(
<CHAIN CHAINID="xxx" NUM_TM="xxx" TYPE="xxx">
  <SEQ>
  </SEQ>
</CHAIN>
)";
            const char* XML_ATTR_CHAINID="CHAINID";
            const char* XML_ATTR_ID="ID";
            const char* XML_ATTR_NEW_CHAINID="NEW_CHAINID";
            const char* XML_ATTR_NUM_TM="NUM_TM";
            const char* XML_ATTR_PDB_BEG="pdb_beg";
            const char* XML_ATTR_PDB_END="pdb_end";
            const char* XML_ATTR_SEQ_BEG="seq_beg";
            const char* XML_ATTR_SEQ_END="seq_end";
            const char* XML_ATTR_TMP="TMP";
            const char* XML_ATTR_TYPE="TYPE";
            const char* XML_ATTR_type="type";
            const char* XML_ATTR_T="T";
            const char* XML_ATTR_X="X";
            const char* XML_ATTR_Y="Y";
            const char* XML_ATTR_Z="Z";
            
            const char* XML_NODE_APPLY_TO_CHAIN="APPLY_TO_CHAIN";
            const char* XML_NODE_BIOMATRIX="BIOMATRIX";
            const char* XML_NODE_CHAIN="CHAIN";
            const char* XML_NODE_CREATE_DATE="CREATE_DATE";
            const char* XML_NODE_DATE="DATE";
            const char* XML_NODE_DELETE="DELETE";
            const char* XML_NODE_DESCR="DESCR";
            const char* XML_NODE_MATRIX="MATRIX";
            const char* XML_NODE_MEMBRANE="MEMBRANE";
            const char* XML_NODE_MODIFICATION="MODIFICATION";
            const char* XML_NODE_NORMAL="NORMAL";
            const char* XML_NODE_ROOT="pdbtm";
            const char* XML_NODE_PDBKWORD="PDBKWORD";
            const char* XML_NODE_PDBKWRES="PDBKWRES";
            const char* XML_NODE_RAWRES="RAWRES";
            const char* XML_NODE_REGION="REGION";
            const char* XML_NODE_ROWX="ROWX";
            const char* XML_NODE_ROWY="ROWY";
            const char* XML_NODE_ROWZ="ROWZ";
            const char* XML_NODE_SEQ="SEQ";
            const char* XML_NODE_SPRES="SPRES";
            const char* XML_NODE_TMATRIX="TMATRIX";
            const char* XML_NODE_TMRES="TMRES";
            const char* XML_NODE_TMTYPE="TMTYPE";

        public:
            
            void read(string& path);
            void write(string& path);
            void create();
            bool getTmp();
            void setTmp(bool& tmp);
            string getCode();
            void setCode(string& code);
            string getCreateDate();
            void setCreateDate(string& date);
            vector<Tmdet::ValueObjects::Modification> getModifications();
            void setModifications(vector<Tmdet::ValueObjects::Modification>& mods);
            double getQvalue();
            void setQvalue(double& q);
            string getTmtype();
            void setTmtype(string& type);
            string getSpres();
            void setSpres(string& type);
            string getPdbkwres();
            void setPdbkwres(string& type);
            Tmdet::ValueObjects::BioMatrix getBioMatrix();
            void setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix);
            vector<Tmdet::ValueObjects::Membrane> getMembranes();
            void setMembranes(vector<Tmdet::ValueObjects::Membrane>& membranes);
            vector<Tmdet::ValueObjects::Chain> getChains();
            void setChains(vector<Tmdet::ValueObjects::Chain>& chains);
            vector<Tmdet::ValueObjects::Region> getRegions(pugi::xml_node& cnode);
            void setRegions(pugi::xml_node& pnode, vector<Tmdet::ValueObjects::Region>& regions);
    
    };
}

#endif