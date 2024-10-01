#ifndef __TMDET_DTOS_XML__
#define __TMDET_DTOS_XML__

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>

namespace Tmdet::DTOs {

    class Xml {
        private:
            pugi::xml_document _doc;
            pugi::xml_node _root;
            Tmdet::ValueObjects::TMatrix getTMatrix(const pugi::xml_node& node) const;
            void setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const;

            const std::string _pdbtm_xml=R"(
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
            const std::string _biomatrix_xml=R"(
<MATRIX ID="xxx">
  <APPLY_TO_CHAIN CHAINID="xxx" NEW_CHAINID="xxx"/>
  <TMATRIX>
    <ROWX X="1.00000000" Y="0.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWY X="0.00000000" Y="1.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWZ X="0.00000000" Y="0.00000000" Z="1.00000000" T="0.00000000"/>
  </TMATRIX>
</MATRIX>
)";

            const std::string _membrane_xml=R"(
<MEMBRANE>
  <NORMAL X="0.0" Y="0.0" Z="15.0"/>
  <TMATRIX>
    <ROWX X="1.00000000" Y="0.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWY X="0.00000000" Y="1.00000000" Z="0.00000000" T="0.00000000"/>
    <ROWZ X="0.00000000" Y="0.00000000" Z="1.00000000" T="0.00000000"/>
  </TMATRIX>
</MEMBRANE>
)";
            const std::string _chain_xml=R"(
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
            
            void read(const std::string& path);
            void write(const std::string& path) const;
            void create();
            bool getTmp() const;
            void setTmp(const bool& tmp) const;
            std::string getCode() const;
            void setCode(const std::string& code) const;
            std::string getCreateDate() const;
            void setCreateDate(const std::string& date) const;
            std::vector<Tmdet::ValueObjects::Modification> getModifications() const;
            void setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods);
            double getQvalue() const;
            void setQvalue(const double& q) const;
            std::string getTmtype() const;
            void setTmtype(const std::string& type) const;
            std::string getSpres() const;
            void setSpres(const std::string& type) const;
            std::string getPdbkwres() const;
            void setPdbkwres(const std::string& type) const;
            Tmdet::ValueObjects::BioMatrix getBioMatrix() const;
            void setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix);
            std::vector<Tmdet::ValueObjects::Membrane> getMembranes() const;
            void setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes);
            void getChains(std::vector<Tmdet::ValueObjects::Chain>& chains);
            void setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains);
            std::vector<Tmdet::ValueObjects::Region> getRegions(const pugi::xml_node& cnode) const;
            void setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const;
    };
}

#endif