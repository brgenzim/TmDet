#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <DTOs/Xml/BaseWriter.hpp>
#include <ValueObjects/BioMatrix.hpp>

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::Xml {

    /**
     * @brief class for writing old (pdbtm_3.0) xml files
     */
    class Writer3 : public BaseWriter {
        private:
            const std::string _pdbtm_xml=R"(
<pdbtm xmlns="https://pdbtm.unitmp.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="https://pdbtm.unitmp.org/data/pdbtm.xsd pdbtm.xsd" ID="xxxx" TMP="unk">
  <COPYRIGHT>
     All  information, data  and  files are copyright.  This document  is
     produced  in  the  Institute  of  Molecular  Life Sciences, Research
     Centre of Natural Sciences,  HUN-REN, Budapest, Hungary.  There  are  
     no restrictions on its use by non-profit institutions as long as its
     content is in no way modified and this statement is not removed from 
     entries. Usage  by  and  for  commercial entities and for commercial
     purpose  requires a license agreement.  It can be purchased from the
     UniTmp home page (https://www.unitmp.org/licences).
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

            /**
             * @brief Set the value of SPRES node
             * 
             * @param type 
             */
            void setSpres(const std::string& type) const;

            /**
             * @brief Set the value of PDBKWRES node
             * 
             * @param type 
             */
            void setPdbkwres(const std::string& type) const;

            /**
             * @brief Set the content of BIOMATRIX node
             * 
             * @param bioMatrix 
             */
            void setBioMatrix(Tmdet::ValueObjects::BioMatrix& bioMatrix);


        protected:
            void setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const;
            void create();
            void setTmp(const bool& tmp) const;
            void setCode(const std::string& code) const;
            void setCreateDate(const std::string& date) const;
            void setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods);
            void setQvalue(const double& q) const;
            void setTmtype(const std::string& type) const;
            void setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes);
            void setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains);
            void setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const;
            void writeXml(Tmdet::ValueObjects::Protein& protein, const std::string& path);
    };
}
