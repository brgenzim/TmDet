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
    class Writer4 : public BaseWriter {
        private:
            const std::string _pdbtm_xml=R"(
<pdbtm xmlns="https://pdbtm.unitmp.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="https://pdbtm.unitmp.org/data/pdbtm.xsd pdbtm.xsd" pdbCode="xxxx" transmembrane="unk">
  <copyright>
     All  information, data  and  files are copyright.  This document  is
     produced  in  the  Institute  of  Molecular  Life Sciences, Research
     Centre of Natural Sciences,  HUN-REN, Budapest, Hungary.  There  are  
     no restrictions on its use by non-profit institutions as long as its
     content is in no way modified and this statement is not removed from 
     entries. Usage  by  and  for  commercial entities and for commercial
     purpose  requires a license agreement.  It can be purchased from the
     UniTmp home page (https://www.unitmp.org/licences).
  </copyright>
  <rawData>
    <created>yyy-mm-dd</created>
    <tmdetVersion>4.0.0</tmdetVersion>
    <assembly>1</assembly>
    <qValue>00.00</qValue>
    <tmType>xxxxx</tmType>
  </rawData>
</pdbtm>
)";

            const std::string _membrane_xml=R"(
<membrane halfThickness="10.0" type="Plain">
  <transformation>
    <translate x="0.00" y="0.00" z="0.00"/>
    <rotate>
        <rowX x="1.00000000" y="0.00000000" z="0.00000000"/>
        <rowY x="0.00000000" y="1.00000000" z="0.00000000"/>
        <rowZ x="0.00000000" y="0.00000000" z="1.00000000"/>
    </rotate>
  </transformation>
</membrane>
)";
            const std::string _chain_xml=R"(
<chain pdb_auth_asym_id="xxx" pdb_label_asym_id numTM="xxx" type="xxx">
  <sequence>
  </sequence>
</chain>
)";
            const std::string _region_xml=R"(
<region>
  <start seq="a" pdb_auth_seq_id="b" padb_label_seq_id="c"/>
  <end seq="a" pdb_auth_seq_id="b" padb_label_seq_id="c"/>
</region>
)";


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
