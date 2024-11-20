#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Protein.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Xml.hpp>

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::XmlRW {

    /**
     * @brief class for writing old (pdbtm_3.0) xml files
     */
    class Writer {
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
            /**
             * @brief pugixml document
             */
            pugi::xml_document _doc;

            /**
             * @brief root node of the xml document
             */
            pugi::xml_node _root;
        
            /**
             * @brief create an empty pugixml object
             * 
             */
            void create();

            /**
             * @brief write xml objects out
             */
            void write(const std::string& path) const;
            
            /**
             * @brief Set the value of tmp attribute
             * 
             * @param tmp 
             */
            void setTmp(const bool& tmp) const;

            /**
             * @brief Set the value of code attribute
             * 
             * @param code 
             */
            void setCode(const std::string& code) const;

            /**
             * @brief Set the value of CREATE_DATE node
             * 
             * @param date 
             */
            void setCreateDate(const std::string& date) const;

			      /**
			       * @brief Set the value of tmdet version
			       * 
			       * @param version 
			       */
            void setVersion(const std::string& version) const;

            /**
             * @brief Set the qValue
             * 
             * @param q 
             */
            void setQvalue(const double& q) const;

            /**
             * @brief Set the  value of transmembrane type
             * 
             * @param type 
             */
            void setTmtype(const std::string& type) const;

            /**
             * @brief set transformation
             * 
             * @param tmatrix 
             */
            void setTMatrix(Tmdet::ValueObjects::TMatrix& tmatrix);

            /**
             * @brief Set the content of MEMBRANE node
             * 
             * @param membranes 
             */
            void setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes);

            /**
             * @brief Set the content of CHAIN node
             * 
             * @param chains 
             */
            void setChains(const std::vector<Tmdet::ValueObjects::XmlChain>& chains);

            /**
             * @brief Set the content of REGION node
             * 
             * @param pnode 
             * @param regions 
             */
            void setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions, const std::vector<Tmdet::ValueObjects::Residue>& residues) const;
          public:
            void writeXml(Tmdet::ValueObjects::Xml& xmlData, const std::string& path);
    };
}
