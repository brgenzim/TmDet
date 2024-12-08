#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <DTOs/XmlRW/BaseReader.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/BioMatrix.hpp>
#include <VOs/Modification.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/Xml.hpp>

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::XmlRW {

    /**
     * @brief class for reading old (pdbtm_3.0) xml files
     */
    class Reader3 : public BaseReader {
        private:
            /**
             * @brief Get the value of SPRES node
             * 
             * @return std::string 
             */
            std::string getSpres() const;

            /**
             * @brief Get the value of PDBKWRES node
             * 
             * @return std::string 
             */
            std::string getPdbkwres() const;

            /**
             * @brief Get the content of BIOMATRIX node
             * 
             * @return Tmdet::VOs::BioMatrix 
             */
            Tmdet::VOs::BioMatrix getBioMatrix() const;

            /**
             * @brief Get the content of MODIFICATIONS node
             * 
             * @return std::vector<Tmdet::VOs::Modification> 
             */
            std::vector<Tmdet::VOs::Modification> getModifications() const;


        protected:
            Tmdet::VOs::TMatrix getTMatrix(const pugi::xml_node& node) const final;
            bool getTmp() const;
            std::string getCode() const;
            std::string getCreateDate() const;
            double getQvalue() const;
            std::string getTmtype() const;
            std::vector<Tmdet::VOs::Membrane> getMembranes() const;
            std::vector<Tmdet::VOs::XmlChain> getChains();
            std::vector<Tmdet::VOs::Region> getRegions(const pugi::xml_node& cnode) const;

        public:

            void readXml(Tmdet::VOs::Xml& xmlData);
            void setRoot(const pugi::xml_document& doc);
            
    };
}
