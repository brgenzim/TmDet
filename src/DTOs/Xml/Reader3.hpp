#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <DTOs/Xml/BaseReader.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Protein.hpp>

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace Xml
 */
namespace Tmdet::DTOs::Xml {

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
             * @return Tmdet::ValueObjects::BioMatrix 
             */
            Tmdet::ValueObjects::BioMatrix getBioMatrix() const;

        protected:
            void read(const std::string& path);
            Tmdet::ValueObjects::TMatrix getTMatrix(const pugi::xml_node& node) const;
            bool getTmp() const;
            std::string getCode() const;
            std::string getCreateDate() const;
            std::vector<Tmdet::ValueObjects::Modification> getModifications() const;
            double getQvalue() const;
            std::string getTmtype() const;
            std::vector<Tmdet::ValueObjects::Membrane> getMembranes() const;
            void getChains(std::vector<Tmdet::ValueObjects::Chain>& chains);
            std::vector<Tmdet::ValueObjects::Region> getRegions(const pugi::xml_node& cnode) const;

        public:

            void readXml(Tmdet::ValueObjects::Protein& protein, const std::string& path);
    };
}
