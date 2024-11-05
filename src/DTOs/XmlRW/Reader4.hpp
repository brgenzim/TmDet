#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <DTOs/XmlRW/BaseReader.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Xml.hpp>

namespace Tmdet::DTOs::XmlRW {

    /**
     * @brief class for reading new (version > 4.0) xml file
     * 
     */
    class Reader4 : public BaseReader {

        protected:
            Tmdet::ValueObjects::TMatrix getTMatrix(const pugi::xml_node& node) const;
            bool getTmp() const;
            std::string getCode() const;
            std::string getCreateDate() const;
            std::string getVersion() const;
            double getQvalue() const;
            std::string getTmtype() const;
            std::vector<Tmdet::ValueObjects::Membrane> getMembranes() const;
            std::vector<Tmdet::ValueObjects::XmlChain> getChains();
            std::vector<Tmdet::ValueObjects::Region> getRegions(const pugi::xml_node& cnode) const;
            

        public:

            void readXml(Tmdet::ValueObjects::Xml& xmlData);
            void setRoot(const pugi::xml_document& doc);
            
    };
}