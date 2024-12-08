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

namespace Tmdet::DTOs::XmlRW {

    /**
     * @brief class for reading new (version > 4.0) xml file
     * 
     */
    class Reader4 : public BaseReader {

        protected:
            Tmdet::VOs::TMatrix getTMatrix(const pugi::xml_node& node) const;
            bool getTmp() const;
            std::string getCode() const;
            std::string getCreateDate() const;
            std::string getVersion() const;
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