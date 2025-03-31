// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

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
 * @namespace XmlRW
 */
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