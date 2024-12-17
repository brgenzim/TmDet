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

            /**
             * @brief Get the content of transformation matrix
             * 
             * @param node 
             * @return Tmdet::VOs::TMatrix 
             */
            Tmdet::VOs::TMatrix getTMatrix(const pugi::xml_node& node) const final;

            /**
             * @brief Get the value of tmp attribute
             * 
             * @return true 
             * @return false 
             */
            bool getTmp() const;

            /**
             * @brief Get the Code of the protein
             * 
             * @return std::string 
             */
            std::string getCode() const;

            /**
             * @brief Get create date
             * 
             * @return std::string 
             */
            std::string getCreateDate() const;

            /**
             * @brief Get qValue
             * 
             * @return double 
             */
            double getQvalue() const;

            /**
             * @brief Get transmembrane type
             * 
             * @return std::string 
             */
            std::string getTmtype() const;

            /**
             * @brief Get membrane definitions
             * 
             * @return std::vector<Tmdet::VOs::Membrane> 
             */
            std::vector<Tmdet::VOs::Membrane> getMembranes() const;

            /**
             * @brief Get chain data (sequence, type, numtm, regions)
             * 
             * @return std::vector<Tmdet::VOs::XmlChain> 
             */
            std::vector<Tmdet::VOs::XmlChain> getChains();

            /**
             * @brief Get region definitions for chain
             * 
             * @param cnode 
             * @return std::vector<Tmdet::VOs::Region> 
             */
            std::vector<Tmdet::VOs::Region> getRegions(const pugi::xml_node& cnode) const;

        public:

            /**
             * @brief read the entire xml document and parse to protein value object
             * 
             * @param xmlData 
             */
            void readXml(Tmdet::VOs::Xml& xmlData);

            /**
             * @brief Set the root of xml document
             * 
             * @param doc 
             */
            void setRoot(const pugi::xml_document& doc);
            
    };
}
