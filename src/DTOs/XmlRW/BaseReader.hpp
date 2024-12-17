// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <VOs/Membrane.hpp>
#include <VOs/Modification.hpp>
#include <VOs/TMatrix.hpp>
#include <VOs/Chain.hpp>
#include <VOs/Protein.hpp>
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
     * @brief base class for reading xml files
     */
    class BaseReader  {
        protected:

            /**
             * @brief root node of the xml document
             */
            pugi::xml_node _root;

            /**
             * @brief read transformation matrix from xml document object
             * 
             * @param node 
             * @return Tmdet::VOs::TMatrix 
             */
            virtual Tmdet::VOs::TMatrix getTMatrix(const pugi::xml_node& node) const = 0;

            /**
             * @brief Get the value of tmp attribute
             * 
             * @return true 
             * @return false 
             */
            virtual bool getTmp() const = 0;

            /**
             * @brief Get the value of code attribute
             * 
             * @return std::string 
             */
            virtual std::string getCode() const = 0;

            /**
             * @brief Get the value of CREATE_DATE node
             * 
             * @return std::string 
             */
            virtual std::string getCreateDate() const = 0;

            /**
             * @brief Get the value of Qvalue 
             * 
             * @return double 
             */
            virtual double getQvalue() const = 0;

            /**
             * @brief Get the value of transmembrane type
             * 
             * @return std::string 
             */
            virtual std::string getTmtype() const = 0;

            /**
             * @brief Get the content of MEMBRANE node
             * 
             * @return std::vector<Tmdet::VOs::Membrane> 
             */
            virtual std::vector<Tmdet::VOs::Membrane> getMembranes() const = 0;

            /**
             * @brief Get the content of CHAIN node
             * 
             * @return std::vector<Tmdet::VOs::XmlChain>
             */
            virtual std::vector<Tmdet::VOs::XmlChain> getChains() = 0;

            /**
             * @brief Get the content of REGION node
             * 
             * @param cnode 
             * @return std::vector<Tmdet::VOs::Region> 
             */
            virtual std::vector<Tmdet::VOs::Region> getRegions(const pugi::xml_node& cnode) const = 0;


        public:

            /**
             * @brief read tmdet xml file into Protein Value Object
             * 
             * @param Tmdet::VOs::Xml& xmlData
             * @param path 
             */
            virtual void readXml(Tmdet::VOs::Xml& xmlData) = 0;

            /**
             * @brief Destroy the BaseReader object
             */
            virtual ~BaseReader()=default;

            /**
             * @brief set the root element of xml document
             * 
             * @param path 
             */
            virtual void setRoot(const pugi::xml_document& doc) = 0;

    };
}
