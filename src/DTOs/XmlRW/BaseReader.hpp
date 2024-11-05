#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/TMatrix.hpp>
#include <ValueObjects/Chain.hpp>
#include <ValueObjects/Protein.hpp>
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
             * @return Tmdet::ValueObjects::TMatrix 
             */
            virtual Tmdet::ValueObjects::TMatrix getTMatrix(const pugi::xml_node& node) const = 0;

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
             * @return std::vector<Tmdet::ValueObjects::Membrane> 
             */
            virtual std::vector<Tmdet::ValueObjects::Membrane> getMembranes() const = 0;

            /**
             * @brief Get the content of CHAIN node
             * 
             * @return std::vector<Tmdet::ValueObjects::XmlChain>
             */
            virtual std::vector<Tmdet::ValueObjects::XmlChain> getChains() = 0;

            /**
             * @brief Get the content of REGION node
             * 
             * @param cnode 
             * @return std::vector<Tmdet::ValueObjects::Region> 
             */
            virtual std::vector<Tmdet::ValueObjects::Region> getRegions(const pugi::xml_node& cnode) const = 0;


        public:

            /**
             * @brief read tmdet xml file into Protein Value Object
             * 
             * @param Tmdet::ValueObjects::Xml& xmlData
             * @param path 
             */
            virtual void readXml(Tmdet::ValueObjects::Xml& xmlData) = 0;

            virtual ~BaseReader()=default;

            /**
             * @brief set the root element of xml document
             * 
             * @param path 
             */
            virtual void setRoot(const pugi::xml_document& doc) = 0;

            
    };
}
