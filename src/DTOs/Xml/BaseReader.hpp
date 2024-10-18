#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
#include <ValueObjects/Membrane.hpp>
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
     * @brief base class for reading xml files
     */
    class BaseReader  {
        protected:

            /**
             * @brief pugixml document
             */
            pugi::xml_document _doc;

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
             * @brief read entire xml file into pugixml doucument object
             * 
             * @param path 
             */
            void read(const std::string& path);

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
             * @brief Get the content of MODIFICATIONS node
             * 
             * @return std::vector<Tmdet::ValueObjects::Modification> 
             */
            virtual std::vector<Tmdet::ValueObjects::Modification> getModifications() const = 0;

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
             * @param chains 
             */
            virtual void getChains(std::vector<Tmdet::ValueObjects::Chain>& chains) = 0;

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
             * @param Tmdet::ValueObjects::Protein
             * @param path 
             */
            virtual void readXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) = 0;

            virtual ~BaseReader()=default;
    };
}
