#pragma once

#include <string>
#include <vector>
#include <pugixml.hpp>
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
     * @brief base class for writing xml files
     */
    class BaseWriter  {
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
             * @brief writing transformation matrix into xml document object
             * 
             * @param node 
             * @param tmatrix 
             */
            virtual void setTMatrix(const pugi::xml_node& node, Tmdet::ValueObjects::TMatrix& tmatrix) const = 0;

            /**
             * @brief writing pugixml object into file
             * 
             * @param path 
             */
            void write(const std::string& path) const;

            /**
             * @brief create an empty pugixml object
             * 
             */
            virtual void create() = 0;

            /**
             * @brief Set the value of tmp attribute
             * 
             * @param tmp 
             */
            virtual void setTmp(const bool& tmp) const = 0;

            /**
             * @brief Set the value of code attribute
             * 
             * @param code 
             */
            virtual void setCode(const std::string& code) const = 0;

            /**
             * @brief Set the value of CREATE_DATE node
             * 
             * @param date 
             */
            virtual void setCreateDate(const std::string& date) const = 0;

            /**
             * @brief Set the content of MODIFICATIONS node
             * 
             * @param mods 
             */
            virtual void setModifications(const std::vector<Tmdet::ValueObjects::Modification>& mods) = 0;

            /**
             * @brief Set the qValue
             * 
             * @param q 
             */
            virtual void setQvalue(const double& q) const = 0;

            /**
             * @brief Set the  value of transmembrane type
             * 
             * @param type 
             */
            virtual void setTmtype(const std::string& type) const = 0;

            /**
             * @brief Set the content of MEMBRANE node
             * 
             * @param membranes 
             */
            virtual void setMembranes(std::vector<Tmdet::ValueObjects::Membrane>& membranes) = 0;

            /**
             * @brief Set the content of CHAIN node
             * 
             * @param chains 
             */
            virtual void setChains(const std::vector<Tmdet::ValueObjects::Chain>& chains) = 0;

            /**
             * @brief Set the content of REGION node
             * 
             * @param pnode 
             * @param regions 
             */
            virtual void setRegions(pugi::xml_node& pnode, const std::vector<Tmdet::ValueObjects::Region>& regions) const = 0;


        public:

            /**
             * @brief write Tmdet Value Objects out in xml
             * 
             * @param protein 
             * @param path 
             */
            virtual void writeXml(Tmdet::ValueObjects::Protein& protein, const std::string& path) = 0;

            virtual ~BaseWriter()=default;
    };
}
