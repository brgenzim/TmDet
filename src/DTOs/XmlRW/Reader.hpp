#pragma once

#include <string>
#include <pugixml.hpp>
#include <DTOs/XmlRW/Reader3.hpp>
#include <DTOs/XmlRW/Reader4.hpp>
#include <VOs/Protein.hpp>

namespace Tmdet::DTOs::XmlRW {

    /**
     * @brief reader detect the version of input
     *        xml file and read accordingly
     * 
     */
    class Reader {
        private:
            
            /**
             * @brief pugixml document
             */
            pugi::xml_document doc;

            /**
             * @brief Reader3 class
             */
            Tmdet::DTOs::XmlRW::Reader3 reader3;

            /**
             * @brief Reader4 class
             */
            Tmdet::DTOs::XmlRW::Reader4 reader4;

            /**
             * @brief read the xml file into pugixml object
             * 
             * @param path 
             */
            void read(const std::string& path);

        public:

            /**
             * @brief read and parse xml file into xml value object
             * 
             * @param xml
             * @param path 
             */
            void readXml(Tmdet::VOs::Xml& xmlData, const std::string& path);

            /**
             * @brief check the version of the xml file is above 4.0
             * 
             * @return true 
             * @return false 
             */
            bool isVersion4();
    };
}