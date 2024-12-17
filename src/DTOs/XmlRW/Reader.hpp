// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <pugixml.hpp>
#include <DTOs/XmlRW/Reader3.hpp>
#include <DTOs/XmlRW/Reader4.hpp>
#include <VOs/Protein.hpp>

/**
 * @brief namespace for tmdet xml data transfer objects
 * 
 * @namespace Tmdet
 * @namespace DTOs
 * @namespace XmlRW
 */
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