#pragma once

#include <string>
#include <VOs/Protein.hpp>
#include <VOs/Xml.hpp>

namespace Tmdet::DTOs {

    class Xml {
        private:
            /**
             * @brief tmdet xml value object
             */
            Tmdet::VOs::Xml xmlData;

            /**
             * @brief copy data from protein value object to xml value object
             * 
             * @param protein 
             */
            void fromProtein(const Tmdet::VOs::Protein& protein);

            /**
             * @brief copy data from xml value object to protein value object
             * 
             * @param protein 
             */
            void toProtein(Tmdet::VOs::Protein& protein) const;

        public:

            /**
             * @brief read xml file
             * 
             * @param xmlPath 
             */
            bool read(const std::string& xmlPath);

            /**
             * @brief read xml file and copy data into protein value object
             * 
             * @param xmlPath 
             * @param protein 
             */
            void read(const std::string& xmlPath, Tmdet::VOs::Protein& protein);

            /**
             * @brief write xml value object to file
             * 
             * @param xmlPath 
             */
            void write(const std::string& xmlPath);

            /**
             * @brief copy protein value object to xml value object and write to file
             * 
             * @param xmlPath 
             * @param protein 
             */
            void write(const std::string& xmlPath, const Tmdet::VOs::Protein& protein);

            /**
             * @brief change xml file to not transmembrane protein
             * 
             * @param xmlPath 
             */
            void notTransmembrane(const std::string& xmlInputPath, const std::string& xmlOutputPath);

            /**
             * @brief Set the file path
             * 
             * @param code 
             * @param x1 
             * @param x2 
             * @return std::string 
             */
            std::string setPath(const std::string& code, const std::string& x1, const std::string& x2) const;
    };
}