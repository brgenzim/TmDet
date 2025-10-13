// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <System/Arguments.hpp>
#include <VOs/Protein.hpp>
#include <VOs/Xml.hpp>

/**
 * @brief namespace for tmdet data transfer objects
 *
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief convert protein value object into xml value object
     * 
     */
    class Xml {
        private:
            /**
             * @brief tmdet xml value object
             */
            Tmdet::VOs::Xml xmlData;

            /**
             * @brief output xml format (v3 or v4, default = v4)
             */
            std::string outXmlFmt = "v4";

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
            void write(const std::string& xmlPath, const Tmdet::System::Arguments& args);

            /**
             * @brief copy protein value object to xml value object and write to file
             * 
             * @param xmlPath 
             * @param protein 
             * @param args
             */
            void write(const std::string& xmlPath, const Tmdet::VOs::Protein& protein, const Tmdet::System::Arguments& args);

            /**
             * @brief change xml file to not transmembrane protein
             * 
             * @param xmlPath 
             */
            void notTransmembrane(const std::string& xmlInputPath, const std::string& xmlOutputPath, const Tmdet::System::Arguments& args);

            /**
             * @brief Set the file path
             * 
             * @param code 
             * @param x1 
             * @param x2 
             * @return std::string 
             */
            std::string setPath(const std::string& code, const std::string& x1, const std::string& x2) const;

            /**
             * @brief set output xml format to v3
             */
            void setV3Fmt();
    };
}
