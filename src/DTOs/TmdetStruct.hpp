#ifndef __TMDET_DTOS_TMDET__
#define __TMDET_DTOS_TMDET__

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

namespace Tmdet::DTOs {

    struct TmdetStruct {

        /**
         * @brief read tmdet xml file into TmdetStruct Value Object
         * 
         * @param tmdetVO 
         * @param path 
         */
        static void readXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path);

        /**
         * @brief write TmdetStruct Value Objects out in xml
         * 
         * @param tmdetVO 
         * @param path 
         */
        static void writeXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path);

        /**
         * @brief write transformed structure into file in cif format
         * 
         * @param tmdetVO 
         * @param path 
         */
        static void writeCif(const Tmdet::ValueObjects::TmdetStruct& tmdetVO, const std::string& path);

        /**
         * @brief parse pdb structure into TmdetStruct Value Object
         * 
         * @param tmdetVO 
         */
        static void parse(Tmdet::ValueObjects::TmdetStruct& tmdetVO);

        /**
         * @brief print the TmdetStruct Value Object into an out stream
         * 
         * @param os 
         * @param tmdetVO 
         */
        static void print(std::ostream& os, const Tmdet::ValueObjects::TmdetStruct& tmdetVO);
        
        /**
         * @brief print the Tmdet Chain Value Object into an out stream
         * 
         * @param os 
         * @param chainVO 
         */
        static void printChain(std::ostream& os, const Tmdet::ValueObjects::Chain& chainVO);
        
        /**
         * @brief print the Tmdet Residue Value Object into an out stream
         * 
         * @param os 
         * @param residueVO 
         */
        static void printResidue(std::ostream& os, const Tmdet::ValueObjects::Residue& residueVO);
        
        /**
         * @brief print the Tmdet Atom Value Object into an out stream
         * 
         * @param os 
         * @param atomVO 
         */
        static void printAtom(std::ostream& os, const Tmdet::ValueObjects::Atom& atomVO);

        /**
         * @brief Get the amino acid sequence of a chain
         * 
         * @param tmdetVO 
         * @param chainVO 
         * @return std::vector<std::string> 
         */
        static std::vector<std::string> getChainSequence(const Tmdet::ValueObjects::TmdetStruct& tmdetVO,
            const gemmi::Chain& chainVO);
    };
}

#endif
