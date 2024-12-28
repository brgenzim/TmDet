// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of tmdet types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    /**
     * @brief definition os secondary structure type
     */
    struct SecStruct {
        /**
         * @brief name of the secondary structure type
         */
        std::string name;

        /**
         * @brief code of the secondary structure type
         */
        char code;

        /**
         * @brief main group of the secondary structure type
         */
        char group;

        /**
         * @brief check if two secondary structure elements are equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator == (const SecStruct &other) {
            return (code == other.code);
        }

        /**
         * @brief check if two secondary structure elements are not equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator != (const SecStruct &other) {
            return (code != other.code);
        }

        /**
         * @brief check if type is a strict alpha type
         * 
         * @return true 
         * @return false 
         */
        bool isStrictAlpha() const {
            return (code == 'H' || code == 'G' || code == 'I');
        }

        /**
         * @brief check if type is alpha type
         * 
         * @return true 
         * @return false 
         */
        bool isAlpha() const {
            return (code == 'H' || code == 'G' || code == 'I' || code == 'T' || code == 'S');
        }

        /**
         * @brief check if type is beta type
         * 
         * @return true 
         * @return false 
         */
        bool isBeta() const {
            return (code == 'E' || code == 'S');
        }

        /**
         * @brief check if type is turn
         * 
         * @return true 
         * @return false 
         */
        bool isTurn() const {
            return (code == 'T' || code == 'B' );
        }

        /**
         * @brief check if type is strict turn
         * 
         * @return true 
         * @return false 
         */
        bool isStrictTurn() const {
            return (code == 'T' || code == 'B' || code == 'S');
        }

        /**
         * @brief check if two secondary structure types are the same type
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool same(const SecStruct &other) const {
            return ((isAlpha() && other.isAlpha()) || (isBeta() && other.isBeta()));
        }

        bool isUnStructured() const {
            return (code == '-');
        }
    };

    namespace SecStructType {
        const SecStruct H = { "Helix", 'H', 'H' };
        const SecStruct G = { "GHelix", 'G', 'H' };
        const SecStruct I = { "IHelix", 'I', 'H' };
        const SecStruct T = { "Turn", 'T', 'H' };
        const SecStruct B = { "Bend", 'B', 'B' };
        const SecStruct E = { "Extended", 'E', 'B' };
        const SecStruct S = { "SBend", 'S', 'B' };
        const SecStruct U = { "UnStructured", '-', '-' };
    }

    const std::unordered_map<char, const SecStruct> SecStructs = {
        { 'H', SecStructType::H },
        { 'G', SecStructType::G },
        { 'I', SecStructType::I },
        { 'T', SecStructType::T },
        { 'B', SecStructType::B },
        { 'E', SecStructType::E },
        { 'S', SecStructType::S },
        { '-', SecStructType::U },
    };
}
