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
     * @brief definition of a region type
     */
    struct Region {
        /**
         * @brief uniq udentifier number
         */
        int id;

        /**
         * @brief name of the region type
         */
        std::string name;

        /**
         * @brief code of the region type
         */
        char code;

        /**
         * @brief description of the region type
         */
        std::string description;

        /**
         * @brief check if two region type are equal
         * 
         * @param other 
         * @return true 
         * @return false 
         */
        bool operator == (const Region& other) const {
            return id == other.id;
        }

        /**
         * @brief check if the region type is annotated transmembrane type
         * 
         * @return true 
         * @return false 
         */
        bool isAnnotatedTransMembraneType() const {
            return code == 'H' || code == 'B';
        }

        /**
         * @brief check if the region type is annotated membrane type
         * 
         * @return true 
         * @return false 
         */
        bool isAnnotatedMembraneType() const {
            return code == 'H' || code == 'B' || code == 'L' || code =='F';
        }

        /**
         * @brief check if region type is not membrane
         * 
         * @return true 
         * @return false 
         */
        bool isNotMembrane() const {
            return code == '1' || code == '2' || code == '3';
        }

        /**
         * @brief check if region type is not annotated membrane type
         * 
         * @return true 
         * @return false 
         */
        bool isNotAnnotatedMembrane() const {
            return code == 'M';
        }

        /**
         * @brief check if region type is membrane inside
         * 
         * @return true 
         * @return false 
         */
        bool isMembraneInside() const {
            return code == 'N';
        }

        /**
         * @brief check if region type is alpha helical transmembrane type
         * 
         * @return true 
         * @return false 
         */
        bool isAlpha() const {
            return code == 'H';
        }

        /**
         * @brief check if region type is beta barrel type
         * 
         * @return true 
         * @return false 
         */
        bool isBeta() const {
            return code == 'B';
        }
    };

    namespace RegionType {
        const Region MEMB = {
            0,
            "membrane", 
            'M',
            "Region within the membrane"
        };
        const Region HELIX = {
            1,
            "Alpha helical transmembrane", 
            'H',
            "Alpha helical region crossing the membrane"
        };
        const Region BETA = {
            2,
            "Beta barrel strain", 
            'B',
            "Beta strain crossing the membrane"
        };
        const Region SIDE1 = {
            3,
            "side1",
            '1',
            "Side one of the membrane"
        };
        const Region SIDE2 = {
            4,
            "side2",
            '2',
            "Other side of the membrane"
        };
        const Region LOOP = {
            5,
            "loop",
            'L',
            "Re-entrant membrane loop"
        };
        const Region TWO_H_LOOP = {
            6,
            "loop with two helices",
            'L',
            "Re-entrant membrane loop with two helices"
        };
        const Region IFH = {
            7,
            "interfacial helix",
            'F',
            "Interfacial helix on membrane surface"
        };
        const Region MEMBINS = {
            8,
            "membins",
            'N',
            "Beta barrel inside element"
        };
        const Region INTERMEMB = {
            9,
            "intermembrane",
            '3',
            "Intermembrane space if protein crosses two membranes"
        };
        const Region INSIDE = {
            10,
            "inside",
            'I',
            "inside (i.e cytoplasmic) in TOPDB"
        };
        const Region OUTSIDE = {
            11,
            "outside",
            'O',
            "outside (i.e extra-cytosolic) in TOPDB"
        };
        const Region PERIPLASM = {
            12,
            "periplasm",
            'E',
            "Periplasmic"
        };
        const Region ERROR_FP = {
            13,
            "false positiv error",
            'P',
            "region in membrane that should not be in the membrane"
        };
        const Region ERROR_FN = {
            14,
            "false negativ error",
            'R',
            "region outside of the membrane that should be in the membrane"
        };
        const Region UNK = {
            15,
            "unknown",
            'U',
            "localisation of the region is unknown"
        };
        
    }

    const std::unordered_map<char, const Region> Regions = {
        { 'M', RegionType::MEMB },
        { 'H', RegionType::HELIX },
        { 'B', RegionType::BETA },
        { '1', RegionType::SIDE1 },
        { '2', RegionType::SIDE2 },
        { 'L', RegionType::LOOP },
        { 'K', RegionType::TWO_H_LOOP},
        { 'F', RegionType::IFH },
        { 'N', RegionType::MEMBINS },
        { '3', RegionType::INTERMEMB },
        { 'I', RegionType::INSIDE },
        { 'O', RegionType::OUTSIDE },
        { 'E', RegionType::PERIPLASM },
        { 'P', RegionType::ERROR_FP },
        { 'R', RegionType::ERROR_FN },
        { 'U', RegionType::UNK }
    };

    /**
     * Regions returned by CCTOP API call.
     */
    const std::unordered_map<std::string, const Region> RegionsByName = {
        { "Membrane", RegionType::MEMB },
        { "Membins", RegionType::MEMBINS },
        { "Re-entrant loop", RegionType::LOOP },
        { "Interfacial helix", RegionType::IFH },
        { "Inside", RegionType::INSIDE },
        { "Outside", RegionType::OUTSIDE },
        { "Periplasm", RegionType::PERIPLASM },
    };
}
