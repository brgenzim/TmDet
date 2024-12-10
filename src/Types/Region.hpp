#pragma once

#include <unordered_map>
#include <string>

/**
 * @brief namespace of types
 * @namespace Tmdet
 * @namespace Types
 */
namespace Tmdet::Types {

    struct Region {
        int id;
        std::string name;
        char code;
        std::string description;

        bool operator == (const Region& other) const {
            return id == other.id;
        }

        bool isAnnotatedTransMembraneType() const {
            return code == 'H' || code == 'B';
        }

        bool isAnnotatedMembraneType() const {
            return code == 'H' || code == 'B' || code == 'L' || code =='F';
        }

        bool isNotMembrane() const {
            return code == '1' || code == '2' || code == '3';
        }

        bool isNotAnnotatedMembrane() const {
            return code == 'M';
        }

        bool isMembraneInside() const {
            return code == 'N';
        }

        bool isAlpha() const {
            return code == 'H';
        }
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
        const Region IFH = {
            6,
            "interfacial helix",
            'F',
            "Interfacial helix on membrane surface"
        };
        const Region MEMBINS = {
            7,
            "membins",
            'N',
            "Beta barrel inside element"
        };
        const Region INTERMEMB = {
            8,
            "intermembrane",
            '3',
            "Intermembrane space if protein crosses two membranes"
        };
        const Region UNK = {
            9,
            "unkown",
            'U',
            "Unknown localization"
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
        const Region ERROR = {
            12,
            "error",
            'R',
            "errorneously placed membrane region"
        };
        const Region PERIPLASM = {
            13,
            "periplasm",
            'E',
            "Periplasmic"
        };
    }

    const std::unordered_map<char, const Region> Regions = {
        { 'M', RegionType::MEMB },
        { 'N', RegionType::MEMBINS },
        { 'H', RegionType::HELIX },
        { 'B', RegionType::BETA },
        { '1', RegionType::SIDE1 },
        { '2', RegionType::SIDE2 },
        { '3', RegionType::INTERMEMB },
        { 'L', RegionType::LOOP },
        { 'F', RegionType::IFH },
        { 'I', RegionType::INSIDE },
        { 'O', RegionType::OUTSIDE },
        { 'E', RegionType::PERIPLASM },
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
