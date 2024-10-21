#ifndef __TMDET_TYPES_REGION__
#define __TMDET_TYPES_REGION__

#include <unordered_map>
#include <string>

namespace Tmdet::Types {

    struct Region {
        std::string name;
        char code;
        std::string description;

        bool operator==(Region other) {
            return code == other.code;
        }
    };

    namespace RegionType {
        const Region MEMB = {
            "transmembrane", 
            'M',
            "Region crossing the membrane"
        };
        const Region HELIX = {
            "Alpha helical transmembrane", 
            'H',
            "Alpha helical region crossing the membrane"
        };
        const Region BETA = {
            "Beta barrel strain", 
            'B',
            "Beta strain crossing the membrane"
        };
        const Region SIDE1 = {
            "side1",
            '1',
            "Side one of the membrane"
        };
        const Region SIDE2 = {
            "side2",
            '2',
            "Other side of the membrane"
        };
        const Region LOOP = {
            "loop",
            'L',
            "Re-entrant membrane loop"
        };
        const Region IFH = {
            "interfacial helix",
            'F',
            "Interfacial helix on membrane surface"
        };
        const Region MEMBINS = {
            "membins",
            'N',
            "Beta barrel inside element"
        };
        const Region INTERMEMB = {
            "intermembrane",
            '3',
            "Intermembrane space if protein crosses two membranes"
        };
        const Region UNK = {
            "unkown",
            'U',
            "Unknown localization"
        };
        const Region INSIDE = {
            "inside",
            'I',
            "inside (i.e cytoplasmic) in TOPDB"
        };
        const Region OUTSIDE = {
            "outside",
            'O',
            "outside (i.e extra-cytosolic) in TOPDB"
        };
        const Region PERIPLASM = {
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

#endif
