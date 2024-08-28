#ifndef __TMDET_TYPES_REGION__
#define __TMDET_TYPES_REGION__

#include <unordered_map>
#include <string>

namespace Tmdet::Types {

    struct Region {
        std::string name;
        char code;
        std::string description;
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
        const Region INSIDE = {
            "inside",
            'I',
            "Beta barrel inside element"
        };
        const Region BEG = {
            "chain begin",
            'N',
            "N-terminal chain segment"
        };
        const Region END = {
            "chain end",
            'C',
            "C-terminal chain segment"
        };
        const Region UNK = {
            "unkown",
            'U',
            "Unknow localization"
        };
    }

    const std::unordered_map<char, const Region> Regions = {
        { 'M', RegionType::MEMB },
        { 'H', RegionType::HELIX },
        { 'B', RegionType::BETA },
        { '1', RegionType::SIDE1 },
        { '2', RegionType::SIDE2 },
        { 'L', RegionType::LOOP },
        { 'F', RegionType::IFH },
        { 'I', RegionType::INSIDE },
        { 'N', RegionType::BEG },
        { 'C', RegionType::END },
        { 'U', RegionType::UNK }
    };
}

#endif
