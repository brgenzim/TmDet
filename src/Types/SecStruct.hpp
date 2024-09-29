#ifndef __TMDET_TYPES_SECSTRUCT__
#define __TMDET_TYPES_SECSTRUCT__

#include <unordered_map>
#include <string>

namespace Tmdet::Types {

    struct SecStruct {
        std::string name;
        char code;
        char group;

        bool operator == (const SecStruct &other) {
            return (code == other.code);
        }

        bool operator != (const SecStruct &other) {
            return (code != other.code);
        }

        bool isAlpha() {
            return (code == 'H' || code == 'G' || code == 'I');
        }

        bool isBeta() {
            return (code == 'E');
        }

        bool isTurn() {
            return (code == 'T' || code == 'B' || code == 'S');
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
        const SecStruct U = { "Unknown", '-', '-' };
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

#endif
