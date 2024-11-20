#pragma once

#include <string>
#include <ValueObjects/Chain.hpp>


namespace Tmdet::DTOs {

    struct Dssp {

        static std::string getSecondaryStructure(const Tmdet::ValueObjects::Chain& chain);
    };
}
