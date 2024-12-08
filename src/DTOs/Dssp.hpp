#pragma once

#include <string>
#include <VOs/Chain.hpp>


namespace Tmdet::DTOs {

    struct Dssp {

        static std::string getSecondaryStructure(const Tmdet::VOs::Chain& chain);
    };
}
