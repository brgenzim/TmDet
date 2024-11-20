#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <ValueObjects/Atom.hpp>

/**
 * @brief namespace for data transfer objects
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    /**
     * @brief writing atom data
     */
    struct Atom {

        /**
         * @brief string representation of the atom
         * 
         * @param atom
         */
        static std::string toString(const Tmdet::ValueObjects::Atom& atom);
    };
}
