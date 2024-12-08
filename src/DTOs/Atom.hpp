#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Atom.hpp>

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
        static std::string toString(const Tmdet::VOs::Atom& atom);
    };
}
