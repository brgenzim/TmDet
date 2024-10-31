#pragma once

#include <string>
#include <format>
#include <ValueObjects/SecStrVec.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    struct SecStrVec {

        /**
         * @brief print value object content for debugging purpose
         */
        static std::string print(Tmdet::ValueObjects::SecStrVec& vec) {
            return std::format(R"(Tmdet::ValueObjects::SecStrVec(
    type: {}
    begin: [{}, {}, {}]
    end: [{}, {}, {}]
    chainIdx: {}
    begResIdx: {}
    endResIdx: {})",
                vec.type.name,vec.begin.x,vec.begin.y,vec.begin.z,
                vec.end.x, vec.end.y, vec.end.z,
                vec.chainIdx,vec.begResIdx,vec.endResIdx);
        }
    };
}