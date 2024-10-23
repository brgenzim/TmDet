#pragma once

#include <string>
#include <vector>
#include <format>
#include <Utils/SecStrVec.hpp>

namespace Tmdet::Helpers::Pymol {

    struct _secStrVec;

    static std::string dumpSecStrVec(std::vector<Tmdet::Utils::_secStrVec>& vectors, std::string color) {
        int counter = 1;
        std::string ret;
        for (auto& vector : vectors) {
            ret += std::format("cgo_arrow [ {}, {}, {}], [ {}, {}, {}], name={}_{}, color={}",
                    vector.begin.x, vector.begin.y, vector.begin.z,
                    vector.end.x, vector.end.y, vector.end.z,
                    vector.type.name, counter, color);
            counter++;
        }
        return ret;
    }
}