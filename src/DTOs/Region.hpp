#pragma once

#include <string>
#include <format>
#include <ValueObjects/Region.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

    struct Region {

        /**
         * @brief print value object content for debugging purpose
         */
        static std::string print(Tmdet::ValueObjects::Region& region) {
            return std::format(R"(Tmdet::ValueObjects::Region(
    beg: {}
    beg_auth_seq_id: {}
    beg_auth_seq_icode: {}
    beg_label_seq_id: {}
    end: {}
    end_auth_seq_id: {}
    end_auth_seq_icode: {}
    end_label_seq_id: {}
    type: {})",
            region.beg,region.beg_auth_seq_id,
            region.beg_auth_seq_icode,region.beg_label_seq_id,
            region.end,region.end_auth_seq_id,
            region.end_auth_seq_icode,region.end_label_seq_id,
            region.type.code);
        }
    };
}