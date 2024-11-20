#include <string>
#include <format>
#include <ValueObjects/Region.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

        std::string toString(Tmdet::ValueObjects::Region& region) {
            return std::format(R"(
    REGION begin: {} end: {} type: {}
)",
                region.beg,region.end,region.type.code);
        }
}
