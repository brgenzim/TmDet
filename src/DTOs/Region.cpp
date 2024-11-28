#include <string>
#include <format>
#include <DTOs/Region.hpp>
#include <ValueObjects/Region.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

        std::string Region::toString(const Tmdet::ValueObjects::Region& region) {
            return std::format(R"(
    REGION begin: {}-{}-{} end: {}-{}-{} type: {}
)",
                region.beg.authId,region.beg.authIcode,region.beg.labelId,
                region.end.authId,region.end.authIcode,region.end.labelId,
                region.type.code);
        }
}
