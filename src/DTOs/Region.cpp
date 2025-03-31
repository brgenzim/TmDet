// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <string>
#include <format>
#include <DTOs/Region.hpp>
#include <VOs/Region.hpp>

/**
 * @namespace Tmdet
 * @namespace DTOs
 */
namespace Tmdet::DTOs {

        std::string Region::toString(const Tmdet::VOs::Region& region) {
            return std::format(R"(
    REGION begin: {}-{}-{} end: {}-{}-{} type: {}
)",
                region.beg.authId,region.beg.authIcode,region.beg.labelId,
                region.end.authId,region.end.authIcode,region.end.labelId,
                region.type.code);
        }
}
