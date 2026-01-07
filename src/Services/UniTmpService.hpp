// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>
#include <vector>
#include <Services/CurlWrapperService.hpp>
#include <ValueObjects/Region.hpp>

using Tmdet::ValueObjects::Region;

/**
 * @brief namespace for tmdet services
 *
 * @namespace Tmdet
 * @namespace Services
 */
namespace Tmdet::Services::UniTmpService {

    std::vector<Region> getCctopPredictionResult(std::string code, CurlWrapperService::Status& resultCode);

}
