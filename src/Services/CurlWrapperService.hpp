// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <string>

/**
 * @brief namespace for tmdet services
 *
 * @namespace Tmdet
 * @namespace Services
 */
namespace Tmdet::Services::CurlWrapperService {

    enum class Status {
        Ok,
        Error
    };

    std::string apiCall(std::string url, Status& resultCode);
    Status download(std::string url, std::string destination);
}

