#ifndef __TMDET_SERVICES_UNITMP__
#define __TMDET_SERVICES_UNIMTP__

#include <string>
#include <vector>
#include <Services/CurlWrapperService.hpp>
#include <ValueObjects/Region.hpp>

using Tmdet::ValueObjects::Region;

namespace Tmdet::Services::UniTmpService {

    std::vector<Region> getCctopPredictionResult(std::string code, CurlWrapperService::Status& resultCode);

}

#endif
