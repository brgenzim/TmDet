#ifndef __TMDET_SERVICES_UNITMP__
#define __TMDET_SERVICES_UNIMTP__

#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <Services/ConfigurationService.hpp>
#include <Services/CurlWrapperService.hpp>
#include <ValueObjects/Region.hpp>

using Tmdet::ValueObjects::Region;

namespace Tmdet::Services::UniTmpService {

    std::vector<Region> getCctopPredictionResult(std::string code, CurlWrapperService::Status& resultCode) {
        std::string url(ConfigurationService::UniTmpSchema);
        url += "cctop." + ConfigurationService::UniTmpDomain + "/api/v1/prediction/Cctop/" + code + ".json";

        std::string response = CurlWrapperService::apiCall(url, resultCode);
        if (resultCode == Tmdet::Services::CurlWrapperService::Status::Error) {
            throw runtime_error("HTTP request to '" + url + "' failed");
        }

        auto json = nlohmann::json::parse(response);

        auto result = std::vector<Region>();
        for (auto& item : json["Regions"]) {
            Region region;
            region.rbeg = item["start"].template get<int>();
            region.rend = item["end"].template get<int>();
            region.beg = region.end = 0;
            region.begi = region.endi = ' ';
            region.type = Tmdet::Types::RegionsByName.at(item["loc"]);
            result.emplace_back(region);
        }

        return result;
    }
}

#endif
