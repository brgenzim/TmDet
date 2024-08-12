#ifndef __TMDET_SERVICES_UNITMP__
#define __TMDET_SERVICES_UNIMTP__

#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <Services/CurlWrapperService.hpp>
#include <ValueObjects/Region.hpp>

using Tmdet::ValueObjects::Region;

namespace Tmdet::Services::UniTmpService {

    // const std::string UNITMP_SCHEMA = "https://";
    // const std::string UNITMP_DOMAIN = "unitmp.org";
    const std::string UNITMP_SCHEMA = "http://";
    const std::string UNITMP_DOMAIN = "unitmp.localhost";
    // http://cctop.unitmp.localhost/api/v1/prediction/Cctop/CCG2_MOUSE.json

    std::vector<Region> getCctopPredictionResult(std::string code, CurlWrapperService::Status& resultCode) {
        // TODO: get values from configuration service
        std::string url(UNITMP_SCHEMA);
        url += "cctop." + UNITMP_DOMAIN + "/api/v1/prediction/Cctop/" + code + ".json";

        std::string response = CurlWrapperService::apiCall(url, resultCode);
        // TODO: check result code

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
