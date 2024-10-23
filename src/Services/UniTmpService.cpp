#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <Config.hpp>
#include <Services/CurlWrapperService.hpp>
#include <ValueObjects/Region.hpp>

namespace TmdetVO = Tmdet::ValueObjects;

namespace Tmdet::Services::UniTmpService {

    std::vector<TmdetVO::Region> getCctopPredictionResult(std::string code, CurlWrapperService::Status& resultCode) {
        std::string url = environment.get("UNITMP_SCHEMA",DEFAULT_UNITMP_SCHEMA)
                    + "cctop." 
                    + environment.get("UNITMP_DOMAIN",DEFAULT_UNITMP_DOMAIN)
                    + "/api/v1/prediction/Cctop/" 
                    + code 
                    + ".json";

        std::string response = CurlWrapperService::apiCall(url, resultCode);
        if (resultCode == Tmdet::Services::CurlWrapperService::Status::Error) {
            throw std::runtime_error("HTTP request to '" + url + "' failed");
        }

        auto json = nlohmann::json::parse(response);

        auto result = std::vector<TmdetVO::Region>();
        for (auto& item : json["Regions"]) {
            TmdetVO::Region region;
            region.beg = item["start"].template get<int>();
            region.end = item["end"].template get<int>();
            region.type = Tmdet::Types::RegionsByName.at(item["loc"]);
            result.emplace_back(region);
        }

        return result;
    }
}
