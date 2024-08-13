#ifndef __TMDET_SERVICES_CURL_WRAPPER__
#define __TMDET_SERVICES_CURL_WRAPPER__

#include <string>

namespace Tmdet::Services::CurlWrapperService {

    enum class Status {
        Ok,
        Error
    };

    std::string apiCall(std::string url, Status& resultCode);
    Status download(std::string url, std::string destination);
}

#endif
