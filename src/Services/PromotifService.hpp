#ifndef __TMDET_SERVICES_PROMOTIF__
#define __TMDET_SERVICES_PROMOTIF__

#include <map>
#include <string>

namespace Tmdet::Services::PromotifService {

    std::map<char, std::string> process(const std::string& cifPath);

}

#endif
