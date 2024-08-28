#ifndef __TMDET_DTOS_TMDET__
#define __TMDET_DTOS_TMDET__

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

namespace Tmdet::DTOS {

    struct TmdetStruct {
        static void readXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, std::string path);
        static void writeXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, std::string path);
        static void writeCif(Tmdet::ValueObjects::TmdetStruct& tmdetVO, std::string path);
        static void parse(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
        static void out(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
        static std::vector<std::string> getChainSequence(const Tmdet::ValueObjects::TmdetStruct& tmdetVO,
            const gemmi::Chain& chainVO);
    };
}

#endif
