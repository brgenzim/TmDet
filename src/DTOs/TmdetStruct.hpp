#ifndef __TMDET_DTOS_TMDET__
#define __TMDET_DTOS_TMDET__

#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

namespace Tmdet::DTOS {

    struct TmdetStruct {
        static void readXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, string path);
        static void writeXml(Tmdet::ValueObjects::TmdetStruct& tmdetVO, string path);
        static void parse(Tmdet::ValueObjects::TmdetStruct& tmdetVO);
        static void out(Tmdet::ValueObjects::TmdetStruct& tmdetVO);

        static bool compareResidues(const gemmi::Residue& res1, const gemmi::Residue& res2);

        private:

        static std::vector<string> getChainSequence(const Tmdet::ValueObjects::TmdetStruct& tmdetVO,
            // const Tmdet::ValueObjects::Chain& chainVO);
            const gemmi::Chain& chainVO);

        static int fillResidueGapIfNeeded(const Tmdet::ValueObjects::TmdetStruct& tmdetVO,
            Tmdet::ValueObjects::Chain& chainVO, gemmi::Residue& residue,
            int from, const std::vector<string> sequence);

        static std::vector<gemmi::Residue> simpleResidueGapFill(gemmi::Chain& chain, gemmi::Residue& residue,
            int residueIndex, const std::vector<string> sequence);

    };
}

#endif