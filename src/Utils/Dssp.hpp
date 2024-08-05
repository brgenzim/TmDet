#ifndef __TMDET_UTILS_DSSP__
#define __TMDET_UTILS_DSSP__

#include <array>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>

using namespace std;

// #define __DSSP_DEBUG 1

namespace Tmdet::Utils {

    class Dssp {
        private:
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;
            void calcDsspOnChain(Tmdet::ValueObjects::Chain& chain);
            void writeDsspOnChain(Tmdet::ValueObjects::Chain& chain);
            void createHydrogenBonds(Tmdet::ValueObjects::Chain& chain);
            void scanNeighbors(Tmdet::ValueObjects::Chain& chain, int r1, const gemmi::Atom* CA, const gemmi::Atom* N, gemmi::Position hcoord);
            void setHydrogenBond(Tmdet::ValueObjects::Residue& donor, Tmdet::ValueObjects::Residue& akceptor, double energy);
            void detectTurns(Tmdet::ValueObjects::Chain& chain, int d, string key);
            bool checkHbond1(Tmdet::ValueObjects::Residue& res, int d);
            bool checkHbond2(Tmdet::ValueObjects::Residue& res, int d);
            void initPbs(Tmdet::ValueObjects::Chain& chain);
            void detectSecStructH(Tmdet::ValueObjects::Chain& chain,string key);
            void detectSecStructG(Tmdet::ValueObjects::Chain& chain,string key);
            void detectSecStructI(Tmdet::ValueObjects::Chain& chain,string key);
            void detectSecStructT(Tmdet::ValueObjects::Chain& chain);
            void detectSecStructS(Tmdet::ValueObjects::Chain& chain);
            void detectSecStructBE(Tmdet::ValueObjects::Chain& chain);
            bool checkPb(Tmdet::ValueObjects::Chain& chain, int i, int j);
            bool checkApb(Tmdet::ValueObjects::Chain& chain, int i, int j);
            bool checkIfAreOther(Tmdet::ValueObjects::Chain& chain, Tmdet::Types::SecStruct ss,int i, int d);
            bool checkIfTurn(Tmdet::ValueObjects::Chain& chain,int pos, int r, string key);

        public:
            Dssp(Tmdet::ValueObjects::TmdetStruct& _tmdetVO) : tmdetVO(_tmdetVO) {} ;
            ~Dssp() {};

            void calcDsspOnStructure();
            void writeDsspOnStructure();
            static string getDsspOfChain(Tmdet::ValueObjects::Chain& chain);
    };
}
#endif
