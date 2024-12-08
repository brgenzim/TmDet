#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>

namespace Tmdet::Utils {

    class Dssp {
        private:
            Tmdet::VOs::Protein& protein;
            void calcDsspOnChain(Tmdet::VOs::Chain& chain);
            void writeDsspOnChain(Tmdet::VOs::Chain& chain);
            void createHydrogenBonds(Tmdet::VOs::Chain& chain);
            void scanNeighbors(Tmdet::VOs::Chain& chain, int r1, const gemmi::Atom* CA, const gemmi::Atom* N, gemmi::Position hcoord);
            void setHydrogenBond(Tmdet::VOs::Residue& donor, Tmdet::VOs::Residue& akceptor, double energy);
            void detectTurns(Tmdet::VOs::Chain& chain, int d, std::string key);
            bool checkHbond1(Tmdet::VOs::Chain& chain, Tmdet::VOs::Residue& res, int d);
            bool checkHbond2(Tmdet::VOs::Chain& chain, Tmdet::VOs::Residue& res, int d);
            void initPbs(Tmdet::VOs::Chain& chain);
            void detectSecStructH(Tmdet::VOs::Chain& chain,std::string key);
            void detectSecStructG(Tmdet::VOs::Chain& chain,std::string key);
            void detectSecStructI(Tmdet::VOs::Chain& chain,std::string key);
            void detectSecStructT(Tmdet::VOs::Chain& chain);
            void detectSecStructS(Tmdet::VOs::Chain& chain);
            void detectSecStructBE(Tmdet::VOs::Chain& chain);
            bool checkPb(Tmdet::VOs::Chain& chain, int i, int j);
            bool checkApb(Tmdet::VOs::Chain& chain, int i, int j);
            bool checkIfAreOther(Tmdet::VOs::Chain& chain, Tmdet::Types::SecStruct ss,int i, int d);
            bool checkIfTurn(Tmdet::VOs::Chain& chain,int pos, int r, std::string key);
            void exec();
            void init();
            void end();

        public:
            explicit Dssp(Tmdet::VOs::Protein& protein) : 
                protein(protein) {
                    exec();
            } ;
            ~Dssp()=default;

    };
}
