#pragma once

#include <string>
#include <gemmi/model.hpp>
#include <ValueObjects/Protein.hpp>

namespace StructVO = Tmdet::ValueObjects;

// #define __DSSP_DEBUG 1

namespace Tmdet::Utils {

    class Dssp {
        private:
            StructVO::Protein& _proteinVO;
            void calcDsspOnChain(StructVO::Chain& chain);
            void writeDsspOnChain(StructVO::Chain& chain);
            void createHydrogenBonds(StructVO::Chain& chain);
            void scanNeighbors(StructVO::Chain& chain, int r1, const gemmi::Atom* CA, const gemmi::Atom* N, gemmi::Position hcoord);
            void setHydrogenBond(StructVO::Residue& donor, StructVO::Residue& akceptor, double energy);
            void detectTurns(StructVO::Chain& chain, int d, std::string key);
            bool checkHbond1(StructVO::Residue& res, int d);
            bool checkHbond2(StructVO::Residue& res, int d);
            void initPbs(StructVO::Chain& chain);
            void detectSecStructH(StructVO::Chain& chain,std::string key);
            void detectSecStructG(StructVO::Chain& chain,std::string key);
            void detectSecStructI(StructVO::Chain& chain,std::string key);
            void detectSecStructT(StructVO::Chain& chain);
            void detectSecStructS(StructVO::Chain& chain);
            void detectSecStructBE(StructVO::Chain& chain);
            bool checkPb(StructVO::Chain& chain, int i, int j);
            bool checkApb(StructVO::Chain& chain, int i, int j);
            bool checkIfAreOther(StructVO::Chain& chain, Tmdet::Types::SecStruct ss,int i, int d);
            bool checkIfTurn(StructVO::Chain& chain,int pos, int r, std::string key);

        public:
            explicit Dssp(StructVO::Protein& proteinVO) : 
                _proteinVO(proteinVO) {
                    exec();
            } ;
            ~Dssp()=default;

            void exec();
            static std::string getSecStructAsString(StructVO::Chain& chain);
    };
}
