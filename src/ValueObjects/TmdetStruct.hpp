#ifndef __TMDET_VALUE_OBJECTS_TMDETSTRUCT__
#define __TMDET_VALUE_OBJECTS_TMDETSTRUCT__

#include <string>
#include <vector>
#include <ValueObjects/Modification.hpp>
#include <ValueObjects/BioMatrix.hpp>
#include <ValueObjects/Membrane.hpp>
#include <ValueObjects/Chain.hpp>
#include <Types/Protein.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>

using namespace std;

namespace Tmdet::ValueObjects {

    struct TmdetStruct {
        string code;
        string inputPath;
        bool tmp;
        string date;
        vector<Modification> modifications;
        double qValue;
        Tmdet::Types::Protein type;
        string spres;
        string pdbkwres;
        BioMatrix bioMatrix;
        vector<Membrane> membranes;
        vector<Chain> chains;
        gemmi::Structure& gemmi;
        gemmi::cif::Document& document;
        gemmi::NeighborSearch neighbors;

        TmdetStruct(TmdetStruct&& other) noexcept :
            code(other.code),
            inputPath(other.inputPath),
            tmp(other.tmp),
            date(other.date),
            modifications(other.modifications),
            qValue(other.qValue),
            type(other.type),
            spres(other.spres),
            pdbkwres(other.pdbkwres),
            bioMatrix(other.bioMatrix),
            membranes(other.membranes),
            chains(other.chains),
            gemmi(other.gemmi),
            document(other.document),
            neighbors(other.neighbors)
        {
            auto& otherChains = other.gemmi.models[0].chains;
            for (int chain = 0;  chain < otherChains.size(); ++chain) {
                auto& residues = otherChains[chain].residues;
                for (int residue = 0; residue < residues.size(); ++residue) {
                    gemmi.models[0].chains[chain].residues[residue].atoms = vector<gemmi::Atom>(residues[residue].atoms);
                }
            }
        }

        TmdetStruct(gemmi::Structure& _gemmi, gemmi::cif::Document& _document) : gemmi(_gemmi), document(_document) {
        }
        ~TmdetStruct() {}
    };
}

#endif