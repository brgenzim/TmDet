// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <iostream>
#include <sstream>
#include <array>
#include <math.h>
#include <gemmi/model.hpp>
#include <gemmi/neighbor.hpp>
#include <Config.hpp>
#include <DTOs/Residue.hpp>
#include <DTOs/Dssp.hpp>
#include <System/Logger.hpp>
#include <VOs/Protein.hpp>
#include <VOs/HBond.hpp>
#include <Utils/Dssp.hpp>

using namespace std;
using namespace gemmi;

namespace Tmdet::Utils {

#define DSSP_Q -27888
#define DSSP_PDB_DL 0.5
#define DSSP_HBLOW -9900
#define DSSP_HBHIGH -500

    void Dssp::exec() {
        init();
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                calcDsspOnChain(chain);
            }
        );
        end();
    }

    void Dssp::init() {
        protein.eachResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.try_emplace("hbond1",std::any(Tmdet::VOs::HBond()));
                residue.temp.try_emplace("hbond2",std::any(Tmdet::VOs::HBond()));
                residue.temp.try_emplace("t3",std::any(' '));
                residue.temp.try_emplace("t4",std::any(' '));
                residue.temp.try_emplace("t5",std::any(' '));
                residue.temp.try_emplace("pb",std::any(-1));
                residue.temp.try_emplace("apb",std::any(-1));
            }
        );
    }

    void Dssp::end() {
        protein.eachResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("hbond2");
                residue.temp.erase("t3");
                residue.temp.erase("t4");
                residue.temp.erase("t5");
                residue.temp.erase("pb");
                residue.temp.erase("apb");
            }
        );
    }

    void Dssp::calcDsspOnChain(Tmdet::VOs::Chain& chain) {
        createHydrogenBonds(chain);
        detectTurns(chain,3,"t3");
        detectTurns(chain,4,"t4");
        detectTurns(chain,5,"t5");
        initPbs(chain);
        detectSecStructH(chain,"t4");
        detectSecStructG(chain,"t3");
        detectSecStructI(chain,"t5");
        detectSecStructT(chain);
        detectSecStructS(chain);
        detectSecStructBE(chain);
    }

    void Dssp::createHydrogenBonds(Tmdet::VOs::Chain& chain) {
        for(int i=1; i<chain.length; i++) {
            if (chain.orderDistance(i-1,i) == 1) {
                auto& gres = chain.residues[i].gemmi;
                auto const& prevGRes = chain.residues[i-1].gemmi;
                auto CA = gres.get_ca();
                auto CE = prevGRes.get_c();
                auto OE = prevGRes.find_atom("O",' ');
                auto N = gres.get_n();
                if (CA != (Atom *)nullptr && CE != (Atom *)nullptr && OE != (Atom *)nullptr && N != (Atom *)nullptr) {
                    double dco = CE->pos.dist(OE->pos);
                    auto hcoord = Position(
                        N->pos.x+(CE->pos.x-OE->pos.x)/dco,
                        N->pos.y+(CE->pos.y-OE->pos.y)/dco,
                        N->pos.z+(CE->pos.z-OE->pos.z)/dco
                    );
                    scanNeighbors(chain, i, CA, N, hcoord);
                }
            }
        }
    }

    void Dssp::scanNeighbors(Tmdet::VOs::Chain& chain, int r1, const gemmi::Atom* CA, const gemmi::Atom* N, gemmi::Position hcoord) {
        for(int r2=0; r2<chain.length; r2++) {
            if (std::abs(chain.orderDistance(r1,r2)) > 1) {
                auto CB = chain.residues[r2].gemmi.get_ca();
                if (CB != (Atom *)nullptr && CA->pos.dist(CB->pos) < 9.0) {
                    auto C = chain.residues[r2].gemmi.get_c();
                    auto O = chain.residues[r2].gemmi.find_atom("O",' ');
                    if (C != (Atom *)nullptr && O != (Atom *)nullptr) {
                        double dho = hcoord.dist(O->pos);
                        double dhc = hcoord.dist(C->pos);
                        double dno = N->pos.dist(O->pos);
                        double dnc = N->pos.dist(C->pos);
                        double energy;
                        if (dho<DSSP_PDB_DL||dhc<DSSP_PDB_DL||dno<DSSP_PDB_DL||dnc<DSSP_PDB_DL) {
                            setHydrogenBond(chain.residues[r2],chain.residues[r1],DSSP_HBLOW);
                        }
                        else if ((energy=DSSP_Q/dho-DSSP_Q/dhc+DSSP_Q/dnc-DSSP_Q/dno+0.5)<DSSP_HBHIGH) {
                            setHydrogenBond(chain.residues[r2],chain.residues[r1],energy);
                        }
                    }
                }
            }
        }
    }

    void Dssp::setHydrogenBond(Tmdet::VOs::Residue& donor, Tmdet::VOs::Residue& akceptor, double energy) {
        if (energy < any_cast<Tmdet::VOs::HBond>(donor.temp["hbond1"]).energy) {
            donor.temp["hbond2"] = donor.temp["hbond1"];
            donor.temp["hbond1"] = std::any(Tmdet::VOs::HBond({energy, akceptor.chainIdx, akceptor.idx}));
        }
        else if (energy< any_cast<Tmdet::VOs::HBond>(donor.temp["hbond2"]).energy) {
            donor.temp["hbond2"] = std::any(Tmdet::VOs::HBond({energy, akceptor.chainIdx, akceptor.idx}));
        }
    }

    void Dssp::detectTurns(Tmdet::VOs::Chain& chain, int d, string key) {
        for(auto i=0; i<chain.length-d; i++) {
            if (chain.orderDistance(i+d,i) == d &&
                (checkHbond1(chain, chain.residues[i],d) || checkHbond2(chain, chain.residues[i],d))) {
                chain.residues[i].temp[key] = (any_cast<char>(chain.residues[i].temp[key]) == '<'?'x':'>');
                for(int j=1; j<d; j++) {
                    if (any_cast<char>(chain.residues[i+j].temp[key]) == ' ') {
                        chain.residues[i+j].temp[key] = '*';
                    }
                }
                chain.residues[i+d].temp[key] = '<';
            }
        }
        for(auto i=1; i<chain.length-1; i++) {
            if (chain.orderDistance(i-1,i) == 1 && chain.orderDistance(i,i+1) == 1) {
                auto& dtmpm = any_cast<char&>(chain.residues[i-1].temp[key]);
                auto& dtmp = any_cast<char&>(chain.residues[i].temp[key]);
                auto& dtmpp = any_cast<char&>(chain.residues[i+1].temp[key]);
                if ((dtmpm != '*' && dtmp == '*' && dtmpp != '*') &&
                    (dtmpm == 'x' || dtmpp == 'x')) {
                        if (dtmpm == '>') { dtmp = '>'; }
                        if (dtmpm == '<') { dtmp = '<'; }
                        if (dtmpp == '>') { dtmp = '>'; }
                        if (dtmpp == '<') { dtmp = '<'; }
                }
            }
        }
    }

    bool Dssp::checkHbond1(Tmdet::VOs::Chain& chain, Tmdet::VOs::Residue& res, int d) {
        return (res.chainIdx == any_cast<Tmdet::VOs::HBond>(res.temp["hbond1"]).toChainIdx && 
            chain.orderDistance(res.idx, any_cast<Tmdet::VOs::HBond>(res.temp["hbond1"]).toResIdx) == d);
    }

    bool Dssp::checkHbond2(Tmdet::VOs::Chain& chain, Tmdet::VOs::Residue& res, int d) {
        return (res.chainIdx == any_cast<Tmdet::VOs::HBond>(res.temp["hbond2"]).toChainIdx && 
            chain.orderDistance(res.idx, any_cast<Tmdet::VOs::HBond>(res.temp["hbond2"]).toResIdx) == d);
    }

    void Dssp::initPbs(Tmdet::VOs::Chain& chain) {
        int n = chain.length;
        for(int i=1; i<n-1; i++) {
            for(int j=1; j<n-1; j++) {
                if (abs(i-j) > 2) {
                    if (checkPb(chain,i,j)) {
                        chain.residues[i].temp["pb"] = any_cast<int>(chain.residues[j].idx);
                    }
                    if (checkApb(chain,i,j)) {
                        chain.residues[i].temp["apb"] = any_cast<int>(chain.residues[j].idx);
                    }
                }
            }
        }
    }

    bool Dssp::checkPb(Tmdet::VOs::Chain& chain, int i, int j) {
        return (((any_cast<Tmdet::VOs::HBond>(chain.residues[i-1].temp["hbond1"]).toResIdx == chain.residues[j].idx ||
                any_cast<Tmdet::VOs::HBond>(chain.residues[i-1].temp["hbond2"]).toResIdx == chain.residues[j].idx) &&
                (any_cast<Tmdet::VOs::HBond>(chain.residues[j].temp["hbond1"]).toResIdx == chain.residues[i+1].idx ||
                any_cast<Tmdet::VOs::HBond>(chain.residues[j].temp["hbond2"]).toResIdx == chain.residues[i+1].idx)) ||
                ((any_cast<Tmdet::VOs::HBond>(chain.residues[j-1].temp["hbond1"]).toResIdx == chain.residues[i].idx || 
                any_cast<Tmdet::VOs::HBond>(chain.residues[j-1].temp["hbond2"]).toResIdx == chain.residues[i].idx) &&
                (any_cast<Tmdet::VOs::HBond>(chain.residues[i].temp["hbond1"]).toResIdx == chain.residues[j+1].idx ||
                any_cast<Tmdet::VOs::HBond>(chain.residues[i].temp["hbond2"]).toResIdx == chain.residues[j+1].idx)));
    }

    bool Dssp::checkApb(Tmdet::VOs::Chain& chain, int i, int j) {
        return  (((any_cast<Tmdet::VOs::HBond>(chain.residues[i].temp["hbond1"]).toResIdx == chain.residues[j].idx ||
                    any_cast<Tmdet::VOs::HBond>(chain.residues[i].temp["hbond2"]).toResIdx == chain.residues[j].idx) &&
                    (any_cast<Tmdet::VOs::HBond>(chain.residues[j].temp["hbond1"]).toResIdx == chain.residues[i].idx ||
                    any_cast<Tmdet::VOs::HBond>(chain.residues[j].temp["hbond2"]).toResIdx == chain.residues[i].idx)) ||
                    ((any_cast<Tmdet::VOs::HBond>(chain.residues[i-1].temp["hbond1"]).toResIdx == chain.residues[j+1].idx ||
                    any_cast<Tmdet::VOs::HBond>(chain.residues[i-1].temp["hbond2"]).toResIdx == chain.residues[j+1].idx) &&
                    (any_cast<Tmdet::VOs::HBond>(chain.residues[i+1].temp["hbond1"]).toResIdx == chain.residues[j-1].idx ||
                    any_cast<Tmdet::VOs::HBond>(chain.residues[i+1].temp["hbond2"]).toResIdx == chain.residues[j-1].idx)));
    }

    void Dssp::detectSecStructH(Tmdet::VOs::Chain& chain,string key) {
        for( int i=1; i<chain.length-4; i++) {
            if (((any_cast<char>(chain.residues[i].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i].temp[key])=='x')) &&
                ((any_cast<char>(chain.residues[i-1].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i-1].temp[key])=='x'))) {
                    for(int j=0; j<4; j++) {
                        chain.residues[i+j].ss = Tmdet::Types::SecStructType::H;
                    }
                }
        }
    }

    void Dssp::detectSecStructG(Tmdet::VOs::Chain& chain,string key) {
        for( int i=1; i<chain.length-3; i++) {
            if (((any_cast<char>(chain.residues[i].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i].temp[key])=='x')) &&
                ((any_cast<char>(chain.residues[i-1].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i-1].temp[key])=='x')) &&
                checkIfAreOther(chain,Tmdet::Types::SecStructType::G,i,3)) {
                    for(int j=0; j<3; j++) {
                        chain.residues[i+j].ss = Tmdet::Types::SecStructType::G;
                    }
                }
        }
    }

    void Dssp::detectSecStructI(Tmdet::VOs::Chain& chain,string key) {
        for( int i=1; i<chain.length-5; i++) {
            if (((any_cast<char>(chain.residues[i].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i].temp[key])=='x')) &&
                ((any_cast<char>(chain.residues[i-1].temp[key])=='>') || 
                (any_cast<char>(chain.residues[i-1].temp[key])=='x')) &&
                checkIfAreOther(chain,Tmdet::Types::SecStructType::I,i,5)) {
                    for(int j=0; j<5; j++) {
                        chain.residues[i+j].ss = Tmdet::Types::SecStructType::I;
                    }
                }
        }
    }

    bool Dssp::checkIfAreOther(Tmdet::VOs::Chain& chain, Tmdet::Types::SecStruct what,int pos, int r) {
        for( int i=0; i<r; i++) {
            if (chain.residues[pos+i].ss != what and
                chain.residues[pos+i].ss != Tmdet::Types::SecStructType::U) {
                    return false;
                }
        }
        return true;
    }

    void Dssp::detectSecStructT(Tmdet::VOs::Chain& chain) {
        for(int i=4; i<chain.length; i++) {
            if (chain.residues[i].ss == Tmdet::Types::SecStructType::U &&
                (checkIfTurn(chain,i,3,"t3") ||
                checkIfTurn(chain,i,4,"t4") ||
                checkIfTurn(chain,i,5,"t5"))) {
                    chain.residues[i].ss = Tmdet::Types::SecStructType::T;
                }
        }
    }

    bool Dssp::checkIfTurn(Tmdet::VOs::Chain& chain,int pos, int r, string key) {
        for(int i =1; i<r; i++) {
            if (any_cast<char>(chain.residues[pos-i].temp[key]) == '>' || 
                any_cast<char>(chain.residues[pos-i].temp[key]) == 'x') {
                    return true;
                }
        }
        return false;
    }

    void Dssp::detectSecStructS(Tmdet::VOs::Chain& chain) {
        for(int i=2; i<chain.length-2; i++) {
            if (chain.residues[i].ss.code == Tmdet::Types::SecStructs.at('-').code) {
                auto prev_ca_atom = chain.residues[i-2].gemmi.get_ca();
                auto this_ca_atom = chain.residues[i].gemmi.get_ca();
                auto next_ca_atom = chain.residues[i+2].gemmi.get_ca();
                if (prev_ca_atom && this_ca_atom && next_ca_atom) {
                    auto u = gemmi::Vec3(
                        this_ca_atom->pos.x - prev_ca_atom->pos.x,
                        this_ca_atom->pos.y - prev_ca_atom->pos.y,
                        this_ca_atom->pos.z - prev_ca_atom->pos.z
                    );
                    auto v = gemmi::Vec3(
                        next_ca_atom->pos.x - this_ca_atom->pos.x,
                        next_ca_atom->pos.y - this_ca_atom->pos.y,
                        next_ca_atom->pos.z - this_ca_atom->pos.z
                    );
                    double q1 = u.x*u.x + u.y*u.y + u.z*u.z;
                    double q2 = v.x*v.x + v.y*v.y + v.z*v.z;
                    double q = u.x*v.x + u.y*v.y + u.z*v.z;
                    double x=q1*q2;
                    double ckap;
                    if (x>0) {
                        ckap = q/sqrt(x);
                    }
                    else {
                        ckap=0;
                    }
                    double skap = sqrt(1-ckap*ckap);
                    double kap = 180.0 * atan2(skap,ckap) / M_PI;
                    if (kap>70.5) {
                        chain.residues[i].ss = Tmdet::Types::SecStructType::S;
                    }
                }
            }
        }
    }

    void Dssp::detectSecStructBE(Tmdet::VOs::Chain& chain) {
        for(int i=1; i<chain.length-1; i++) {
            if (chain.orderDistance(i-1,i) == 1 && chain.orderDistance(i,i+1) == 1) {
                if (any_cast<int>(chain.residues[i].temp["pb"]) >= 0) {
                    if (any_cast<int>(chain.residues[i-1].temp["pb"]) >=0 || 
                        any_cast<int>(chain.residues[i+1].temp["pb"]) >= 0) {
                        chain.residues[i].ss = Tmdet::Types::SecStructType::E;
                    }
                    else {
                        chain.residues[i].ss = Tmdet::Types::SecStructType::B;
                    }
                }
                if (any_cast<int>(chain.residues[i].temp["apb"]) > 0) {
                    if (any_cast<int>(chain.residues[i-1].temp["apb"]) >= 0 || 
                        any_cast<int>(chain.residues[i+1].temp["apb"]) >= 0) {
                        chain.residues[i].ss = Tmdet::Types::SecStructType::E;
                    }
                    else {
                        chain.residues[i].ss = Tmdet::Types::SecStructType::B;
                    }
                }
            }
        }
        for(int i=1; i<chain.length-2; i++) {
            if (chain.orderDistance(i,i+1) == 1 && chain.orderDistance(i+1,i+2) == 1) {
                if ((chain.residues[i].ss == Tmdet::Types::SecStructType::B ||
                    chain.residues[i].ss == Tmdet::Types::SecStructType::E) &&
                    chain.residues[i+1].ss == Tmdet::Types::SecStructType::U &&
                    (chain.residues[i+2].ss == Tmdet::Types::SecStructType::B ||
                    chain.residues[i+2].ss == Tmdet::Types::SecStructType::E)) {
                        chain.residues[i].ss = 
                        chain.residues[i+1].ss = 
                        chain.residues[i+2].ss = Tmdet::Types::SecStructType::E;
                    }
            }
        }
    }
}
