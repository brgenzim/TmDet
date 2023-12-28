#include <iostream>
#include <array>
#include <gemmi/model.hpp>
#include "pdbDist.hpp"
#include "pdbDssp.hpp"

namespace UniTmp::PdbLib::Utils {

#define DSSP_Q -27888
#define DSSP_PDB_DL 0.5 /*0.5*/
#define DSSP_HBLOW -9900
#define DSSP_HBHIGH -500 /*-500*/

    void pdbCalcDsspOnStructure(gemmi::Structure& pdb) {
        for (gemmi::Model& model : pdb.models) {
            pdbCalcDsspOnModel(model);
        }
    }

    void pdbCalcDsspOnModel(gemmi::Model& model) {
        for (gemmi::Chain& chain : model.chains) {
            pdbCalcDsspOnChain(chain);
        }
    }

    void pdbCalcDsspOnChain(gemmi::Chain& chain) {
        pdbCreateDsspTemp(chain);
        pdbCreateHydrogenBonds(chain);
        for(int d = 0; d < 3; d++) {
            pdbDetectTurns(chain,d);
        }
    }

    void pdbWriteDsspOnStructure(gemmi::Structure& pdb) {
        for (gemmi::Model& model : pdb.models) {
            pdbWriteDsspOnModel(model);
        }
    }

    void pdbWriteDsspOnModel(gemmi::Model& model) {
        for (gemmi::Chain& chain : model.chains) {
            pdbWriteDsspOnChain(chain);
        }
    }

    void pdbWriteDsspOnChain(gemmi::Chain& chain) {
        for (gemmi::Residue& res : chain.residues) {
            std::cout << res.name << " " << DSSPTEMP(res).to[0] << std::endl;
        }
    }

    void pdbCreateDsspTemp(gemmi::Chain& chain) {
        for (gemmi::Residue& res : chain.residues) {
            res.any.insert({{dsspTempIdx, std::any_cast<dsspTemp>(&res)}});
        }
    }

    void pdbCreateHydrogenBonds(gemmi::Chain& chain) {
        for(auto it = std::next(chain.residues.begin()); it != std::prev(chain.residues.end()); ++it) {
            auto const& res = (*it);
            auto prevRes = (*(std::prev(it)));
            auto nextRes = (*(std::next(it)));
            auto CA = res.get_ca();
            auto CE = prevRes.get_c();
            auto OE = prevRes.get_o();
            auto N = res.get_n();
            if (CA != (gemmi::Atom *)nullptr && CE != (gemmi::Atom *)nullptr && OE != (gemmi::Atom *)nullptr && N != (gemmi::Atom *)nullptr) {
                double dco = CE->pos.dist(OE->pos);
                auto hcoord = gemmi::Vec3(
                    N->pos.x+(CE->pos.x-OE->pos.x)/dco,
                    N->pos.y+(CE->pos.y-OE->pos.y)/dco,
                    N->pos.z+(CE->pos.z-OE->pos.z)/dco
                );
                for(auto const& [otherRes, distance] : DISTTEMP(res).neighbours) {
                    if (abs((int)(otherRes.seqid.num - res.seqid.num)) > 1) {
                        auto C = otherRes.get_c();
                        auto O = otherRes.get_o();
                        if (C != (gemmi::Atom *)nullptr && O != (gemmi::Atom *)nullptr) {
                            double dho = hcoord.dist(O->pos);
                            double dhc = hcoord.dist(C->pos);
                            double dno = N->pos.dist(O->pos);
                            double dnc = N->pos.dist(C->pos);
                            double energy;
                            if (dho<DSSP_PDB_DL||dhc<DSSP_PDB_DL||dno<DSSP_PDB_DL||dnc<DSSP_PDB_DL) {
			                	pdbSetHydrogenBond(otherRes,res,DSSP_HBLOW);
			                }
			                else if ((energy=DSSP_Q/dho-DSSP_Q/dhc+DSSP_Q/dnc-DSSP_Q/dno+0.5)<DSSP_HBHIGH) {
				            	pdbSetHydrogenBond(otherRes,res,energy);
				            }
                        }
                    }
                }
            }
        }
    }

    void pdbSetHydrogenBond(gemmi::Residue donor, gemmi::Residue akceptor, double energy) {
        auto dtmp = DSSPTEMP(donor);
        if (energy<dtmp.energy[0]) {
            dtmp.energy[1] = dtmp.energy[0];
            dtmp.to[1] = dtmp.to[0];
            dtmp.energy[0] = energy;
            dtmp.to[0] = &akceptor;
        }
        if (energy<dtmp.energy[1]) {
            dtmp.energy[1] = energy;
            dtmp.to[1] = &akceptor;
        }
    }

    void pdbDetectTurns(gemmi::Chain& chain, int d) {
        int di = d+3;
        for(auto i=0; i<(int)chain.residues.size()-di; i++) {
            auto dtmp = DSSPTEMP(chain.residues[i]);
            if( (int)(dtmp.to[0]->seqid.num - chain.residues[i].seqid.num) == di ||
                (int)(dtmp.to[1]->seqid.num - chain.residues[i].seqid.num) == di) {

                dtmp.ts[d] = (dtmp.ts[d] == '<'?'x':'>');
                for(int j=1; j<di; j++) {
                    if (DSSPTEMP(chain.residues[i+j]).ts[d] == ' ') {
                        DSSPTEMP(chain.residues[i+j]).ts[d] = '*';
                    }
                }
                DSSPTEMP(chain.residues[i+di]).ts[d] = '<';
            }
        }
        for(auto i=1; i<(int)chain.residues.size()-1; i++) {
            auto dtmpm = DSSPTEMP(chain.residues[i-1]);
            auto dtmp = DSSPTEMP(chain.residues[i]);
            auto dtmpp = DSSPTEMP(chain.residues[i+1]);
            if ((dtmpm.ts[d] != '*' && dtmp.ts[d] == '*' && dtmpp.ts[d] != '*') &&
                (dtmpm.ts[d] == 'x' || dtmpp.ts[d] == 'x')) {
                    if (dtmpm.ts[d] == '>') { dtmp.ts[d] = '>'; }
                    if (dtmpm.ts[d] == '<') { dtmp.ts[d] = '<'; }
                    if (dtmpp.ts[d] == '>') { dtmp.ts[d] = '>'; }
                    if (dtmpp.ts[d] == '<') { dtmp.ts[d] = '<'; }

            }
        }
    }

}
