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
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <System/Logger.hpp>
#include <VOs/Protein.hpp>
#include <Utils/MyDssp.hpp>

namespace Tmdet::Utils {

    void MyDssp::exec() {
        setCO();
        setAngle();
        //end();
    }

    void MyDssp::end() {
        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("co");
            }
        );
    }

    void MyDssp::setCO() {
        protein.eachSelectedResidue(
            [](Tmdet::VOs::Residue& residue) -> void {
                auto C = residue.gemmi.get_c();
                auto O = residue.gemmi.find_atom("O",' ');
                auto Ca = residue.getCa();
                if ( C!= nullptr && O != nullptr && Ca != nullptr) {
                    residue.temp.try_emplace("co",std::any(gemmi::Vec3(C->pos-O->pos)));
                    residue.temp.try_emplace("ca",std::any(gemmi::Vec3(Ca->pos)));
                }
            }
        );
    }

    void MyDssp::setAngle() {
        protein.eachSelectedChain(
            [](Tmdet::VOs::Chain& chain) -> void {
                for(int i=0; i<chain.length; i++) {
                    int pi = (i==0?i+1:i-1);
                    int ni = (i==chain.length-1?i-1:i+1);
                    if (chain.residues[pi].temp.contains("co")
                        && chain.residues[i].temp.contains("co")
                        && chain.residues[ni].temp.contains("co")
                        && chain.residues[pi].temp.contains("ca")
                        && chain.residues[i].temp.contains("ca")
                        && chain.residues[ni].temp.contains("ca")) {
                        double pco = Tmdet::Helpers::Vector::angle(
                            any_cast<gemmi::Vec3>(chain.residues[pi].temp["co"]),
                            any_cast<gemmi::Vec3>(chain.residues[i].temp["co"])
                        );
                        double nco = Tmdet::Helpers::Vector::angle(
                            any_cast<gemmi::Vec3>(chain.residues[i].temp["co"]),
                            any_cast<gemmi::Vec3>(chain.residues[ni].temp["co"])
                        );
                        double caa = (i!=pi&&i!=ni?Tmdet::Helpers::Vector::angle(
                            any_cast<gemmi::Vec3>(chain.residues[i].temp["ca"]) - any_cast<gemmi::Vec3>(chain.residues[pi].temp["ca"]),
                            any_cast<gemmi::Vec3>(chain.residues[ni].temp["ca"]) - any_cast<gemmi::Vec3>(chain.residues[i].temp["ca"])
                        ):0.0);
                        //if (pco>100&&nco>100&&(pco>135||nco>135)) {
                        if ((caa<85&&((pco>125||nco>125)||(pco>105&&nco>105))) || caa<55) {
                            chain.residues[i].ss = Tmdet::Types::SecStructType::E;
                        }
                        if (pco<100&&nco<100&&(pco<55||nco<55)) {
                            chain.residues[i].ss = Tmdet::Types::SecStructType::S;
                        }

                        INFO_LOG("Angle: {}:{}:{} {:.2f} {:.2f} {:.2f}",
                            chain.id,chain.residues[i].authId,
                            chain.residues[i].ss.code,pco,nco,caa);
                    }
                }
            }
        );
    }
}
