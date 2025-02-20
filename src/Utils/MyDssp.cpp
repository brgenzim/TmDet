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
#include <DTOs/Dssp.hpp>
#include <Helpers/Vector.hpp>
#include <System/Logger.hpp>
#include <VOs/Protein.hpp>
#include <Utils/MyDssp.hpp>

namespace Tmdet::Utils {

    void MyDssp::exec() {
        setCO();
        setAngle();
        setHelix();
        setExtended();
        //removeMins("S");
        //removeMins("E");
        end();
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                DEBUG_LOG("MYDSSP: {}:{}",chain.id,Tmdet::DTOs::Dssp::getSecondaryStructure(chain));
            }
        );
    }

    void MyDssp::end() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.temp.erase("S");
                residue.temp.erase("E");

            }
        );
    }

    void MyDssp::setCO() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                if (auto Ca = residue.getCa();  Ca != nullptr) {
                    residue.temp.try_emplace("ca",std::any(gemmi::Vec3(Ca->pos)));
                }
                auto C = residue.gemmi.get_c();
                auto O = residue.gemmi.find_atom("O",' ');
                if ( C!= nullptr && O != nullptr) {
                    residue.temp.try_emplace("co",std::any(gemmi::Vec3(C->pos-O->pos)));
                }
                residue.temp.try_emplace("S",std::any(0));
                residue.temp.try_emplace("E",std::any(0));
            }
        );
    }

    void MyDssp::setAngle() {
        protein.eachSelectedChain(
            [](Tmdet::VOs::Chain& chain) -> void {
                for(int i=0; i<chain.length-4; i++) {
                    if (chain.residues[i].temp.contains("ca")
                        && chain.residues[i+1].temp.contains("ca")
                        && chain.residues[i+2].temp.contains("ca")
                        && chain.residues[i+3].temp.contains("ca")) {
                            gemmi::Vec3 ca1 = any_cast<gemmi::Vec3>(chain.residues[i].temp["ca"])
                                                - any_cast<gemmi::Vec3>(chain.residues[i+1].temp["ca"]);
                            gemmi::Vec3 ca2 = any_cast<gemmi::Vec3>(chain.residues[i+1].temp["ca"])
                                                - any_cast<gemmi::Vec3>(chain.residues[i+2].temp["ca"]);
                            gemmi::Vec3 ca3 = any_cast<gemmi::Vec3>(chain.residues[i+2].temp["ca"])
                                                - any_cast<gemmi::Vec3>(chain.residues[i+3].temp["ca"]);
                            gemmi::Vec3 ax1 = ca1.cross(ca2);
                            gemmi::Vec3 ax2 = ca2.cross(ca3);
                            double a1 = Tmdet::Helpers::Vector::angle(ca1,ca2);
                            double a2 = Tmdet::Helpers::Vector::angle(ca2,ca3);
                            double a3 = Tmdet::Helpers::Vector::angle(ax1,ax2);
                            if (std::abs(90-a1) < 13 && std::abs(90-a2) < 13 && std::abs(52-a3) < 15) {
                                chain.residues[i].temp["S"] = std::any(any_cast<int>(chain.residues[i].temp["S"]) + 1);
                                chain.residues[i+1].temp["S"] = std::any(any_cast<int>(chain.residues[i+1].temp["S"]) + 1);
                                chain.residues[i+2].temp["S"] = std::any(any_cast<int>(chain.residues[i+2].temp["S"]) + 1);
                                chain.residues[i+3].temp["S"] = std::any(any_cast<int>(chain.residues[i+3].temp["S"]) + 1);
                            }
                            else if (std::abs(55-a1) < 35 && std::abs(55-a2) < 35 && std::abs(160-a3) < 20) {
                                chain.residues[i].temp["E"] = std::any(any_cast<int>(chain.residues[i].temp["E"]) + 1);
                                chain.residues[i+1].temp["E"] = std::any(any_cast<int>(chain.residues[i+1].temp["E"]) + 1);
                                chain.residues[i+2].temp["E"] = std::any(any_cast<int>(chain.residues[i+2].temp["E"]) + 1);
                                //chain.residues[i+3].temp["E"] = std::any(any_cast<int>(chain.residues[i+3].temp["E"]) + 1);
                            }
                            if (a3<55) {
                                chain.residues[i].temp["E"] = std::any(0);
                            }
                            else if (a3<80 && any_cast<int>(chain.residues[i].temp["E"]) == 1) {
                                chain.residues[i].temp["E"] = std::any(0);
                            }

                        INFO_LOG("Angle: {}:{}:{} {:.2f} {:.2f} {:.2f} {} {}:{}",
                            chain.id,chain.residues[i].authId,
                            chain.residues[i].ss.code,a1,a2,a3,
                            any_cast<int>(chain.residues[i].temp["S"]),
                            any_cast<int>(chain.residues[i].temp["E"]),
                            (any_cast<int>(chain.residues[i].temp["E"])>0?'E':' '));
                    }
                }
            }
        );
    }

    void MyDssp::setExtended() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int numOne = 0;
                for(int i=0; i<chain.length; i++) {
                    if (chain.residues[i].selected) {
                        if (any_cast<int>(chain.residues[i].temp["E"]) == 1) {
                            numOne += 1;
                        }
                        else {   
                            if (numOne > 4) {
                                for (int j=1; j<=numOne; j++) {
                                    chain.residues[i-j].temp["E"] = std::any(0);
                                }
                            }
                            numOne = 0;
                        }
                    }
                }
            }
        );
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.ss = any_cast<int>(residue.temp["E"])>0?
                    Tmdet::Types::SecStructType::E:
                    residue.ss;
            }
        );
    }

    void MyDssp::setHelix() {
        protein.eachSelectedResidue(
            [&](Tmdet::VOs::Residue& residue) -> void {
                residue.ss = (residue.type.a1 == 'P'?
                    Tmdet::Types::SecStructType::U:
                    residue.ss);
                //residue.ss = any_cast<int>(residue.temp["S"])>0?
                //    Tmdet::Types::SecStructType::S:
                //    residue.ss; //Tmdet::Types::SecStructType::U;
            }
        );
    }

    void MyDssp::removeMins(std::string what) {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int i=1;
                while(i<chain.length-1) {
                    int sp = any_cast<int>(chain.residues[i-1].temp[what]);
                    int s = any_cast<int>(chain.residues[i].temp[what]);
                    int min = s;
                    int beg = i;
                    if (s>0 && sp>s) {
                        while(s<=min && i<chain.length-1 && min>0) {
                            int s = any_cast<int>(chain.residues[i].temp[what]);
                            if (min<s) {
                                min=s;
                            }
                            i++;
                        }
                        if (min>0) {
                            for (int j=beg; j<i; j++) {
                                if (any_cast<int>(chain.residues[j].temp[what]) == min) {
                                    chain.residues[j].ss = Tmdet::Types::SecStructType::U;
                                }
                            }
                        }
                    }
                    else {
                        i++;
                    }
                }
            }
        );
    }
}
