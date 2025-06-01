// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <any>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/SideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>

namespace Tmdet::Engine {

    std::string RegionHandler::toString(std::string what) {
        std::string ret = "\n";
        ret += what;
        ret += "\n";
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                std::string ts = "";
                chain.eachSelectedResidue(
                    [&](Tmdet::VOs::Residue& residue) -> void {
                        ts += any_cast<Tmdet::Types::Region>(residue.temp.at(what)).code;
                    }
                );
                ts += "\n";
                ret += std::format("Chain: {} {}",chain.id, ts);
            }
        );
        return ret;
    }

    template <typename T>
    bool RegionHandler::getNext(Tmdet::VOs::Chain& chain, int& begin, int& end, std::string what) const {
        return (getNextDefined(chain, begin, what) && getNextSame<T>(chain, begin, end, what));
    }
    template bool RegionHandler::getNext<int>(Tmdet::VOs::Chain& chain, int& begin, int& end, std::string what) const;
    template bool RegionHandler::getNext<Tmdet::Types::Region>(Tmdet::VOs::Chain& chain, int& begin, int& end, std::string what) const;

    bool RegionHandler::getNextDefined(Tmdet::VOs::Chain& chain, int& begin, std::string what) const {
        while(begin < (int)chain.residues.size() && !chain.residues[begin].temp.contains(what)) {
            begin++;
        }
        return (begin < (int)chain.residues.size());
    }

    template <typename T>
    bool RegionHandler::getNextSame(Tmdet::VOs::Chain& chain, const int& begin, int& end, std::string what) const {
        end = begin + 1;
        while(end < (int)chain.residues.size() && chain.residues[end].temp.contains(what)
            && chain.orderDistance(end-1,end) == 1
            && std::any_cast<T>(chain.residues[begin].temp.at(what)) 
                        == std::any_cast<T>(chain.residues[end].temp.at(what))) {
            end++;
        }
        return true;
    }
    template bool RegionHandler::getNextSame<int>(Tmdet::VOs::Chain& chain, const int& begin, int& end, std::string what) const;
    template bool RegionHandler::getNextSame<Tmdet::Types::Region>(Tmdet::VOs::Chain& chain, const int& begin, int& end, std::string what) const;

    void RegionHandler::replace(Tmdet::VOs::Chain& chain, int beg, int end, Tmdet::Types::Region regionType, std::string what, bool check, Tmdet::Types::Region checkType) {
        for (int i = beg; i<= end; i++) {
            if (!check || (check && std::any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at(what)) == checkType)) {
                chain.residues[i].temp.at(what) = std::any(regionType);
            }
        }
    }

    template <typename T>
    std::vector<simpleRegion> RegionHandler::getAll(Tmdet::VOs::Chain& chain, std::string what) {
        int beg = 0;
        int end = 0;
        std::vector<simpleRegion> ret;
        while(getNext<T>(chain,beg,end,what)) {
            ret.emplace_back(beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")));
            beg = end;
        }
        return ret;
    }
    template std::vector<simpleRegion> RegionHandler::getAll<Tmdet::Types::Region>(Tmdet::VOs::Chain& chain, std::string what);

    template <typename T>
    int RegionHandler::finalize() {
        int ret = 0;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int beg = 0;
                int end = 0;
                while(getNext<T>(chain,beg,end,"type")) {
                    auto begType = any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type"));
                    // 1111MM..MM111 or 111BB..BB1111 or 111HH..HH111
                    if ((begType.isAnnotatedTransMembraneType() 
                            || begType.isNotAnnotatedMembrane())
                        && beg>0
                        && end<chain.length-1
                        && chain.residues[beg-1].selected
                        && chain.residues[end].selected
                        && !chain.isGapBetween(beg-1,beg)
                        && !chain.isGapBetween(end-1,end)
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isNotMembrane()
                        && any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isNotMembrane()
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("ztype")).code ==
                            any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("ztype")).code) {

                        replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("ztype")));
                    }
                    //short B or H or M at the end of the chain
                    if ((begType.isAnnotatedTransMembraneType()
                            || begType.isNotAnnotatedMembrane())
                        && (beg==0 || end == chain.length)
                        && end-beg < (begType.isBeta()?3:0)
                        && chain.residues[beg].selected
                        && chain.residues[end-1].selected ) {
                        replace(chain,beg,end-1,(beg==0?any_cast<Tmdet::Types::Region>(chain.residues[end-1].temp.at("ztype")):
                            any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype"))));
                    }
                    // any MMM anywhere that not handled so far
                    if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane()) {
                        replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype")));
                        double q = any_cast<double>(chain.residues[beg].temp.at("z")) * any_cast<double>(chain.residues[end-1].temp.at("z"));
                        if (end-beg>4 &&  q < -10
                            && std::abs(any_cast<double>(chain.residues[beg].temp.at("z"))) > 10
                            && std::abs(any_cast<double>(chain.residues[end-1].temp.at("z"))) > 10) {
                            ret += (end-beg);
                        }
                    }
                    beg = end;
                }
            }
        );
        return ret;
    }
    template int RegionHandler::finalize<int>();
    template int RegionHandler::finalize<Tmdet::Types::Region>();

    template <typename T>
    void RegionHandler::store() {
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
               if (chain.signalP[1] > 0) {
                    Tmdet::VOs::Region region = {
                        {chain.residues[0].authId, chain.residues[0].authIcode,chain.residues[0].labelId,0},
                        {chain.residues[chain.signalP[1]-1].authId, chain.residues[chain.signalP[1]-1].authIcode,chain.residues[chain.signalP[1]-1].labelId,chain.signalP[1]-1},
                        Tmdet::Types::RegionType::SIGNAL
                    };
                    chain.regions.push_back(region);
                    DEBUG_LOG(">>Region stored: {}:{}-{}-{}",chain.id,region.beg.authId,region.end.authId,region.type.code);
                }
                chain.numtm = 0;
                int begin = chain.signalP[1];
                int end = 0;
                while(getNext<T>(chain,begin,end,"type")) {
                    if (end - begin > 0) {
                        Tmdet::VOs::Region region = {
                            {chain.residues[begin].authId, chain.residues[begin].authIcode,chain.residues[begin].labelId,begin},
                            {chain.residues[end-1].authId, chain.residues[end-1].authIcode,chain.residues[end-1].labelId,end-1},
                            std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type"))
                        };
                        if (std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type")).isAnnotatedTransMembraneType()) {
                            chain.numtm++;
                        }
                        chain.regions.push_back(region);
                    }
                    begin = end;
                }
                if (chain.numtm == 0) {
                    chain.type = Tmdet::Types::ChainType::NON_TM;
                }
            }
        );
    }

    template void RegionHandler::store<int>();
    template void RegionHandler::store<Tmdet::Types::Region>();
}
