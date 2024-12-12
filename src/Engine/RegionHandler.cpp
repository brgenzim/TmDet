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

    bool RegionHandler::getNext(Tmdet::VOs::Chain& chain, int& begin, int& end, std::string what) const {
        return (getNextDefined(chain, begin, what) && getNextSame(chain, begin, end, what));
    }

    bool RegionHandler::getNextDefined(Tmdet::VOs::Chain& chain, int& begin, std::string what) const {
        while(begin < (int)chain.residues.size() && !chain.residues[begin].temp.contains(what)) {
            begin++;
        }
        return (begin < (int)chain.residues.size());
    }

    bool RegionHandler::getNextSame(Tmdet::VOs::Chain& chain, const int& begin, int& end, std::string what) const {
        end = begin;
        while(end < (int)chain.residues.size() && chain.residues[end].temp.contains(what)
            && std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at(what)) 
                        == std::any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at(what))) {
            end++;
        }
        return true;
    }

    void RegionHandler::replace(Tmdet::VOs::Chain& chain, int beg, int end, Tmdet::Types::Region regionType, std::string what, bool check, Tmdet::Types::Region checkType) {
        DEBUG_LOG("Processing RegionHandler::replace({} {} {} {} --> {})",
            chain.id,chain.residues[beg].authId,chain.residues[end].authId,
            any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at(what)).name,
            regionType.name);
        for (int i = beg; i<= end; i++) {
            if (!check || (check && std::any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at(what)) == checkType)) {
                chain.residues[i].temp.at(what) = std::any(regionType);
            }
        }
        DEBUG_LOG(" Processed RegionHandler::replace()");
    }

    bool RegionHandler::sameDirection(Tmdet::VOs::Chain& chain, int p1, int p2) {
        if (!chain.residues[p1].temp.contains("direction") || !chain.residues[p2].temp.contains("direction")) {
            return false;
        }
        double d1 = any_cast<double>(chain.residues[p1].temp.at("direction"));
        double d2 = any_cast<double>(chain.residues[p2].temp.at("direction"));
        DEBUG_LOG("sameDirection: {} {} {} {} {}",chain.id,chain.residues[p1].authId,chain.residues[p2].authId,d1,d2);
        return (d1*d2>0 && std::abs(d1)>8 && std::abs(d2)>8);
    }

    bool RegionHandler::notSameDirection(Tmdet::VOs::Chain& chain, int p1, int p2) {
        if (!chain.residues[p1].temp.contains("direction") || !chain.residues[p2].temp.contains("direction")) {
            return false;
        }
        double d1 = any_cast<double>(chain.residues[p1].temp.at("direction"));
        double d2 = any_cast<double>(chain.residues[p2].temp.at("direction"));
        return (d1*d2<0);
    }

    int RegionHandler::finalize() {
        DEBUG_LOG("Processing regionHandler::finalize()");
        int ret = 0;
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int beg = 0;
                int end = 0;
                while(getNext(chain,beg,end,"type")) {
                    auto begType = any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type"));
                    DEBUG_LOG("Finalize checks: {}:{}-{}:{}",chain.id,chain.residues[beg].authId,chain.residues[end-1].authId,begType.code);
                    if ((begType.isAnnotatedTransMembraneType() 
                            || begType.isNotAnnotatedMembrane())
                        && beg>0
                        && end<chain.length-1
                        && chain.residues[beg-1].selected
                        && chain.residues[end].selected
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isNotMembrane()
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).code ==
                            any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).code) {
                        replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("ztype")));
                    }
                    if ((begType.isAnnotatedTransMembraneType()
                            || begType.isNotAnnotatedMembrane())
                        && (beg==0 || end == chain.length-1)
                        && end-beg < (begType.isBeta()?5:10)
                        && chain.residues[beg].selected
                        && chain.residues[end].selected ) {
                        replace(chain,beg,end-1,(beg==0?any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("ztype")):
                                                    any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype"))));
                    }
                    if (begType.isMembraneInside()
                        && end - beg < 6
                        && beg>0
                        && end<chain.length-1
                        && chain.residues[beg-1].selected
                        && chain.residues[end].selected
                        && any_cast<Tmdet::Types::Region>(chain.residues[beg-1].temp.at("type")).isBeta()
                        && any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at("type")).isBeta()) {
                        replace(chain,beg,end-1,Tmdet::Types::RegionType::BETA);
                    }
                    if (any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("type")).isNotAnnotatedMembrane()) {
                        replace(chain,beg,end-1,any_cast<Tmdet::Types::Region>(chain.residues[beg].temp.at("ztype")));
                        if (end-beg>4) {
                            ret += (end-beg);
                        }
                    }
                    beg = end;
                }
            }
        );
        return ret;
        DEBUG_LOG(" Processed RegionHandler::finalize({})",ret);
    }

    void RegionHandler::store() {
        DEBUG_LOG("Processing: RegionHandler::store()");
        protein.eachSelectedChain(
            [&](Tmdet::VOs::Chain& chain) -> void {
                int begin = 0;
                int end = 0;
                while(getNext(chain,begin,end,"type")) {
                    if (end - begin > 0) {
                        Tmdet::VOs::Region region = {
                            {chain.residues[begin].authId, chain.residues[begin].authIcode,chain.residues[begin].labelId,begin},
                            {chain.residues[end-1].authId, chain.residues[end-1].authIcode,chain.residues[end-1].labelId,end},
                            std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type"))
                        };
                        if (std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at("type")).isAnnotatedTransMembraneType()) {
                            chain.numtm++;
                        }
                        DEBUG_LOG("Region stored: {}:{}-{}-{}",chain.id,region.beg.authId,region.end.authId,region.type.code);
                        chain.regions.push_back(region);
                    }
                    begin = end;
                }
                if (chain.numtm == 0) {
                    chain.type = Tmdet::Types::ChainType::NON_TM;
                }
            }
        );
        DEBUG_LOG(" Processed: RegionHandler::store()");
    }
}