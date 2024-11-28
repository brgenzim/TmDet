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
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                std::string ts = "";
                chain.eachSelectedResidue(
                    [&](Tmdet::ValueObjects::Residue& residue) -> void {
                        ts += any_cast<Tmdet::Types::Region>(residue.temp.at(what)).code;
                    }
                );
                ts += "\n";
                ret += std::format("Chain: {} {}",chain.id, ts);
            }
        );
        return ret;
    }

    bool RegionHandler::getNext(Tmdet::ValueObjects::Chain& chain, int& begin, int& end, std::string what) const {
        return (getNextDefined(chain, begin, what) && getNextSame(chain, begin, end, what));
    }

    bool RegionHandler::getNextDefined(Tmdet::ValueObjects::Chain& chain, int& begin, std::string what) const {
        while(begin < (int)chain.residues.size() && !chain.residues[begin].temp.contains(what)) {
            begin++;
        }
        return (begin < (int)chain.residues.size());
    }

    bool RegionHandler::getNextSame(Tmdet::ValueObjects::Chain& chain, const int& begin, int& end, std::string what) const {
        end = begin;
        while(end < (int)chain.residues.size() && chain.residues[end].temp.contains(what)
            && std::any_cast<Tmdet::Types::Region>(chain.residues[begin].temp.at(what)) 
                        == std::any_cast<Tmdet::Types::Region>(chain.residues[end].temp.at(what))) {
            end++;
        }
        DEBUG_LOG("getNextSame: {} {} {}",chain.id,begin,end-1);
        return true;
    }

    void RegionHandler::replace(Tmdet::ValueObjects::Chain& chain, int beg, int end, Tmdet::Types::Region regionType, std::string what, bool check, Tmdet::Types::Region checkType) {
        DEBUG_LOG("Processing RegionHandler::replace({} {} {} --> {})",
            chain.id,beg,end,regionType.name);
        for (int i = beg; i<= end; i++) {
            if (!check || (check && std::any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at(what)) == checkType)) {
                chain.residues[i].temp.at(what) = std::any(regionType);
            }
        }
        DEBUG_LOG(" Processed RegionHandler::replace()");
    }

    void RegionHandler::finalize() {
        DEBUG_LOG("Processing regionHandler::finalize()");
        //todo rethink this function
        /*protein.eachSelectedChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                for(int i=1; i<chain.length-1; i++) {
                    if (std::any_cast<Tmdet::Types::Region>(chain.residues[i].temp.at("type")) == Tmdet::Types::RegionType::MEMB) {
                        int end=i;
                        getNext(chain,i,end); end--;
                        auto prev = std::any_cast<Tmdet::Types::Region>(chain.residues[i-1].temp.at("type"));
                        auto next = std::any_cast<Tmdet::Types::Region>(chain.residues[end+1].temp.at("type"));
                        DEBUG_LOG("Finalize MEMB: {}({}) - {}({})",i,prev.code,end,next.code);
                        if (prev.isAnnotatedMembraneType() && next.isAnnotatedMembraneType()) {
                            //1, 2 or 3 between two helices
                            auto type = Tmdet::Engine::SiteDetector::getSideByZ(chain.residues[i].getCa()->pos.z * 2);
                            replaceRegion(chain,i,end,type);
                            DEBUG_LOG("\t\t==> {}",type.code);
                        }
                        else if (prev.isNotMembrane() && next.isAnnotatedMembraneType()) {
                            //set type of next
                            replaceRegion(chain,i,end,next);
                            DEBUG_LOG("\t\t==> {}",next.code);
                        }
                        else if (prev.isAnnotatedMembraneType() && next.isNotMembrane()) {
                            //set type of prev
                            replaceRegion(chain,i,end,prev);
                            DEBUG_LOG("\t\t==> {}",prev.code);
                        }
                        else if (prev.isNotMembrane() && next.isNotMembrane()) {
                            //alpha or beta
                            auto type = (chain.type==Tmdet::Types::ChainType::ALPHA?
                                            (end-i>12?Tmdet::Types::RegionType::HELIX:prev):
                                            (end-i>6?Tmdet::Types::RegionType::BETA:prev)
                                        );
                            
                            replaceRegion(chain,i,end,type);
                            DEBUG_LOG("\t\t==> {}",type.code);
                        }
                    }
                }
            }
        );*/
        DEBUG_LOG(" Processed RegionHandler::finalize()");
    }

    void RegionHandler::store() {
        DEBUG_LOG("Processing: RegionHandler::store()");
        protein.eachSelectedChain(
            [&](Tmdet::ValueObjects::Chain& chain) -> void {
                int begin = 0;
                int end = 0;
                while(getNext(chain,begin,end,"type")) {
                    if (end - begin > 0) {
                        Tmdet::ValueObjects::Region region = {
                            {chain.residues[begin].authId, chain.residues[begin].authIcode,chain.residues[begin].labelId},
                            {chain.residues[end-1].authId, chain.residues[end-1].authIcode,chain.residues[end-1].labelId},
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
            }
        );
        DEBUG_LOG(" Processed: RegionHandler::store()");
    }
}