#pragma once

#include <any>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/SideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>

namespace Tmdet::Engine {

    class RegionHandler {
        private:
            Tmdet::ValueObjects::Protein& protein;

            /**
             * @brief get next residue that has defined region
             * @param chain 
             * @param begin 
             * @param what 
             * @return bool
             */
            bool getNextDefined(Tmdet::ValueObjects::Chain& chain, int& begin, std::string what) const;

            /**
             * @brief get next region that has the same type than the first one
             * @param chain 
             * @param begin 
             * @param end 
             * @param what 
             * @return bool
             */
            bool getNextSame(Tmdet::ValueObjects::Chain& chain, const int& begin, int& end, std::string what) const;

        public:
            /**
             * @brief Construct a new Region Handler object
             * @param protein 
             */
            explicit RegionHandler(Tmdet::ValueObjects::Protein& protein) :
                protein(protein) {}

            /**
             * @brief Destroy the Region Helper object
             */
            ~RegionHandler()=default;

            /**
             * @brief string representation of regions in the chain (for debugging purpose)
             * @param what
             * @return std::string 
             */
            std::string toString(std::string what);

            /**
             * @brief get next region
             * @param chain 
             * @param begin 
             * @param end 
             * @param what 
             * @return bool
             */
            bool getNext(Tmdet::ValueObjects::Chain& chain, int& begin, int& end, std::string what) const;

            /**
             * @brief replace type of a region
             * @param chain 
             * @param beg 
             * @param end 
             * @param regionType 
             * @param check 
             * @param checkType 
             */
            void replace(Tmdet::ValueObjects::Chain& chain, int beg, int end, Tmdet::Types::Region regionType, std::string what = "type", bool check = false, Tmdet::Types::Region checkType = Tmdet::Types::RegionType::MEMB);

            /**
             * @brief store regions from residue type to chain regions
             */
            void store();

            /**
             * @brief finalize regions, remediate possible errors
             */
            int finalize();

            bool sameDirection(Tmdet::ValueObjects::Chain& chain, int pos1, int pos2);
    };
}
