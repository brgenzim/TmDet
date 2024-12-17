// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <any>
#include <Config.hpp>
#include <Helpers/Vector.hpp>
#include <Engine/RegionHandler.hpp>
#include <Engine/SideDetector.hpp>
#include <System/Logger.hpp>
#include <Types/Chain.hpp>
#include <Types/Region.hpp>

/**
 * @brief namespace for tmdet engine
 *
 * @namespace Tmdet
 * @namespace Engine
 */
namespace Tmdet::Engine {

    class RegionHandler {
        private:
            /**
             * @brief structure and tmdet data containing protein value object
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief get next residue that has defined region
             * @param chain 
             * @param begin 
             * @param what 
             * @return bool
             */
            bool getNextDefined(Tmdet::VOs::Chain& chain, int& begin, std::string what) const;

            /**
             * @brief get next region that has the same type than the first one
             * @param chain 
             * @param begin 
             * @param end 
             * @param what 
             * @return bool
             */
            template <typename T>
            bool getNextSame(Tmdet::VOs::Chain& chain, const int& begin, int& end, std::string what) const;

        public:
            /**
             * @brief Construct a new Region Handler object
             * @param protein 
             */
            explicit RegionHandler(Tmdet::VOs::Protein& protein) :
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
            template <typename T>
            bool getNext(Tmdet::VOs::Chain& chain, int& begin, int& end, std::string what) const;

            /**
             * @brief replace type of a region
             * @param chain 
             * @param beg 
             * @param end 
             * @param regionType 
             * @param check 
             * @param checkType 
             */
            void replace(Tmdet::VOs::Chain& chain, int beg, int end, Tmdet::Types::Region regionType, std::string what = "type", bool check = false, Tmdet::Types::Region checkType = Tmdet::Types::RegionType::MEMB);

            /**
             * @brief store regions from residue type to chain regions
             */
            template <typename T>
            void store();

            /**
             * @brief finalize regions, remediate possible errors
             */
            template <typename T>
            
            /**
             * @brief finalize regions, detect and remediate errors
             * 
             * @return int 
             */
            int finalize();

            /**
             * @brief check if directions are the same for two residues
             * 
             * @param chain 
             * @param pos1 
             * @param pos2 
             * @return true 
             * @return false 
             */
            bool sameDirection(Tmdet::VOs::Chain& chain, int pos1, int pos2);

            /**
             * @brief check if direction are not the same for two residues
             * 
             * @param chain 
             * @param pos1 
             * @param pos2 
             * @return true 
             * @return false 
             */
            bool notSameDirection(Tmdet::VOs::Chain& chain, int pos1, int pos2);
    };
}
