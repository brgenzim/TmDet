// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#pragma once

#include <array>
#include <any>
#include <gemmi/model.hpp>
#include <VOs/Protein.hpp>

#define VDW(a) (any_cast<double>(a.temp["vdw"]))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/**
 * @brief namespace for tmdet utils
 */
namespace Tmdet::Utils {

    /**
     * @brief temporary geometry data for neighboring atom
     */
    struct surfNeighbor {
        Tmdet::VOs::Atom atom;
        double d;
        double d2;
        double beta;
    };

    /**
     * @brief temporary data for neighbors
     */
    struct surfTemp {
        std::vector<surfNeighbor> neighbors;
        std::vector<double> arc1;
        std::vector<double> arc2;
        std::vector<int> sorted;
    };

    /**
     * @brief the bounding box containing the molecule
     */
    struct boundingBox {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        int l1;
        int l2;
        int l3;
        int op;
        std::vector<std::vector<gemmi::Position>> frame;
        std::vector<std::vector<Tmdet::VOs::Atom *>> closestAtoms;
        std::vector<double> closestDist;
    };

    /**
     * @brief surface cache data
     */
    class surfaceCache {
        private:
            /**
            * @brief the cache container
            */
            std::vector<double> cache;

            /**
            * @brief Convert cache data back to protein value object
            * 
            * @param protein 
            */
            void proteinFromCache(Tmdet::VOs::Protein& protein);

            /**
             * @brief Convert cache data back to chain value object
             * 
             * @param chain 
             * @param index 
             */
            void chainFromCache(Tmdet::VOs::Chain& chain, unsigned int& index);

            /**
             * @brief Convert protein value object to cache data
             * 
             * @param protein 
             */
            void proteinToCache(const Tmdet::VOs::Protein& protein);
        
            /**
             * @brief Convert chain value object to cache data
             * 
             * @param chain 
             */
            void chainToCache(const Tmdet::VOs::Chain& chain);

        public:

            /**
            * @brief Read cache data from file and store it in protein value object
            * 
            * @param protein 
            * @return bool
            */
            bool read(Tmdet::VOs::Protein& protein);

            /**
            * @brief Write cache data from protein value objects to file
            * 
            * @param protein 
            */
            void write(const Tmdet::VOs::Protein& protein);
    };

    /**
     * @brief class for calculating solvent accessible
     *        surface of protein
     */
    class Surface {
        private:

            /**
             * @brief protein value objects containing the structure
             */
            Tmdet::VOs::Protein& protein;

            /**
             * @brief flag for do not use cache
             */
            bool noCache = false;

            /**
             * @brief initialize temporary datat containers
             */
            void initTempData();

            /**
             * @brief Set contacts and summarize surface
             */
            void setContacts();

            /**
             * @brief Set contacts for one atom
             * 
             * @param a_atom 
             */
            void setContactsOfAtom(Tmdet::VOs::Atom& a_atom);

            /**
             * @brief Set the neighbors of atoms
             * 
             * @param a_atom 
             * @param b_atom 
             * @param st 
             */
            void setNeighbor(const Tmdet::VOs::Atom& a_atom, const Tmdet::VOs::Atom& b_atom, surfTemp& st);

            /**
             * @brief Calculate surface of one atom
             * 
             * @param atom 
             * @param st 
             */
            void calcSurfaceOfAtom(Tmdet::VOs::Atom& atom,  surfTemp& st);

            /**
             * @brief Calculate arcs of overlaping atoms in the given plane
             * 
             * @param a_atom 
             * @param st 
             * @param z 
             * @return true 
             * @return false 
             */
            bool calcArcsOfAtom(Tmdet::VOs::Atom& a_atom, surfTemp& st, double z);

            /**
             * @brief Sum up overlaping arcs
             * 
             * @param atom 
             * @param st 
             * @param ss 
             * @return double 
             */
            double calcSumArcsOfAtom(const Tmdet::VOs::Atom& atom, surfTemp& st, bool ss) const;

            /**
             * @brief Set up parameters of the bounding box object
             * 
             * @param box 
             */
            void setBoundingBox(boundingBox& box) const;

            /**
             * @brief Initialize frame for setting outside atoms
             * 
             * @param box 
             */
            void initFrame(boundingBox& box);

            /**
             * @brief Set up frame for setting outside atoms
             * 
             * @param box 
             */
            void setFrame(boundingBox& box);

            /**
             * @brief Smoothing frame data
             * 
             * @param box 
             */
            void smoothFrame(boundingBox& box);

            /**
             * @brief Find the closest atom to the frame
             * 
             * @param box 
             */
            void findClosestAtoms(boundingBox& box);

            /**
             * @brief Run the solvent accessible surface calculation
             */
            void run();

            /**
             * @brief Calculate outside surface (i.e. the surface that is accessible
             *        from outside, needed for beta barrels)
             */
            void setOutsideSurface();

        public:

            /**
             * @brief Construct a new Surface object
             * 
             * @param protein
             */
            explicit Surface(Tmdet::VOs::Protein& protein, bool noCache) : 
                protein(protein),
                noCache(noCache) {
                    run();
            }
            
            /**
             * @brief Destroy the Surface object
             */
            ~Surface()=default;

            std::string toString();
    };
}
