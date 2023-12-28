#ifndef __UNITMP_PDBLIB_UTILS_PDBDIST__
#define __UNITMP_PDBLIB_UTILS_PDBDIST__

#include <vector>
#include <string>
#include <any>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>

#define DISTTEMP(obj) (std::any_cast<distTempRes>((obj).any.at(distTempIdx)))

namespace UniTmp::PdbLib::Utils {
    
    const std::string distTempIdx = "distTemp";

    struct distTempRes {
        gemmi::Residue *residue; //parent residue
        std::vector<std::pair<gemmi::Residue,double>> neighbours;  //list of residues close to this residue with the distance

        explicit distTempRes(gemmi::Residue *res) :
            residue(res) {};
    };

    struct distTempAtom {
        gemmi::Atom *atom; //parent atom
        std::vector<std::pair<gemmi::Atom,double>> neighbours;  //list of atoms close to this atom with the distance

        explicit distTempAtom(gemmi::Atom *_atom) :
            atom(_atom) {};
    };

    using boxType = std::vector<std::vector<std::vector<std::vector<std::pair<gemmi::Residue, gemmi::Vec3 >>>>>;

    /**
     * This class is for rapid determination of contacting atoms for surface
     * calculation and/or determining residues close to each other for dssp 
     * and surface calulation algorithms
     */
    class pdbDistance {
        private:
            /**
             * @var gemmi::Structure& _pdb : the protein structure
             */
            gemmi::Structure& _pdb;
            
            /**
             * @var gemmi::Vec3 _min : box lower corner
             */
            gemmi::Vec3 _min = gemmi::Vec3(10000,10000,10000);
            
            /**
             * @var gemmi::Vec3 _max : box upper corner
             */
            gemmi::Vec3 _max = gemmi::Vec3(-10000,-10000,-10000);
            
            /**
             * @var gemmi::Vec3 _unit : box unit
             */
            gemmi::Vec3 _unit;
            
            /**
             * @var boxType _box : the box used to store CA atoms in space
             */
            boxType _box;
            
            /**
             * @var double _maxCaDist : the distance limit used to define
             *                          two C alpha atoms as neighbours
             */
            double _maxCaDist;
            
            /**
             * allocate the temporary structure for each residue and atom
             * @return void
             */
            void _createTemp();
            
            /**
             * delete the temporary structs
             * @return void
             */
            void _deleteTemp();
            
            /**
             * calculate the box size and unit that can incorporate the molecule
             * @return void
             */
            void _unitBox();

            /**
             * allocate memory for the box that can incorporate the molecule
             * @return void
             */
            void _createBox();
            
            /**
             * link C alpha atoms to the appropriate element of the box
             * @return void
             */
            void _fillBox();

            /**
             * iterate through connected boxes for searching for neighbours
             * @param const gemmi::Atom *CA : the C alpha atom
             * @param disTempRes& dtr : the temporary disTempRes struct linked to the residue
             * @return void
             */
            void _getBoxesForCaNeighbours(const gemmi::Atom* CA, distTempRes& dtr);
            
            /**
             * iterate through C alpha atoms in the given box
             * @param gemmi::Atom *CA : the C alpha atom
             * @param disTempRes& dtr : the temporary disTempRes struct linked to the residue
             * @param std::pair<gemmi::Residue, gemmi::Vec3> boxElem : one elem of _box
             * @return void
             */
            void _searchForCaNeighbours(const gemmi::Atom* CA, distTempRes& dtr,
                    std::vector<std::pair<gemmi::Residue, gemmi::Vec3>> boxElem) const;
    

        public:
            pdbDistance(gemmi::Structure& pdb, double maxCaDist);
            
            /**
             * search for C alpha atoms closer than the given limit (setted in constructor)
             * @return void
             */
            void setNeighboursForResiduesByCa();

            /**
             * list atoms in the boxes
             * @return void
             */
            void listBoxElements();

            /**
             * list C alpha neighbours
             * @return void
             */
            void listCaNeighbours();
            
            /**
             * set interacting atom pairs using the calculates C alpha pairs
             * @return void
             */
            void setNeighBoursForAtoms();
    };
    
}
#endif