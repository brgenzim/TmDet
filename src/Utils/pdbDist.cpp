#include <iostream>
#include <vector>
#include <utility>
#include <any>
#include <gemmi/model.hpp>
#include <gemmi/math.hpp>
#include "pdbDist.hpp"

namespace UniTmp::PdbLib::Utils {

    /**
     * Constuctor
     * @param gemmi::Structure& pdb : the protein structrure
     * @param double maxCaDist : maximum of C alpha atoms considered as interacted pair
     */
    pdbDistance::pdbDistance(gemmi::Structure& pdb, double maxCaDist)  : 
        _pdb(pdb),
        _maxCaDist(maxCaDist) {
        _createTemp();
        _unitBox();
        _createBox();
        _fillBox();
    }

    void pdbDistance::_createTemp() {
        for (gemmi::Model& model = _pdb.first_model(); gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                auto dtr = std::any_cast<distTempRes>(&res);
                res.any.insert({{distTempIdx, std::any_cast<distTempRes>(dtr)}});
                for(gemmi::Atom& atom : res.atoms) {
                    auto dta = std::any_cast<distTempAtom>(&atom);
                    atom.any.insert({{distTempIdx, std::any_cast<distTempAtom>(dta)}});
                }
            }
        }
    }

    void pdbDistance::_unitBox() {
        for (gemmi::Model model = _pdb.first_model(); gemmi::Chain chain : model.chains) {
            for (gemmi::Residue res : chain.residues) {
                auto CA = res.get_ca();
                if (_min.x > CA->pos.x) { _min.x = CA->pos.x; }
                if (_min.y > CA->pos.y) { _min.y = CA->pos.y; }
                if (_min.z > CA->pos.z) { _min.z = CA->pos.z; }
                if (_max.x < CA->pos.x) { _max.x = CA->pos.x; }
                if (_max.y < CA->pos.y) { _max.y = CA->pos.y; }
                if (_max.z < CA->pos.z) { _max.z = CA->pos.z; }
            }
        }
        _unit = gemmi::Vec3(
                    ( _max.x - _min.x ) / _maxCaDist,
                    ( _max.y - _min.y ) / _maxCaDist,
                    ( _max.z - _min.z ) / _maxCaDist);
        std::cerr << "unitBox calculated: " << "x=" << _unit.x << " y=" << _unit.y << " z=" << _unit.z << std::endl;
    }

    void pdbDistance::_createBox() {
        _box.clear();
        for(int x=0; x<_unit.x+1; x++) {
            std::vector<std::vector<std::vector<std::pair<gemmi::Residue, gemmi::Vec3>>>> yv;
            for(int y=0; y<_unit.y+1; y++) {
                std::vector<std::vector<std::pair<gemmi::Residue, gemmi::Vec3>>> zv;
                for(int z=0; z<_unit.z+1; z++) {
                    std::vector<std::pair<gemmi::Residue, gemmi::Vec3>> be;
                    zv.push_back(be);
                }
                yv.push_back(zv);
            }
            _box.push_back(yv);
        }
       std::cerr << "box created" << std::endl;
    }

    void pdbDistance::_fillBox() {
        for (gemmi::Model model = _pdb.first_model(); gemmi::Chain chain : model.chains) {
            for (gemmi::Residue res : chain.residues) {
                auto CA = res.get_ca();
                auto p = gemmi::Vec3i(
                    (int)(( CA->pos.x - _min.x ) / _maxCaDist),
                    (int)(( CA->pos.y - _min.y ) / _maxCaDist),
                    (int)(( CA->pos.z - _min.z ) / _maxCaDist));
                _box[p.x][p.y][p.z].emplace_back(res,CA->pos);
            }
        }
        std::cerr << "box filled" << std::endl;
    }

    void pdbDistance::setNeighboursForResiduesByCa() {
        gemmi::Model model = _pdb.first_model();
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                auto CA = res.get_ca();
                auto dtr = any_cast<distTempRes>(res.any.at(distTempIdx));
                _getBoxesForCaNeighbours(CA,dtr);
            }
        }
    }

    void pdbDistance::listBoxElements() {
        for(int xx = 0; xx<_unit.x; xx++ ) {
            for(int yy = 0; yy<_unit.y; yy++ ) {
                for(int zz = 0; zz<_unit.z; zz++ ) {
                    std::cout << xx << ":" << yy << ":" << zz << std::endl;
                    for(auto const& [otherRes, otherPos] : _box[xx][yy][zz]) {
                        std::cout << otherPos.str() << "; ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    void pdbDistance::_getBoxesForCaNeighbours(const gemmi::Atom* CA, distTempRes& dtr) {
        auto p = gemmi::Vec3i(
            (int)(( CA->pos.x - _min.x ) / _maxCaDist),
            (int)(( CA->pos.y - _min.y ) / _maxCaDist),
            (int)(( CA->pos.z - _min.z ) / _maxCaDist));
        for(int xx = (p.x>0?p.x - 1:0); (xx<p.x+2 && xx<_unit.x); xx++ ) {
            for(int yy = (p.y>0?p.y - 1:0); (yy<p.y+2 && yy<_unit.y); yy++ ) {
                for(int zz = (p.z>0?p.z - 1:0); (zz<p.z+2 && zz<_unit.z); zz++ ) {
                    _searchForCaNeighbours(CA, dtr, _box[xx][yy][zz]);
                }
            }
        }
    }

    void pdbDistance::_searchForCaNeighbours(const gemmi::Atom* CA, distTempRes& dtr,
        std::vector<std::pair<gemmi::Residue, gemmi::Vec3>> boxElem) const {
        for(auto [otherRes, otherPos] : boxElem) {
            double dist = CA->pos.dist(otherPos);
            if (dist < _maxCaDist) {
                dtr.neighbours.emplace_back(otherRes, dist);
            }
        }
    }

    void pdbDistance::listCaNeighbours() {
        gemmi::Model model = _pdb.first_model();
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                std::cout << res.name << "---------" <<  std::endl;
                auto dtr = any_cast<distTempRes>(res.any.at(distTempIdx));
                for(const auto& [otherRes, distance] : dtr.neighbours) {
                    std::cout << res.str() << "-" << otherRes.str() << ":" << distance << std::endl;
                }
            }
        }
    }

    void pdbDistance::setNeighBoursForAtoms() {
        gemmi::Model model = _pdb.first_model();
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& res : chain.residues) {
                auto dtr = any_cast<distTempRes>(res.any.at(distTempIdx));
                for(gemmi::Atom& atom : res.atoms) {
                    //TODO iterate through atoms of neighbour residues
                }
            }
        }
    }

        
}
