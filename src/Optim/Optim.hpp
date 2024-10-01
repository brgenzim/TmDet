#ifndef __TMDET_UTILS_OPTIM__
#define __TMDET_UTILS_OPTIM__

#include <array>
#include <string>
#include <vector>
#include <gemmi/model.hpp>
#include <ValueObjects/TmdetStruct.hpp>
#include <ValueObjects/Membrane.hpp>

namespace Tmdet::Optim {

    struct _slice {
        int numStraight;
        int numTurn;
        int numCA;
        double surf;
        double voronota;
    };

    class Optim {
        private:
            bool run;
            double min;
            double max;
            std::vector<_slice> slices;
            Tmdet::ValueObjects::Membrane membraneVO;
            Tmdet::ValueObjects::TmdetStruct& tmdetVO;

            void init();
            void end();
            void setDistances();
            void setAtomDistances(Tmdet::ValueObjects::Residue& residue);
            void setBoundaries();
            void sumupSlices();
            void residueToSlice(Tmdet::ValueObjects::Residue& residue);
            double getQValueForSlice(const _slice& s);
            std::vector<double> getQValueForSlices();
            double smoothQValues(std::vector<double> qs);

        public:
            explicit Optim(Tmdet::ValueObjects::TmdetStruct& tmdetVO) : tmdetVO(tmdetVO) {}
            ~Optim()=default;

            double getQValue();


    };
}
#endif