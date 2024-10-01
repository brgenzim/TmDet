#ifndef __TMDET_UTILS_ROTATOR__
#define __TMDET_UTILS_ROTATOR__

#include <vector>
#include <map>
#include <set>

namespace Tmdet::Optim {

    class Rotator {
        private:
            const double BALL_DIST = 1.0;
            double alpha_dist;
            double alpha=0;
            double alpha_step;
            double beta=0;
            double beta_step;
            double q;
            double qq;

            void nextAlpha();
        public:
            explicit Rotator() {
                alpha_dist = BALL_DIST / 2;
            }
            ~Rotator()=default;

            bool next(gemmi::Vec3& normal);
    };
}

#endif
