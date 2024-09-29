#include <cmath>
#include <gemmi/model.hpp>
#include <Utils/Rotator.hpp>

namespace Tmdet::Utils {

    /**
     * @brief rotate a normal vector around the 4PI
     *        give the next rotated vector
     *        return false if it ended
     * 
     * @param normal 
     * @return bool
     */
    bool Rotator::next(gemmi::Vec3& normal) {
        if (alpha>M_PI) {
            return false;
        }
        normal.x = cos(beta) * q;
        normal.y = sin(beta) * q;
        normal.z = qq;
        beta += beta_step;
        if (beta > 2*M_PI) {
            beta=0;
            nextAlpha();
        }
        return true;
    }

    /**
     * @brief calculate next alpha angle
     * 
     */
    void Rotator::nextAlpha() {
        alpha += alpha_step;
        q = sin(alpha);
        qq = cos(alpha);
        if (q>1e-10) {
            beta_step = BALL_DIST/q;
        }
        else {
            beta_step = 2* M_PI; 
            q = 0;
        }
    }
}
