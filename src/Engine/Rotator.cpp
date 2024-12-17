// © 2003-2024 Gabor E. Tusnady <tusnady.gabor@ttk.hu> and TmDet developer team
//             Protein Bioinformatics Research Group 
//             Research Center of Natural Sciences, HUN-REN
//
// License:    CC-BY-NC-4.0, see LICENSE.txt

#include <cmath>
#include <gemmi/model.hpp>
#include <Config.hpp>
#include <Engine/Rotator.hpp>

namespace Tmdet::Engine {

    Rotator::Rotator() {
        alpha_step = std::stof(environment.get("TMDET_BALL_DIST",DEFAULT_TMDET_BALL_DIST)) / 2;
        logger.debug("Init rotator. alpha_step: {}",alpha_step);
    }

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
        logger.debug("Next normal: alpha {} beta {} normal {} {} {}",
            alpha, beta, normal.x, normal.y, normal.z);
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
            beta_step = std::stof(environment.get("TMDET_BALL_DIST",DEFAULT_TMDET_BALL_DIST))/q;
        }
        else {
            beta_step = 2* M_PI; 
            q = 0;
        }
        logger.debug("nextAlpha: alpha: {} beta_step: {} q: {} qq: {}",
            alpha,beta_step,q,qq);
    }
}
