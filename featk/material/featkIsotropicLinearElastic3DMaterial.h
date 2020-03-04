/*==========================================================================

  Program:   Finite Element Analysis Toolkit
  Module:    featkIsotropicLinearElastic3DMaterial.h

  Copyright (c) Corentin Martens
  All rights reserved.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
     NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR
     ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR
     OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING
     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
     OTHER DEALINGS IN THE SOFTWARE.

==========================================================================*/

/**
 *
 * @class featkIsotropicLinearElastic3DMaterial
 *
 * @brief Hookean linear elastic material in 3D.
 *
 * featkIsotropicLinearElastic3DMaterial is a 3D linear elastic
 * featkMaterialBase governed by Hooke's law.
 * featkIsotropicLinearElastic3DMaterial generates and stores the
 * non-symmetrized 4th-order Hooke's tensor (i.e. a
 * Eigen::Matrix<double, 9, 9>), given its Young's modulus \f$E\f$ and
 * Poisson's ratio \f$\nu\f$. For generalization purpose of the solving
 * routines, the Hooke's law is not explicitly symmetrized in featk, hence
 * Hooke's tensor is re-arranged as follows:
 *
 * \f[
 *
 * \begin{bmatrix}
 * \sigma_{xx}\\
 * \sigma_{xy}\\
 * \sigma_{xz}\\
 * \sigma_{yx}\\
 * \sigma_{yy}\\
 * \sigma_{yz}\\
 * \sigma_{zx}\\
 * \sigma_{zy}\\
 * \sigma_{zz}\\
 * \end{bmatrix}
 *
 * =
 *
 * \begin{bmatrix}
 * 2\mu+\lambda & 0 & 0 & 0 & \lambda & 0 & 0 & 0 & \lambda\\
 * 0 & \mu & 0 & \mu & 0 & 0 & 0 & 0 & 0\\
 * 0 & 0 & \mu & 0 & 0 & 0 & \mu & 0 & 0\\
 * 0 & \mu & 0 & \mu & 0 & 0 & 0 & 0 & 0\\
 * \lambda & 0 & 0 & 0 & 2\mu+\lambda & 0 & 0 & 0 & \lambda\\
 * 0 & 0 & 0 & 0 & 0 & \mu & 0 & \mu & 0\\
 * 0 & 0 & \mu & 0 & 0 & 0 & \mu & 0 & 0\\
 * 0 & 0 & 0 & 0 & 0 & \mu & 0 & \mu & 0\\
 * \lambda & 0 & 0 & 0 & \lambda & 0 & 0 & 0 & 2\mu+\lambda
 * \end{bmatrix}
 *
 * \begin{bmatrix}
 * u_{x,x}\\
 * u_{x,y}\\
 * u_{x,z}\\
 * u_{y,x}\\
 * u_{y,y}\\
 * u_{y,z}\\
 * u_{z,x}\\
 * u_{z,y}\\
 * u_{z,z}\\
 * \end{bmatrix} \quad \lambda = \frac{\nu E}{(1+\nu)(1-2\nu)}, \quad \mu=\frac{E}{2(1+\nu)}
 *
 * \f]
 *
 */

#ifndef FEATKISOTROPICLINEARELASTIC3DMATERIAL_H
#define FEATKISOTROPICLINEARELASTIC3DMATERIAL_H

#include <featk/material/featkMaterialBase.h>

class featkIsotropicLinearElastic3DMaterial final : public featkMaterialBase<3> {

    public:

        featkIsotropicLinearElastic3DMaterial(double youngModulus, double poissonRatio);
        ~featkIsotropicLinearElastic3DMaterial();
};

#endif // FEATKISOTROPICLINEARELASTIC3DMATERIAL_H
