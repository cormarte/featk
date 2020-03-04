#include <featk/core/featkDefines.h>
#include <featk/material/featkIsotropicLinearElastic3DMaterial.h>


featkIsotropicLinearElastic3DMaterial::featkIsotropicLinearElastic3DMaterial(double youngModulus, double poissonRatio) {

    double l = (youngModulus*poissonRatio)/((1.0+poissonRatio)*(1.0-2.0*poissonRatio));
    double m = youngModulus/(2.0*(1.0+poissonRatio));

    this->constitutiveMatrix = (AttributeValueType<3, 4>() << 2.0*m+l,     0.0,     0.0,     0.0,       l,     0.0,     0.0,     0.0,       l,
                                                                  0.0,       m,     0.0,       m,     0.0,     0.0,     0.0,     0.0,     0.0,
                                                                  0.0,     0.0,       m,     0.0,     0.0,     0.0,       m,     0.0,     0.0,
                                                                  0.0,       m,     0.0,       m,     0.0,     0.0,     0.0,     0.0,     0.0,
                                                                    l,     0.0,     0.0,     0.0, 2.0*m+l,     0.0,     0.0,     0.0,       l,
                                                                  0.0,     0.0,     0.0,     0.0,     0.0,       m,     0.0,       m,     0.0,
                                                                  0.0,     0.0,       m,     0.0,     0.0,     0.0,       m,     0.0,     0.0,
                                                                  0.0,     0.0,     0.0,     0.0,     0.0,       m,     0.0,       m,     0.0,
                                                                    l,     0.0,     0.0,     0.0,       l,     0.0,     0.0,     0.0, 2.0*m+l).finished();
}

featkIsotropicLinearElastic3DMaterial::~featkIsotropicLinearElastic3DMaterial() {

}
