/*
 * cl_Geometry_AD.hpp
 *
 *  Created on: Oct 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_GEOMETRY_AD_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_GEOMETRY_AD_HPP_

#include "cl_Geometry.hpp"
#include <Sacado.hpp>

namespace xtk
{
class Geometry_AD: public Geometry
{
public:
    Geometry_AD(){}

    virtual
    ~Geometry_AD()
    {
    }

    /*
     * This is the direct ad sensitivity evaulation function. allowing everything
     * to stay as a Sacado DFad until last step where it is cast to matrix
     */
    virtual
    void
    ad_sensitivity_evaluation(moris::size_t const & aRowIndex,
                              moris::Matrix< moris::DDRMat > const & aCoordinates,
                              Sacado::Fad::DFad<moris::real> & aLsVal) const
    {
        MORIS_ASSERT(0,"Automatic differentiation sensitivity evaluation not implemented");
    }


    moris::Matrix< moris::DDRMat >
    evaluate_sensitivity_dphi_dp_with_coordinate( moris::size_t const & aRowIndex,
                                                  moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {

        // compute automatic differentiated derivatives
        Sacado::Fad::DFad<moris::real> tLsValWSensSacado;
        this->ad_sensitivity_evaluation(aRowIndex,aCoordinates,tLsValWSensSacado);

        // cast from sacado to moris
        moris::Matrix< moris::DDRMat > tdphidp;
        cast_fad_derivatives_to_matrix(tLsValWSensSacado, tdphidp);

        return tdphidp;
    }


    void
    cast_fad_derivatives_to_matrix(Sacado::Fad::DFad<moris::real> const & aADValue,
                                   moris::Matrix<moris::DDRMat> & aDerivMatrix) const
    {
        aDerivMatrix.resize(1,aADValue.size());
        for(int i = 0; i < aADValue.size(); i++)
        {
            aDerivMatrix(i) = aADValue.dx(i);
        }

    }
};
}



#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_GEOMETRY_AD_HPP_ */
