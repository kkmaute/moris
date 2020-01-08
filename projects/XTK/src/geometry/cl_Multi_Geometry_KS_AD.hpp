/*
 * cl_Multi_Geometry_KS_AD.hpp
 *
 *  Created on: Oct 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_KS_AD_HPP_
#define PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_KS_AD_HPP_

#include "cl_Geometry_AD.hpp"

namespace xtk
{
class Multi_Geometry_KS_AD : public xtk::Geometry_AD
{
public:
    Multi_Geometry_KS_AD(){}

    Multi_Geometry_KS_AD(moris::Cell<xtk::Geometry_AD*> const & aGeomVector,
                         moris::real                    aBeta,
                         moris::uint                    aNumDesVar) :
                             mGeometries(aGeomVector),
                             mBeta(aBeta),
                             mNumDesVars(aNumDesVar)
    {
    }

    bool is_analytic() const
    {
        return true;
    }


    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 4;
    }

    moris::uint get_num_des_vars() const
    {
        return mNumDesVars;
    }

    moris::real evaluate_field_value_with_coordinate(moris::size_t const & aRowIndex,
                                                     moris::Matrix< moris::DDRMat > const & aCoordinates) const
    {
        real tSummationTerm = 0;
        Sacado::Fad::DFad<moris::real> tBuffer = 0;

        for(moris::uint i =0; i <mGeometries.size(); i++)
        {
            mGeometries(i)->ad_sensitivity_evaluation(aRowIndex,aCoordinates,tBuffer);

            tSummationTerm = tSummationTerm + std::exp(mBeta*tBuffer.val());
        }

        return (1/mBeta) * std::log(tSummationTerm);
    }

    void
    ad_sensitivity_evaluation(moris::size_t const & aRowIndex,
                              moris::Matrix< moris::DDRMat > const & aCoordinates,
                              Sacado::Fad::DFad<moris::real> & aLsVal) const
    {
        Sacado::Fad::DFad<moris::real> tSummationTerm = 0;
        Sacado::Fad::DFad<moris::real> tBuffer = 0;

        for(moris::uint i =0; i <mGeometries.size(); i++)
        {
            mGeometries(i)->ad_sensitivity_evaluation(aRowIndex,aCoordinates,tBuffer);

            tSummationTerm = tSummationTerm + std::exp(mBeta*tBuffer);
        }

        aLsVal = (1/mBeta) * std::log(tSummationTerm);
    }

private:
    moris::Cell<xtk::Geometry_AD*> mGeometries;
    moris::real mBeta;
    moris::uint mNumDesVars;
};
}


#endif /* PROJECTS_XTK_SRC_GEOMETRY_CL_MULTI_GEOMETRY_KS_AD_HPP_ */
