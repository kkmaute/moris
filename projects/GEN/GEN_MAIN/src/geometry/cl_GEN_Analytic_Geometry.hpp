/*
 * cl_GEN_Analytic_Geometry.hpp
 *
 *  Created on: Feb 18, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_ANALYTIC_GEOMETRY_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_ANALYTIC_GEOMETRY_HPP_

#include "cl_GEN_Geometry.hpp"


namespace moris
{
namespace ge
{
    typedef std::function< void (       moris::real         & aReturnValue,
                                  const Matrix< DDRMat >    & aPoint,
                                  const moris::Cell< real > & aConstant ) > AnalyticFunction;
//------------------------------------------------------------------------------

    class Analytic_Geometry : public GEN_Geometry
    {
    private:
        AnalyticFunction mAnalyticFunction;
        AnalyticFunction mAnalyticDerFunction;

    //------------------------------------------------------------------------------
    public:
        Analytic_Geometry( AnalyticFunction aAnalyticFunction )
        : mAnalyticFunction( aAnalyticFunction )
        {
        }

        Analytic_Geometry( AnalyticFunction aAnalyticFunction,
                           AnalyticFunction aAnalyticDerFunction )
        : mAnalyticFunction( aAnalyticFunction ),
          mAnalyticDerFunction( aAnalyticDerFunction )
        {
        }

        void set_analytic_function( AnalyticFunction aAnalyticFunction )
        {
            mAnalyticFunction = aAnalyticFunction;
        }

        //------------------------------------------------------------------------------
        bool is_analytic() const
        {
            return true;
        }
        //------------------------------------------------------------------------------
        void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
        {
            MORIS_ERROR( false,"Analytic_Geometry::get_dphi_dp_size() - not implemented" );
        }
        //------------------------------------------------------------------------------
        void eval( moris::real                      & aReturnValue,
                   const moris::Matrix< DDRMat >    & aPoint,
                   const moris::Cell< moris::real > & aConst )
        {
            mAnalyticFunction( aReturnValue, aPoint, aConst );
        }
        //------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
}   /* end ge namespace */
}   /* end moris namespace */


#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMETRY_CL_GEN_ANALYTIC_GEOMETRY_HPP_ */
