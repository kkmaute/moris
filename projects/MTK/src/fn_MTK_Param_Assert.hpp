#ifndef MORIS_CL_MTK_PARAM_ASSERT_HPP_
#define MORIS_CL_MTK_PARAM_ASSERT_HPP_

#include "cl_Matrix.hpp"    // moris::Matrix, DDRMat  (LINALG)
#include "fn_assert.hpp"    // MORIS_ASSERT           (MRS/ASR)

namespace moris::mtk
{
    /**
     * @brief Checks that the parametric coordinates lie within a quad with local coordinates
     * 
     * @param aParam parametric coordinates of a point in 3D space (2D parametric). size: 2x1
     * @param aContext Error message printing message to determine where the check occurred
     */
    static void
    assert_param_in_bounds_box(
            const Matrix< DDRMat >& aParam,
            const char*             aContext )
    {
#ifdef MORIS_HAVE_DEBUG
        if ( aParam.n_rows() >= 1 )
        {
            MORIS_ASSERT(
                    aParam( 0 ) >= -1.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                    "MTK parametric coordinate out of bounds (param[0]=%e, context=%s)",
                    aParam( 0 ),
                    aContext );
        }
        if ( aParam.n_rows() >= 2 )
        {
            MORIS_ASSERT(
                    aParam( 1 ) >= -1.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK parametric coordinate out of bounds (param[1]=%e, context=%s)",
                    aParam( 1 ),
                    aContext );
        }
#endif
    }

        /**
     * @brief Checks that the parametric coordinates lie within a triangle with local coordinates
     * 
     * @param aParam parametric coordinates of a point in 3D space (2D parametric). size: 2x1
     * @param aContext Error message printing message to determine where the check occurred
     */
    static void
    assert_param_in_bounds_simplex(
            const Matrix< DDRMat >& aParam,
            const char*             aContext )
    {
#ifdef MORIS_HAVE_DEBUG
        if ( aParam.n_rows() >= 1 )
        {
            MORIS_ASSERT(
                    aParam( 0 ) >= 0.0 - 1e-12 && aParam( 0 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (eta=%e, context=%s)",
                    aParam( 0 ),
                    aContext );
        }
        if ( aParam.n_rows() >= 2 )
        {
            MORIS_ASSERT(
                    aParam( 1 ) >= 0.0 - 1e-12 && aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (zeta=%e, context=%s)",
                    aParam( 1 ),
                    aContext );
            MORIS_ASSERT(
                    aParam( 0 ) + aParam( 1 ) <= 1.0 + 1e-12,
                    "MTK TRI parametric coordinate out of bounds (eta+zeta=%e, context=%s)",
                    aParam( 0 ) + aParam( 1 ),
                    aContext );
        }
#endif
    }
}    // namespace moris::mtk

#endif