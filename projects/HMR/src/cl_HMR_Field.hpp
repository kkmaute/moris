/*
 * cl_HMR_Data.hpp
 *
 *  Created on: Jun 25, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_FIELD_HPP_
#define SRC_HMR_CL_HMR_FIELD_HPP_

#include <string>

#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_T_Matrix.hpp" //HMR/src

namespace moris
{
    namespace hmr
        {
//------------------------------------------------------------------------------


            class Field
            {
                //! ref to settings object
                const Parameters       * mParameters;

                //! name of field
                std::string        mLabel;

                //! pointer to Background Mesh
                Background_Mesh_Base * mBackgroundMesh;

                //! pointer to B-Spline Mesh
                BSpline_Mesh_Base    * mBSplineMesh;

                //! pointer to Lagrange Mesh
                Lagrange_Mesh_Base   * mLagrangeMesh;

                //! Matrix containing B-Spline Values
                Mat< real > mBSplineValues;

                //! Matrix containing Lagrange Values
                Mat< real > mLagrangeValues;

                //! T-Matrix object for field
                T_Matrix    mTMatrix;

                const luint mNumberOfBasis;
                const luint mNumberOfNodes ;


                //! field dimension
                uint        mDimension = 1;


//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                Field(  const Parameters     * aParameters,
                        const std::string    & aLabel,
                        Background_Mesh_Base * aBackgroundMesh,
                        BSpline_Mesh_Base    * aBSplineMesh,
                        Lagrange_Mesh_Base   * aLagrangeMesh );
//------------------------------------------------------------------------------

                /**
                 * for testing
                 */
                void
                set_bspline_values( const luint & aIndex,
                             const real  & aValue )
                {
                    mBSplineValues( aIndex ) = aValue;
                }

//------------------------------------------------------------------------------

                /**
                 * for testing
                 */
                void
                set_lagrange_values( const Mat< real >  & aValues )
                {
                    mLagrangeValues = aValues;
                }

//------------------------------------------------------------------------------

                void
                calculate_lagrange_values();

//-------------------------------------------------------------------------------

                void
                append_to_mtk_object( MTK * aMTK );

//-------------------------------------------------------------------------------

                auto
                get_order() const -> decltype( mLagrangeMesh->get_order() )
                {
                    return mLagrangeMesh->get_order();
                }

//-------------------------------------------------------------------------------

                auto
                get_label() const -> decltype( mLabel )
                {
                    return mLabel;
                }

//-------------------------------------------------------------------------------

                auto
                get_data() const -> decltype( mLagrangeValues )
                {
                    return mLagrangeValues;
                }

//-------------------------------------------------------------------------------
            };
//------------------------------------------------------------------------------
        } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FIELD_HPP_ */
