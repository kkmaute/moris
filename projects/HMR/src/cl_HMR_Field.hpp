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
//#include "cl_HMR_T_Matrix.hpp" //HMR/src

namespace moris
{
    namespace hmr
        {

//------------------------------------------------------------------------------
            class HMR;
//------------------------------------------------------------------------------


            class Field
            {
                //! ref to settings object
                const Parameters   * mParameters;

                //! index of lagrange mesh
                const uint           mMeshIndex;

                //! index of this field in HMR parent
                const uint           mFieldIndex;

                //! pointer to Lagrange Mesh
                Lagrange_Mesh_Base * mMesh;

                //! Matrix containing Lagrange Values
                Mat< real >        & mNodeValues;

                //! name of field
                std::string        & mLabel;

                //! Matrix containing B-Spline Values
                //Mat< real > mBSplineValues;

                //! field dimension
                uint                 mDimension = 1;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                /**
                 * constructor
                 */
                Field(  HMR                * aHMR,
                        const std::string  & aLabel,
                        const uint         & aLagrangeMeshIndex );

//------------------------------------------------------------------------------

               /**
                * returns the index of the linked Lagrange mesh
                */
                auto
                get_lagrange_index() const -> decltype ( mMeshIndex )
                {
                    return mMeshIndex;
                }

//-------------------------------------------------------------------------------

                /**
                 * returns the interpolation order of the underlying mesh
                 */
                auto
                get_order() const -> decltype( mMesh->get_order() )
                {
                    return mMesh->get_order();
                }

//-------------------------------------------------------------------------------

                std::string
                get_label() const
                {
                    return mLabel;
                }

//-------------------------------------------------------------------------------

                Mat< real > &
                get_data()
                {
                    return mNodeValues;
                }

//-------------------------------------------------------------------------------

                const Mat< real > &
                get_data() const
                {
                    return mNodeValues;
                }

//-------------------------------------------------------------------------------

                /**
                 * returns the number of nodes of this field
                 */
                auto
                get_number_of_nodes() const
                    -> decltype ( mMesh->get_number_of_nodes_on_proc() )
                {
                    return mMesh->get_number_of_nodes_on_proc();
                }

//-------------------------------------------------------------------------------

                /**
                 * evaluates a function on the mesh and writes data into field
                 */
                void
                evaluate_function(  real (*aFunction)( const Mat< real > & aPoint ) );

//-------------------------------------------------------------------------------

                /**
                 * expose mesh pointer
                 */
                auto
                get_mesh() -> decltype ( mMesh )
                {
                    return mMesh;
                }

//-------------------------------------------------------------------------------

                auto
                get_number_of_dimensions() const -> decltype( mDimension )
                {
                    return mDimension;
                }

            };
//-------------------------------------------------------------------------------
        } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FIELD_HPP_ */
