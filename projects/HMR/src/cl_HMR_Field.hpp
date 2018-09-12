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
            class HMR;
//------------------------------------------------------------------------------


            class Field
            {
                //! ref to settings object
                const Parameters   * mParameters;

                //! ref to HMR object
                HMR                * mHMR;

                //! index of lagrange mesh
                const uint           mMeshIndex;

                //! index of this field in HMR parent
                //uint           mFieldIndex = MORIS_UINT_MAX;

                //! pointer to Lagrange Mesh
                Lagrange_Mesh_Base * mMesh;

                //! pointer to T-Matrix object
                T_Matrix           * mTMatrix;

                //! Matrix containing Lagrange Values
                Mat< real >          mNodeValues;

                //! Matrix containing B-Spline coefficients
                Mat< real >          mCoefficients;

                //! name of field
                std::string          mLabel;

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

                Mat< real > &
                get_coefficients()
                {
                    return mCoefficients;
                }

//-------------------------------------------------------------------------------

                const Mat< real > &
                get_coefficients() const
                {
                    return mCoefficients;
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

//-------------------------------------------------------------------------------

                /**
                 * calculates the node values form the saved coefficients
                 */
                void
                evaluate_node_values();

//-------------------------------------------------------------------------------

                /**
                 * calculates the node values form passed coefficients
                 */
                void
                evaluate_node_values(  const Mat< real > & aCoefficients );

//-------------------------------------------------------------------------------

                /**
                 * performs an L2 projection in order to calculate coefficients
                 */
                void
                l2_project_coefficients();

//-------------------------------------------------------------------------------

                /**
                 * performs an L2 projection in order to calculate coefficients
                 * contains error and exact function for testing purposes
                 */
                void
                l2_project_coefficients(
                        real & aIntegrationError,
                        real (*aFunction)( const Mat< real > & aPoint ) );

//-------------------------------------------------------------------------------

                /**
                 * assigns memory for node values
                 */
                void
                allocate_node_values();

//-------------------------------------------------------------------------------

                /**
                 * returns the activation pattern of the Lagrange mesh
                 */
                auto
                get_activation_pattern() const
                    -> decltype( mMesh->get_activation_pattern() )
                {
                    return  mMesh->get_activation_pattern();
                }

//-------------------------------------------------------------------------------
            };
//-------------------------------------------------------------------------------
        } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FIELD_HPP_ */
