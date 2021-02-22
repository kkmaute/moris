/*
 * cl_MTK_Field.hpp
 *
 *  Created on: Jan 19, 2021
 *      Author: schmidt
 */

#ifndef PROJECTS_HMR_SRC_CL_MTK_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_MRK_FIELD_HPP_

#include <memory>

#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "st_MTK_Mesh_Pair.hpp"


namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        class Mesh;
        class Mesh_Manager;

        //------------------------------------------------------------------------------

        class Field
        {
            protected:

                mtk::Mesh_Pair * mMeshPair = nullptr;

                //! Discretization Index
                moris_index mDiscretizationMeshIndex = -1;

                std::string mLabel;

                //! index of field in mesh
                uint mFieldIndex = MORIS_UINT_MAX;

                //! //FIXME right now only scalar field
                const uint mNumberOfDimensions = 1;

                //! Nodal field values
                Matrix< DDRMat > mNodalValues;

                //! Coefficients values
                Matrix< DDRMat > mCoefficients;

                bool mFieldIsLocked = false;

                //------------------------------------------------------------------------------
            public :
                //------------------------------------------------------------------------------

                Field()
                {};

                //------------------------------------------------------------------------------

                Field(
                        mtk::Mesh_Pair * aMeshPair,
                        uint const     & aDiscretizationMeshIndex =0 )
                : mMeshPair( aMeshPair ),
                  mDiscretizationMeshIndex( aDiscretizationMeshIndex )
                {};

                Field(
                        uint        const & aDiscretizationMeshIndex,
                        std::string const & aName)
                : mDiscretizationMeshIndex( aDiscretizationMeshIndex ),
                  mLabel( aName )
                {};

                //------------------------------------------------------------------------------

                virtual ~Field();

                //------------------------------------------------------------------------------

                Mesh_Pair * get_mesh_pair();

                //------------------------------------------------------------------------------

                void set_mesh_pair( Mesh_Pair * aMeshPair);

                //------------------------------------------------------------------------------

                virtual void compute_nodal_values()
                {
                    MORIS_ERROR( false, "Field::compute_nodal_values(), Child implementation missing. ");
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the dimensionality of the field
                 */
                uint get_number_of_dimensions() const
                {
                    //FIXME right now only scalar field
                    return mNumberOfDimensions;
                }

                //------------------------------------------------------------------------------

                moris::real get_field_value( const uint & aNodeIndex )
                {
                    return mNodalValues( aNodeIndex );
                }

                //------------------------------------------------------------------------------

                void set_field_value( const uint & aFieldIndex,
                                      const real & aFieldValue )
                {
                    mNodalValues( aFieldIndex ) = aFieldValue;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation order of the Lagrange Mesh
                 */
                virtual uint get_lagrange_order() const
                {
                    return mMeshPair->mInterpolationMesh->get_order();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the discretazion order. When using HMR this is the BSpline order
                 */
                virtual uint get_discretization_order() const
                {
                    return mMeshPair->
                            mInterpolationMesh->
                            get_discretization_order( mDiscretizationMeshIndex );
                };

                //------------------------------------------------------------------------------

                const std::string & get_label() const
                {
                    return mLabel;
                };

                //------------------------------------------------------------------------------

                void set_label( const std::string & aLabel )
                {
                    mLabel = aLabel;
                };

                //------------------------------------------------------------------------------

                /**
                 * Expert user function. Allows for the unerlying mesh as well as the nodal and coefficients values to be changed.
                 * Such a change can result in an unwanted behavior.
                 */
                void unlock_field()
                {
                    mFieldIsLocked = false;
                };

                //------------------------------------------------------------------------------

                virtual const Matrix< DDRMat > & get_nodal_values() const
                {
                    return mNodalValues;
                };

                //------------------------------------------------------------------------------

                virtual void set_nodal_values( const Matrix< DDRMat > & aNodalValues )
                {
                    this->error_if_locked();

                    mNodalValues = aNodalValues;

                    mFieldIsLocked = true;
                };

                //------------------------------------------------------------------------------

                virtual const Matrix< DDRMat > & get_coefficients() const
                {
                    return mCoefficients;
                };

                //------------------------------------------------------------------------------

                virtual void set_coefficients( const Matrix< DDRMat > & aCoefficients )
                {
                    this->error_if_locked();

                    mCoefficients = aCoefficients;

                    mFieldIsLocked = true;
                };

                //------------------------------------------------------------------------------

                virtual uint get_discretization_mesh_index() const
                {
                    MORIS_ASSERT( mDiscretizationMeshIndex != -1, "get_discretization_mesh_index() Discretization index not set");
                    return mDiscretizationMeshIndex;
                }

                //------------------------------------------------------------------------------
                //
                //void get_element_local_node_values(
                //        const moris_index  aElementIndex,
                //        Matrix< DDRMat > & aValues );
                //
                //------------------------------------------------------------------------------

                /**
                 * return the field index on the linked mesh
                 */
                uint get_field_index() const
                {
                    return mFieldIndex;
                }

                //------------------------------------------------------------------------------

                void save_field_to_hdf5( const std::string & aFilePath, const bool aCreateNewFile=true );

                //------------------------------------------------------------------------------

                void save_node_values_to_hdf5( const std::string & aFilePath, const bool aCreateNewFile=true );

                //------------------------------------------------------------------------------

                void load_field_from_hdf5( const std::string & aFilePath,
                        const uint          aBSplineOrder=0 );

                //------------------------------------------------------------------------------

                void save_bspline_coeffs_to_binary( const std::string & aFilePath );

                //------------------------------------------------------------------------------

                void save_node_values_to_binary( const std::string & aFilePath );

                //------------------------------------------------------------------------------

                void save_field_to_exodus( const std::string & aFileName );

                //------------------------------------------------------------------------------

                void error_if_locked(  ) const;

        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_MRK_FIELD_HPP_ */
