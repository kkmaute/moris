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
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

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
                std::shared_ptr<mtk::Mesh_Manager> mMeshManager = nullptr;

                //! Mesh Index for Mesh Manager
                moris_index mMeshIndex = -1;

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

                //------------------------------------------------------------------------------
            public :
                //------------------------------------------------------------------------------

                Field()
                {};

                //------------------------------------------------------------------------------

                Field(
                        std::shared_ptr<mtk::Mesh_Manager>   aMeshManager,
                        uint const                         & aMeshIndex,
                        uint const                         & aDiscretizationMeshIndex =0)
                : mMeshManager( aMeshManager ),
                  mMeshIndex( aMeshIndex ),
                  mDiscretizationMeshIndex( aDiscretizationMeshIndex )
                {};

                Field(
                        uint const     & aDiscretizationMeshIndex,
                        std::string const    & aName)
                : mDiscretizationMeshIndex( aDiscretizationMeshIndex ),
                  mLabel( aName )
                {};

                //------------------------------------------------------------------------------

                virtual ~Field();

                //------------------------------------------------------------------------------

                std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > get_mesh_pair()
                {
                    MORIS_ERROR( mMeshIndex != -1, " Field::get_mesh_pair()(), Mesh pair index not set" );
                    MORIS_ERROR( mMeshManager != nullptr, " Field::get_mesh_pair()(), Mesh_Manager not set" );

                    return std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> >( mMeshIndex, mMeshManager );
                };

                //------------------------------------------------------------------------------
                /**
                 * returns the dimensionality of the field
                 */
                uint get_number_of_dimensions() const
                {
                    //FIXME right now only scalar field
                    return mNumberOfDimensions;
                }

                virtual moris::real get_field_value(  uint aNodeIndex, const Matrix< DDRMat > & aCoordinates)
                {
                    return mNodalValues( aNodeIndex );
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the interpolation order of the Lagrange Mesh
                 */
                virtual uint get_lagrange_order() const
                {
                    return mMeshManager->get_interpolation_mesh( mMeshIndex )->get_order();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the discretazion order. When using HMR this is the BSpline order
                 */
                virtual uint get_discretization_order() const
                {
                    return mMeshManager->
                            get_interpolation_mesh( mMeshIndex )->
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

                virtual Matrix< DDRMat > & get_node_values()
                {
                    return mNodalValues;
                };

                //------------------------------------------------------------------------------

                virtual const Matrix< DDRMat > & get_node_values() const
                {
                    return mNodalValues;
                };

                //------------------------------------------------------------------------------

                virtual Matrix< DDRMat > & get_coefficients()
                {
                    return mCoefficients;
                };

                //------------------------------------------------------------------------------

                virtual const Matrix< DDRMat > & get_coefficients() const
                {
                    return mCoefficients;
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
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_MRK_FIELD_HPP_ */
