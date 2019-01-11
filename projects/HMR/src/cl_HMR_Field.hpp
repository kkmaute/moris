/*
 * cl_HMR_Field.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_

#include <memory>

#include "../../../HMR/src/cl_HMR_Database.hpp"
#include "../../../HMR/src/cl_HMR_Lagrange_Mesh_Base.hpp"
#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------
        class Mesh;
//------------------------------------------------------------------------------

        class Field
        {
            //! pointer to mesh or block object this field refers to
            const std::shared_ptr< mtk::Mesh > mMesh;

            //! pointer to database
            std::shared_ptr< Database > mDatabase;

            //! mesh that holds data
            Lagrange_Mesh_Base * mLagrangeMesh;

            uint mInputBSplineOrder = 0;
            uint mOutputBSplineOrder = 0;

            // index of field in mesh
            uint mFieldIndex;

            // parameter for minumum volume refinement
            uint mMinVolumeLevel = 0;

            // parameter for minimun surface refinement
            uint mMinSurfaceLevel = 0;

            // parameter for maximum volume refinement
            uint mMaxVolumeLevel = gMaxNumberOfLevels;

            // parameter for maximum surface refinement
            uint mMaxSurfaceLevel = gMaxNumberOfLevels;

            //! Dimensionality of the field, currently fixed to 1
            const uint     mNumberOfDimensions = 1;

            //! id of field, if set
            moris_id mID = gNoID;

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            Field(  const std::string             & aLabel,
                    std::shared_ptr< mtk::Mesh >    aMesh,
                    const uint                    & aBSplineOrder,
                    std::shared_ptr< Database >     aDatabase,
                    Lagrange_Mesh_Base *            aLagrangeMesh );

//------------------------------------------------------------------------------

            Field( const std::string             & aLabel,
                   std::shared_ptr< Mesh >         aMesh,
                   const std::string             & aHdf5FilePath );

//------------------------------------------------------------------------------

            ~Field();

//------------------------------------------------------------------------------

            /**
             * returns the dimensionality of the field
             */
            uint
            get_number_of_dimensions() const
            {
                return mNumberOfDimensions;
            }

//------------------------------------------------------------------------------
            /**
             * returns the interpolation order of the Lagrange Mesh
             */
            uint
            get_lagrange_order() const
            {
                return mLagrangeMesh->get_order();
            }

//------------------------------------------------------------------------------
            /**
             * returns the interpolation order of the B-Splines
             */
            uint
            get_bspline_order() const
            {
                return mLagrangeMesh->get_real_scalar_field_bspline_order( mFieldIndex );
            }

//------------------------------------------------------------------------------
            uint
            get_bspline_output_order() const
            {
                return mOutputBSplineOrder;
            }

//------------------------------------------------------------------------------

            void
            set_bspline_output_order( const uint & aOrder )
            {
                mOutputBSplineOrder = aOrder;
            }

//------------------------------------------------------------------------------

            void
            set_id( const moris_id & aID )
            {
                mID = aID;
            }

//------------------------------------------------------------------------------

            moris_id
            get_id() const
            {
                return mID;
            }

//------------------------------------------------------------------------------

            // parameter copied from input settings
            void
            set_min_surface_level( const uint & aLevel );

//------------------------------------------------------------------------------

            // parameter copied from input settings
            void
            set_min_volume_level( const uint & aLevel );

//------------------------------------------------------------------------------

            // parameter copied from input settings
            void
            set_max_surface_level( const uint & aLevel );

//------------------------------------------------------------------------------

            // parameter copied from input settings
            void
            set_max_volume_level( const uint & aLevel );

//------------------------------------------------------------------------------

            uint
            get_min_surface_level() const;

//------------------------------------------------------------------------------

            // parameter copied from input settings
            uint
            get_min_volume_level() const;

//------------------------------------------------------------------------------

            // parameter copied from input settings
            uint
            get_max_surface_level() const;

//------------------------------------------------------------------------------

            // parameter copied from input settings
            uint
            get_max_volume_level() const;

//------------------------------------------------------------------------------

            const std::string &
            get_label() const;

//------------------------------------------------------------------------------

            void
            set_label( const std::string & aLabel );

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_node_values();

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_node_values() const;

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_coefficients();
//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_coefficients() const;

//------------------------------------------------------------------------------

            /**
             * sets the pointer of the mesh to another mesh
             * this is needed by the mapper
             */
            void
            change_mesh( Lagrange_Mesh_Base * aMesh, const uint aFieldIndex );

//------------------------------------------------------------------------------

            /**
             * returns the pointer of the underlying mesh
             */

            Lagrange_Mesh_Base *
            get_mesh()
            {
                return mLagrangeMesh;
            }

            const Lagrange_Mesh_Base *
            get_mesh() const
            {
                return mLagrangeMesh;
            }

//------------------------------------------------------------------------------

            void
            get_element_local_node_values(
                    const moris_index  aElementIndex,
                    Matrix< DDRMat > & aValues );

//------------------------------------------------------------------------------

            /**
             * return the field index on the linked mesh
             */
            uint
            get_field_index() const
            {
                return mFieldIndex;
            }

//------------------------------------------------------------------------------

            mtk::Interpolation_Order
            get_interpolation_order() const
            {
                // assume that all elements on mesh have same order
                return mMesh->get_mtk_cell( 0 ).get_interpolation_order();
            }

//------------------------------------------------------------------------------

            /**
             * returns the rank of the B-Spline interpolation
             */
            EntityRank
            get_bspline_rank() const;

//------------------------------------------------------------------------------

            void
            evaluate_scalar_function(
                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

            void
            evaluate_node_values();

//------------------------------------------------------------------------------

            void
            evaluate_node_values( const Matrix< DDRMat > & aCoefficients );

//------------------------------------------------------------------------------

            void
            save_field_to_hdf5( const std::string & aFilePath, const bool aCreateNewFile=true );

//------------------------------------------------------------------------------

            void
            load_field_from_hdf5(
                    const std::string & aFilePath,
                    const uint          aBSplineOrder=0 );

//------------------------------------------------------------------------------

            void
            save_bspline_coeffs_to_binary( const std::string & aFilePath );

//------------------------------------------------------------------------------

            void
            save_node_values_to_binary( const std::string & aFilePath );

//------------------------------------------------------------------------------

            void
            set_bspline_order( const uint & aOrder );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_ */
