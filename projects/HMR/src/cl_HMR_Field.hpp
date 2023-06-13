/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Field.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_

#include <memory>

#include "cl_HMR_Database.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris::hmr
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

            uint mInputBSplineIndex = 0;
            uint mOutputBSplineOrder = 0;

            // index of field in mesh
            uint mFieldIndex = MORIS_UINT_MAX;

            // parameter for minimum volume refinement
            uint mMinVolumeLevel = 0;

            // parameter for minimum surface refinement
            uint mMinSurfaceLevel = 0;

            // parameter for maximum volume refinement
            uint mMaxVolumeLevel = gMaxNumberOfLevels;

            // parameter for maximum surface refinement
            uint mMaxSurfaceLevel = gMaxNumberOfLevels;

            //! Dimensionality of the field, currently fixed to 1
            const uint mNumberOfDimensions = 1;

            //! id of field, if set
            moris_id mID = gNoID;

            //------------------------------------------------------------------------------
        public :
            //------------------------------------------------------------------------------

            Field(
                    const std::string            & aLabel,
                    std::shared_ptr< mtk::Mesh >   aMesh,
                    uint aBSplineIndex,
                    std::shared_ptr< Database >    aDatabase,
                    Lagrange_Mesh_Base           * aLagrangeMesh );

            //------------------------------------------------------------------------------

            Field( const std::string        & aLabel,
                    std::shared_ptr< Mesh >   aMesh,
                    const std::string       & aHdf5FilePath );

            //------------------------------------------------------------------------------

            ~Field();

            //------------------------------------------------------------------------------

            /**
             * returns the dimensionality of the field
             */
            uint get_number_of_dimensions() const
            {
                return mNumberOfDimensions;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation order of the Lagrange Mesh
             */
            uint get_lagrange_order() const
            {
                return mLagrangeMesh->get_order();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the interpolation order of the B-Splines
             */
            uint get_bspline_order() const
            {
                return mLagrangeMesh->get_real_scalar_field_bspline_order( mFieldIndex );
            }

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation order of the Lagrange Mesh
             */
            uint get_lagrange_pattern() const
            {
                return mLagrangeMesh->get_activation_pattern();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the mesh index of the Lagrange Mesh
             */
            uint get_lagrange_mesh_index() const
            {
                return mLagrangeMesh->get_index();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the interpolation order of the B-Splines
             */
            uint get_bspline_pattern() const
            {
                MORIS_ASSERT( false, "does not retun pattern, returns order");
                return mLagrangeMesh->get_real_scalar_field_bspline_order( mFieldIndex );
            }

            //------------------------------------------------------------------------------

            uint get_bspline_output_order() const
            {
                return mOutputBSplineOrder;
            }

            //------------------------------------------------------------------------------

            void set_bspline_output_order( uint aOrder )
            {
                mOutputBSplineOrder = aOrder;
            }

            //------------------------------------------------------------------------------

            void set_id( const moris_id & aID )
            {
                mID = aID;
            }

            //------------------------------------------------------------------------------

            moris_id get_id() const
            {
                MORIS_ASSERT( mID != gNoID, "Field ID not initialized");

                return mID;
            }

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            void set_min_surface_level( uint aLevel );

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            void set_min_volume_level( uint aLevel );

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            void set_max_surface_level( uint aLevel );

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            void set_max_volume_level( uint aLevel );

            //------------------------------------------------------------------------------

            uint get_min_surface_level() const;

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            uint get_min_volume_level() const;

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            uint get_max_surface_level() const;

            //------------------------------------------------------------------------------

            // parameter copied from input settings
            uint get_max_volume_level() const;

            //------------------------------------------------------------------------------

            const std::string & get_label() const;

            //------------------------------------------------------------------------------

            void set_label( const std::string & aLabel );

            //------------------------------------------------------------------------------

            Matrix< DDRMat > & get_node_values();

            //------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_node_values() const;

            //------------------------------------------------------------------------------

            Matrix< DDRMat > & get_coefficients();
            //------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_coefficients() const;

            //------------------------------------------------------------------------------

            /**
             * sets the pointer of the mesh to another mesh
             * this is needed by the mapper
             */
            void change_mesh(
                    Lagrange_Mesh_Base * aMesh,
                    const uint           aFieldIndex );

            //------------------------------------------------------------------------------

            /**
             * returns the pointer of the underlying mesh
             */

            Lagrange_Mesh_Base * get_mesh()
            {
                return mLagrangeMesh;
            }

            const Lagrange_Mesh_Base * get_mesh() const
            {
                return mLagrangeMesh;
            }

            //------------------------------------------------------------------------------

            void get_element_local_node_values(
                    const moris_index  aElementIndex,
                    Matrix< DDRMat > & aValues );

            //------------------------------------------------------------------------------

            /**
             * return the field index on the linked mesh
             */
            uint get_field_index() const
            {
                return mFieldIndex;
            }

            //------------------------------------------------------------------------------

            mtk::Interpolation_Order get_interpolation_order() const
            {
                // assume that all elements on mesh have same order
                return mMesh->get_mtk_cell( 0 ).get_interpolation_order();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the rank of the B-Spline interpolation
             */
            EntityRank get_bspline_rank() const;

            //------------------------------------------------------------------------------

            template<typename T>
            void evaluate_scalar_function( T aLambda )
            {
                // get pointer to node values
                Matrix< DDRMat > & tNodeValues = this->get_node_values();

                // get number of nodes on block
                uint tNumberOfVertices = mMesh->get_num_nodes();

                // set size of node values
                tNodeValues.set_size( tNumberOfVertices, 1 );

                // loop over all vertices
                for( uint k=0; k<tNumberOfVertices; ++k )
                {
                    // evaluate function at vertex coordinates
                    tNodeValues( k ) = aLambda( mMesh->get_mtk_vertex( k ).get_coords() );
                }
            }

            //------------------------------------------------------------------------------

            void put_scalar_values_on_field( const Matrix< DDRMat > & aValues );

            //------------------------------------------------------------------------------

            void evaluate_nodal_values();

            //------------------------------------------------------------------------------

            void evaluate_nodal_values( const Matrix< DDRMat > & aCoefficients );

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

            void set_bspline_order( uint aOrder );

            //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_ */

