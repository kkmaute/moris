/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Field.cpp
 *
 */

#include "cl_HMR_Field.hpp"

#include <iostream>
#include <cstdio>

#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_HMR_Mesh.hpp"
#include "HMR_Tools.hpp"
// HD5 c-interface
#include "hdf5.h"

#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "HDF5_Tools.hpp"

#include "cl_MTK_Mesh.hpp"
#include "fn_dot.hpp"

namespace moris::hmr
{
    //------------------------------------------------------------------------------

    Field::Field(
            const std::string&           aLabel,
            std::shared_ptr< mtk::Mesh > aMesh,
            uint                         aBSplineMeshIndex,
            std::shared_ptr< Database >  aDatabase,
            Lagrange_Mesh_Base *         aLagrangeMesh )
    : mMesh( aMesh ),
      mDatabase( aDatabase ),
      mLagrangeMesh( aLagrangeMesh ),
      mFieldIndex( aLagrangeMesh->create_real_scalar_field_data( aLabel ) )
    {
        this->set_label( aLabel );
        mInputBSplineIndex = aBSplineMeshIndex;         //FIXME Index not order

        // assume input and output order are the same
        mOutputBSplineOrder = mInputBSplineIndex;

        //get Lagrange mesh than bspline order
        uint tBSplineOrder = mLagrangeMesh->get_bspline_order( aBSplineMeshIndex );

        aLagrangeMesh->set_real_scalar_field_bspline_order( mFieldIndex, tBSplineOrder );
    }

    //------------------------------------------------------------------------------

    Field::Field(
            const std::string&      aLabel,
            std::shared_ptr< Mesh > aMesh,
            const std::string&      aHdf5FilePath  )
    : mMesh( aMesh )
    {
        // link to database
        mDatabase = aMesh->get_database();

        mLagrangeMesh = aMesh->get_lagrange_mesh();

        // create data on target mesh
        mFieldIndex = mLagrangeMesh->create_real_scalar_field_data( aLabel );

        // set my label
        this->set_label( aLabel );

        // load values into field
        herr_t tStatus = 0;
        hid_t tHDF5File = open_hdf5_file( aHdf5FilePath );
        load_matrix_from_hdf5_file( tHDF5File, aLabel, this->get_coefficients(), tStatus );
        close_hdf5_file( tHDF5File );

        // get number of coeffs
        uint tNumberOfBSplines = this->get_coefficients().length();

        // find out order
        mInputBSplineIndex = 0;
        for( uint k=1; k<=3; ++k )
        {
            if( aMesh->get_max_num_coeffs_on_proc( k ) == tNumberOfBSplines )
            {
                mInputBSplineIndex = k;
                break;
            }
        }

        // assume input and output order are the same
        mOutputBSplineOrder = mInputBSplineIndex;

        // make sure that we have found an order
        MORIS_ERROR( mInputBSplineIndex != 0, "Could not find corresponding B-Spline mesh for passed coefficients" );

        // set order of this field
        mLagrangeMesh->set_real_scalar_field_bspline_order( mFieldIndex, mInputBSplineIndex );

        this->evaluate_nodal_values();
    }

    //------------------------------------------------------------------------------

    Field::~Field()
    {
        mLagrangeMesh = nullptr;
    }

    //------------------------------------------------------------------------------

    // parameter copied from input settings
    void Field::set_min_surface_level( uint aLevel )
    {
        mMinSurfaceLevel = aLevel;
    }

    //------------------------------------------------------------------------------

    void Field::set_min_volume_level( uint aLevel )
    {
        mMinVolumeLevel = aLevel;
    }

    //------------------------------------------------------------------------------

    void Field::set_max_surface_level( uint aLevel )
    {
        mMaxSurfaceLevel = aLevel;
    }

    //------------------------------------------------------------------------------

    void Field::set_max_volume_level( uint aLevel )
    {
        mMaxVolumeLevel = aLevel;
    }

    //------------------------------------------------------------------------------

    uint Field::get_min_surface_level() const
    {
        return mMinSurfaceLevel;
    }

    //------------------------------------------------------------------------------

    // parameter copied from input settings
    uint Field::get_min_volume_level() const
    {
        return mMinVolumeLevel;
    }

    //------------------------------------------------------------------------------

    // parameter copied from input settings
    uint Field::get_max_surface_level() const
    {
        return mMaxSurfaceLevel;
    }

    //------------------------------------------------------------------------------

    // parameter copied from input settings
    uint Field::get_max_volume_level() const
    {
        return mMaxVolumeLevel;
    }

    //------------------------------------------------------------------------------

    const std::string & Field::get_label() const
    {
        return mLagrangeMesh->get_real_scalar_field_label( mFieldIndex );
    }

    //------------------------------------------------------------------------------

    void Field::set_label( const std::string & aLabel )
    {
        mLagrangeMesh->set_real_scalar_field_label( mFieldIndex, aLabel );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > & Field::get_node_values()
    {
        return mLagrangeMesh->get_real_scalar_field_data( mFieldIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > & Field::get_node_values() const
    {
        return mLagrangeMesh->get_real_scalar_field_data( mFieldIndex );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > & Field::get_coefficients()
    {
        return mLagrangeMesh->get_real_scalar_field_coeffs( mFieldIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat > & Field::get_coefficients() const
    {
        return mLagrangeMesh->get_real_scalar_field_coeffs( mFieldIndex );
    }

    //------------------------------------------------------------------------------

    // change field to the mesh which is deliverd as an input
    void Field::change_mesh(
            Lagrange_Mesh_Base* aMesh,
            const uint          aFieldIndex )
    {
        // set order of B-Spline
        mLagrangeMesh->set_real_scalar_field_bspline_order(
                mFieldIndex,
                this->get_bspline_order() );

        // change mesh pointer
        mLagrangeMesh = aMesh;

        // change mesh index
        mFieldIndex   = aFieldIndex;
    }
    //------------------------------------------------------------------------------

    void Field::get_element_local_node_values(
            const moris_index    aElementIndex,
            Matrix< DDRMat >   & aValues )
    {
        // get pointer to element
        Element * tElement = mLagrangeMesh->get_element( aElementIndex );

        // get number of nodes
        uint tNumberOfNodes = tElement->get_number_of_vertices();

        // allocate output matrix
        aValues.set_size( tNumberOfNodes, 1 );

        // write values into matrix
        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            aValues( k ) = this->get_node_values()( tElement->get_basis_function( k )->get_index() );
        }
    }

    //------------------------------------------------------------------------------

    void Field::evaluate_nodal_values( const Matrix< DDRMat > & aCoefficients )
    {
        // ask mesh for number of nodes
        uint tNumberOfNodes = mMesh->get_num_nodes();

        Matrix< DDRMat > & tNodeValues = this->get_node_values();

        // allocate memory for matrix
        tNodeValues.set_size( tNumberOfNodes, mNumberOfDimensions );

        MORIS_ERROR( mNumberOfDimensions == 1, "currently, only scalar fields are supported" );

        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            // get pointer to node
            auto tNode = &mMesh->get_mtk_vertex( k );

            // get PDOFs from node
            auto tBSplines = tNode->get_interpolation( mInputBSplineIndex )->get_coefficients();

            // get T-Matrix
            const Matrix< DDRMat > & tTMatrix = *tNode->get_interpolation( mInputBSplineIndex )->get_weights();

            // get number of coefficients
            uint tNumberOfCoeffs = tTMatrix.length();

            MORIS_ASSERT( tNumberOfCoeffs > 0, "No coefficients defined for node" ) ;

            // fill coeffs vector
            Matrix< DDRMat > tCoeffs( tNumberOfCoeffs, 1 );
            for( uint i = 0; i < tNumberOfCoeffs; ++i )
            {
                tCoeffs( i ) = aCoefficients( tBSplines( i )->get_index() );
            }

            // write value into solution
            tNodeValues( k ) = dot( tTMatrix, tCoeffs );
        }
    }

    //------------------------------------------------------------------------------

    void Field::evaluate_nodal_values()
    {
        this->evaluate_nodal_values( this->get_coefficients() );
    }

    //------------------------------------------------------------------------------

    void Field::put_scalar_values_on_field( const Matrix< DDRMat > & aValues )
    {
        // get pointer to node values
        Matrix< DDRMat > & tNodeValues = this->get_node_values();

        // get number of nodes on block
        uint tNumberOfVertices = mMesh->get_num_nodes();

        // set size of node values
        tNodeValues.set_size( tNumberOfVertices, 1 );
        tNodeValues.fill(0);

        // loop over all vertices
        for( uint k=0; k<tNumberOfVertices; ++k )
        {
            tNodeValues( k ) = aValues( k );
        }
    }

    //------------------------------------------------------------------------------

    void Field::save_field_to_hdf5(
            const std::string & aFilePath,
            const bool          aCreateNewFile )
    {
        // test if file exists
        std::string tFilePath = make_path_parallel( aFilePath );

        // test if file exists
        std::ifstream tFile( tFilePath );
        bool tFileExists;
        if( tFile )
        {
            tFileExists = true;
        }
        else
        {
            tFileExists = false;
        }

        tFile.close();

        // delete file if it exists and user does not want to keep it
        if( aCreateNewFile && tFileExists )
        {
            std::remove( tFilePath.c_str() );
            tFileExists = false;
        }

        hid_t tFileID;

        if( tFileExists )
        {
            tFileID = open_hdf5_file( aFilePath );
        }
        else
        {
            tFileID = create_hdf5_file( aFilePath );
        }

        herr_t tStatus;

        save_matrix_to_hdf5_file( tFileID,
                this->get_label(),
                this->get_coefficients(),
                tStatus );

        // close file
        tStatus = close_hdf5_file( tFileID );
    }

    //------------------------------------------------------------------------------

    void Field::save_node_values_to_hdf5(
            const std::string & aFilePath,
            const bool          aCreateNewFile )
    {
        // test if file exists
        std::string tFilePath = make_path_parallel( aFilePath );

        // test if file exists
        std::ifstream tFile( tFilePath );
        bool tFileExists;
        if( tFile )
        {
            tFileExists = true;
        }
        else
        {
            tFileExists = false;
        }

        tFile.close();

        // delete file if it exists and user does not want to keep it
        if( aCreateNewFile && tFileExists )
        {
            std::remove( tFilePath.c_str() );
            tFileExists = false;
        }

        hid_t tFileID;

        if( tFileExists )
        {
            tFileID = open_hdf5_file( aFilePath );
        }
        else
        {
            tFileID = create_hdf5_file( aFilePath );
        }

        herr_t tStatus;

        save_matrix_to_hdf5_file( tFileID,
                this->get_label(),
                this->get_node_values(),
                tStatus );

        // close file
        tStatus = close_hdf5_file( tFileID );
    }

    //------------------------------------------------------------------------------

    void Field::load_field_from_hdf5(
            const std::string & aFilePath,
            const uint          aBSplineOrder )
    {
        hid_t tFile    = open_hdf5_file( aFilePath );
        herr_t tStatus = 0;
        load_matrix_from_hdf5_file( tFile,
                this->get_label(),
                this->get_coefficients(),
                tStatus );

        tStatus = close_hdf5_file( tFile );
    }

    //------------------------------------------------------------------------------

    void Field::save_node_values_to_binary( const std::string & aFilePath )
    {
        // make path parallel
        std::string tFilePath = parallelize_path( aFilePath );

        save_matrix_to_binary_file( this->get_node_values(), tFilePath );
    }

    //------------------------------------------------------------------------------

    void Field::save_bspline_coeffs_to_binary( const std::string & aFilePath )
    {
        // make path parallel
        std::string tFilePath = parallelize_path( aFilePath );

        save_matrix_to_binary_file( this->get_coefficients(), tFilePath );
    }

    //------------------------------------------------------------------------------

    void Field::set_bspline_order( uint aOrder )
    {
        mLagrangeMesh->set_real_scalar_field_bspline_order( mFieldIndex, aOrder );
    }

    //------------------------------------------------------------------------------

    mtk::EntityRank Field::get_bspline_rank() const
    {
        return mtk::EntityRank::BSPLINE;
    }

    //------------------------------------------------------------------------------
} /* namespace moris */
