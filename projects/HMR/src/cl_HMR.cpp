

/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

// see http://pubs.opengroup.org/onlinepubs/7908799/xsh/dlfcn.h.html
#include "dlfcn.h"

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_times.hpp" //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_eye.hpp" //LINALG/src
#include "fn_unique.hpp" //LINALG/src
#include "fn_print.hpp" //LINALG/src
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//#include "cl_DLA_Solver_Factory.hpp"
//#include "cl_DLA_Linear_Solver_Aztec.hpp"
//#include "cl_Vector.hpp"
//
//#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
//#include "cl_NLA_Newton_Solver.hpp"
//#include "cl_NLA_Nonlinear_Problem.hpp"
//#include "cl_NLA_Solver_Interface_Proxy.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "HDF5_Tools.hpp"
#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_STK.hpp" //HMR/src
#include "cl_HMR_File.hpp" //HMR/src

#include "cl_HMR_Field.hpp"          //HMR/src
#include "cl_HMR_Background_Element_Base.hpp"

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        HMR::HMR ( Parameters * aParameters ) :
                mParameters( aParameters )
        {
            mDatabase   = std::make_shared< Database >( aParameters );

            this->create_input_and_output_meshes();

            mDatabase ->calculate_t_matrices_for_input();
        }

// -----------------------------------------------------------------------------

        // alternative constuctor that converts ref to a pointer
        HMR::HMR ( Parameters & aParameters ) :
                                HMR( & aParameters )
        {

        }

// -----------------------------------------------------------------------------

        // alternative constuctor that uses parameter list
        HMR::HMR ( ParameterList & aParameterList )
            : HMR( new Parameters( aParameterList ) )
        {
            mDatabase->set_parameter_owning_flag();
        }

// -----------------------------------------------------------------------------

        HMR::HMR( const std::string & aPath )
        {
            mDatabase = std::make_shared< Database >( aPath );
            // set shared pointer of database to itself

            mDatabase->set_parameter_owning_flag();

            // set parameters of HMR object
            mParameters = mDatabase->get_parameters();

            this->create_input_and_output_meshes();

            mDatabase->calculate_t_matrices_for_input();
        }

// -----------------------------------------------------------------------------

        HMR::HMR( const std::string & aInPath, const std::string & aOutPath )
        {
            mDatabase = std::make_shared< Database >( aInPath, aOutPath );

            // set shared pointer of database to itself
            mDatabase->set_parameter_owning_flag();

            // set parameters of HMR object
            mParameters = mDatabase->get_parameters();

            //mDatabase->calculate_t_matrices_for_input();

            // create union of input and output
            mDatabase->create_union_pattern();

            mDatabase->set_activation_pattern( mParameters->get_bspline_output_pattern() );

            // update database
            mDatabase->update_bspline_meshes();
            mDatabase->update_lagrange_meshes();

            // finalize database
            this->finalize();

            this->create_input_and_output_meshes();
        }

// -----------------------------------------------------------------------------

        void
        HMR::finalize()

        {
            // if mesh has not been refined, copy input to output before finalizing
            if( ! mDatabase->have_refined_at_least_one_element() )
            {
                // copy input to output
                mDatabase->get_background_mesh()->copy_pattern(
                        mParameters->get_bspline_input_pattern(),
                        mParameters->get_bspline_output_pattern() );

                // copy input to output
                mDatabase->get_background_mesh()->copy_pattern(
                        mParameters->get_lagrange_input_pattern(),
                        mParameters->get_lagrange_output_pattern() );

                mDatabase->get_background_mesh()->copy_pattern(
                        mParameters->get_bspline_input_pattern(),
                        mParameters->get_union_pattern() );

                // select output pattern
                mDatabase->set_activation_pattern(
                        mParameters->get_lagrange_output_pattern()  );

                // update database
                mDatabase->update_bspline_meshes();
                mDatabase->update_lagrange_meshes();
            }

            // finish database
            mDatabase->finalize();

        }
// -----------------------------------------------------------------------------
        void
        HMR::save_to_exodus( const std::string & aPath, const double aTimeStep,  const uint aOutputOrder )
        {

            uint tOutputOrder = MORIS_UINT_MAX;
            uint tIndex = MORIS_UINT_MAX;

            if( aOutputOrder == 0 )
            {
                tOutputOrder = mParameters->get_lagrange_orders().max();
            }
            else
            {
                tOutputOrder = aOutputOrder;
            }

            if( tOutputOrder > 2 )
            {
                // we cant output a cubic mesh
                tIndex = mParameters->get_refined_output_mesh();
            }
            else
            {
                tIndex = this->get_mesh_index( tOutputOrder, mParameters->get_lagrange_output_pattern() );
            }


            MORIS_ERROR( tIndex != MORIS_UINT_MAX,
                    "Something went wrong while trying to find mesh for exodus file" );

            MORIS_ASSERT( mDatabase->get_lagrange_mesh_by_index( tIndex )->get_order() == tOutputOrder,
                    "Picked wrong mesh for output");

            this->save_to_exodus(
                    tIndex,
                    aPath,
                    aTimeStep );
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_last_step_to_exodus(
                            const std::string & aPath,
                            const double aTimeStep,
                            const uint aOutputOrder )
        {
            MORIS_ERROR( ! mUpdateRefinementCalled,
                    "HMR does not feel comfortable with you calling save_last_step_to_exodus() after you have overwritten the input pattern using update_refinement_pattern()");

            uint tOutputOrder = MORIS_UINT_MAX;
            uint tIndex = MORIS_UINT_MAX;

            if( aOutputOrder == 0 )
            {
                tOutputOrder = mParameters->get_lagrange_orders().max();
            }
            else
            {
                tOutputOrder = aOutputOrder;
            }

            if( tOutputOrder > 2 )
            {
                // we cant output a cubic mesh
                tIndex = mParameters->get_refined_output_mesh();
            }
            else
            {
                tIndex = this->get_mesh_index( tOutputOrder, mParameters->get_lagrange_input_pattern() );
            }

            MORIS_ERROR( tIndex != MORIS_UINT_MAX,
                    "Something went wrong while trying to find mesh for exodus file" );

            this->save_to_exodus( tIndex,
                                   aPath,
                                   aTimeStep );
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_to_exodus(
                const uint        & aMeshIndex,
                const std::string & aPath,
                const double aTimeStep  )
        {
            STK * tSTK = mDatabase->get_lagrange_mesh_by_index( aMeshIndex )
                    ->create_stk_object( aTimeStep );

            // save MTK to exodus
            tSTK->save_to_file( aPath );

            // delete file
            delete tSTK;

        }

// -----------------------------------------------------------------------------

        void
        HMR::save_to_hdf5( const std::string & aPath )
        {
            // create file object
            File tHDF5;

            // create file on disk
            tHDF5.create( aPath );

            // store settings object
            tHDF5.save_settings( mParameters );

            // get pointer to background mesh
            Background_Mesh_Base * tBackgroundMesh = mDatabase->get_background_mesh();

            // remember active pattern
            auto tActivePattern = tBackgroundMesh->get_activation_pattern();


            // save output pattern into file
            tHDF5.save_refinement_pattern( tBackgroundMesh );

            if( tActivePattern != tBackgroundMesh->get_activation_pattern() )
            {
                tBackgroundMesh->set_activation_pattern( tActivePattern );
            }

            // close hdf5 file
            tHDF5.close();
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_coeffs_to_hdf5_file( const std::string & aFilePath )
        {
            // fixme: at the moment, this function only stores the coeffs
            // of one order. This function will be modified soon so that
            // multiple refinement orders are supported

            // MORIS_ERROR( false, "save_coeffs_to_hdf5_file() is currently out of order " );


            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            // get pointer to output mesh
            Lagrange_Mesh_Base * tMesh = nullptr;

            for( uint k=0; k<mDatabase->get_number_of_lagrange_meshes(); ++k )
            {
                tMesh = mDatabase->get_lagrange_mesh_by_index( k );
                if( tMesh->get_activation_pattern() == mParameters->get_lagrange_output_pattern() )
                {
                    // cancel the loop. We only save one mesh and one order
                    break;
                }
            }

            // Create a new file using default properties
           herr_t tFileID = H5Fcreate(
                    tFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // error handler
            herr_t tStatus;

            // get number of nodes of this mesh
            uint tNumberOfNodes = tMesh->get_number_of_nodes_on_proc();

            // allocate matrix with ids
            Matrix< IdMat > tIDs( tNumberOfNodes, 1 );

            // populate matrix
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                tIDs( k ) = tMesh->get_node_by_index( k )->get_id();
            }

            // save ids to file
            save_matrix_to_hdf5_file(
                    tFileID,
                    "NodeID",
                    tIDs,
                    tStatus );

            // loop over all B-Spline meshes
            uint tNumberOfBSplineMeshes = tMesh->get_number_of_bspline_meshes();

            for ( uint m=0; m<tNumberOfBSplineMeshes; ++m )
            {
                // get pointer to mesh
                BSpline_Mesh_Base * tBMesh = tMesh->get_bspline_mesh( m );

                if ( tBMesh != NULL )
                {
                    // get order of mesh
                    uint tOrder = tBMesh->get_order();


                    // generate label
                    std::string tLabel = "NumberOfCoefficients_" + std::to_string( tOrder );

                    // count number of coefficients per node
                    Matrix< DDUMat > tNumberOfCoeffs( tNumberOfNodes, 1, 0 );

                    // populate matrix
                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        tNumberOfCoeffs( k ) = tMesh->get_node_by_index( k )
                                                    ->get_interpolation( tOrder )
                                                    ->get_number_of_coefficients();
                    }

                    // save number of coeffs to file
                    save_matrix_to_hdf5_file(
                            tFileID,
                            tLabel,
                            tNumberOfCoeffs,
                            tStatus );

                    // get max number of coeffs
                    uint tMaxNumCoeffs = tNumberOfCoeffs.max();

                    Matrix< IdMat > tCoeffIDs( tNumberOfNodes, tMaxNumCoeffs, gNoID );
                    Matrix< DDRMat >  tWeights( tNumberOfNodes, tMaxNumCoeffs, 0.0 );

                    // populate matrix
                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        // get max number of dofs
                        uint tMaxI = tNumberOfCoeffs( k );

                        // get pointer to interpolation object
                        mtk::Vertex_Interpolation * tInterp = tMesh
                                ->get_node_by_index( k )
                                ->get_interpolation( tOrder );


                        Matrix< IdMat >    tLocalIDs = tInterp->get_ids();
                        const Matrix< DDRMat > & tLocalWeights = *tInterp->get_weights();

                        // copy data into global matrix
                        for( uint i=0; i<tMaxI; ++i )
                        {
                            tCoeffIDs( k, i ) = tLocalIDs( i );
                            tWeights( k, i ) = tLocalWeights( i );
                        }
                    }
                    // generate label
                    tLabel = "BSplineIDs_" + std::to_string( tOrder );

                    // save ids to file
                    save_matrix_to_hdf5_file(
                            tFileID,
                            tLabel,
                            tCoeffIDs,
                            tStatus );

                    // generate  label
                    tLabel = "InterpolationWeights_" + std::to_string( tOrder );

                    // save weights to file
                    save_matrix_to_hdf5_file(
                            tFileID,
                            tLabel,
                            tWeights,
                            tStatus );

                }
            }

            // close file
            tStatus = H5Fclose( tFileID );
        }

// -----------------------------------------------------------------------------

        void
        HMR::flag_elements(
                      Cell< mtk::Cell* > & aElements,
                const uint                  aMinRefinementLevel )
        {

            // get  working pattern
            uint tWorkingPattern = mParameters->get_working_pattern();

            // get pointer to background mesh
            Background_Mesh_Base * tBackgroundMesh = mDatabase->get_background_mesh();


            // use Lagrange pattern for flagging
            tBackgroundMesh->set_activation_pattern( mParameters->get_lagrange_input_pattern() );

            // pick any lagrange mesh ( it really doesn't matter which one )
            Lagrange_Mesh_Base * tLagrangeMesh = mDatabase->get_lagrange_mesh_by_index( 0 );

            // loop over all active elements
            for( mtk::Cell* tCell :  aElements )
            {
                // get pointer to Background Element
                // ( the input elements are general mtk::cells, so they might not have any )
                Background_Element_Base * tElement =
                        tLagrangeMesh->get_element_by_memory_index(
                                tCell->get_memory_index_of_background_element() )
                                ->get_background_element();

                // put this element on the list
                tElement->set_refined_flag( tWorkingPattern );

                // set the minumum refinement level, which is inherited to children
                tElement->update_min_refimenent_level( aMinRefinementLevel );

                // also flag all parents
                while( tElement->get_level() > 0 )
                {
                    // get parent of this element
                    tElement = tElement->get_parent();

                    // set flag for parent of element
                    tElement->set_refined_flag( tWorkingPattern );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::perform_refinement( const bool aRefinementMode )
        {
            // refine database and remember flag
            mDatabase->perform_refinement( aRefinementMode, ! mPerformRefinementCalled );

            // remember that refinement has been called
            mPerformRefinementCalled = true;
        }

// -----------------------------------------------------------------------------

        void
        HMR::update_refinement_pattern()
        {
            mDatabase->copy_pattern(
                    mParameters->get_bspline_output_pattern(),
                    mParameters->get_bspline_input_pattern() );

            mDatabase->copy_pattern(
                    mParameters->get_lagrange_output_pattern(),
                    mParameters->get_lagrange_input_pattern() );

            // get number of bspline meshes
            uint tNumberOfBsplineMeshes = mDatabase->get_number_of_bspline_meshes();

            // update bspline meshes
            for( uint k = 0; k < tNumberOfBsplineMeshes; ++k )
            {
                // get pointer to bspline mesh
                BSpline_Mesh_Base * tMesh = mDatabase->get_bspline_mesh_by_index( k );

                if( tMesh->get_activation_pattern() == mParameters->get_bspline_input_pattern() )
                {
                    tMesh->update_mesh();
                }
            }

            // get number of bspline meshes
            uint tNumberOfLagrangeMeshes = mDatabase->get_number_of_lagrange_meshes();

            // update lagrange meshes
            for( uint k = 0; k < tNumberOfLagrangeMeshes; ++k )
            {
                // get pointer to bspline mesh
                Lagrange_Mesh_Base * tMesh = mDatabase->get_lagrange_mesh_by_index( k );

                if( tMesh->get_activation_pattern() == mParameters->get_lagrange_input_pattern() )
                {
                    tMesh->update_mesh();
                }
            }

            // set flag that this function has been called
            mUpdateRefinementCalled = true;
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Mesh >
        HMR::create_mesh()
        {
            return std::make_shared< Mesh >( mDatabase,
                                             mParameters->get_lagrange_orders().max(),
                                             mParameters->get_lagrange_output_pattern() );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Mesh >
        HMR::create_mesh( const uint & aLagrangeOrder )
        {
            return std::make_shared< Mesh >( mDatabase,
                                             aLagrangeOrder,
                                             mParameters->get_lagrange_output_pattern() );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Mesh >
        HMR::create_mesh( const uint & aLagrangeOrder, const uint & aPattern )
        {
            return std::make_shared< Mesh >( mDatabase,
                                             aLagrangeOrder,
                                             aPattern );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::create_field( const std::string & aLabel )
        {
            return this->create_field( aLabel,
                                       mParameters->get_lagrange_orders().max(),
                                       mParameters->get_bspline_orders().max() );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::create_field( const std::string & aLabel,
                           const uint        & aLagrangeOrder,
                           const uint        & aBSplineOrder )
        {
            //return mInputMesh->create_field( aLabel );
            uint tFieldIndex = mFields.size();

            // add a new field to the list
            mFields.push_back( mInputMeshes( mLagrangeOrderToInputMeshIndexMap( aLagrangeOrder ) )
                                                ->create_field( aLabel, aBSplineOrder ) );

            // return the pointer
            return mFields( tFieldIndex );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::create_field( ParameterList & aParameters )
        {
            if( mParameters->is_verbose() )
            {

                std::fprintf( stdout,"%s Loading field %s from file %s.\n\n",
                        proc_string().c_str(),
                        aParameters.get< std::string >("label").c_str(),
                        aParameters.get< std::string >("source").c_str() );
            }

            // load the field from an exodos or hdf file
            std::shared_ptr< Field > aField = this->load_field_from_file(
                    aParameters.get< std::string >("label"),
                    aParameters.get< std::string >("source"),
                    aParameters.get< sint >( "lagrange_order" ),
                    aParameters.get< sint >( "bspline_order" ) );

            // set refinement levels of field
            aField->set_min_volume_level( aParameters.get< sint >( "min_volume_refinement_level" ) );
            aField->set_max_volume_level( aParameters.get< sint >( "max_volume_refinement_level" ) );
            aField->set_min_surface_level( aParameters.get< sint >( "min_surface_refinement_level" ) );
            aField->set_max_surface_level( aParameters.get< sint >( "max_surface_refinement_level" ) );

            // return the field pointer
            return aField;
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_background_mesh_to_vtk( const std::string & aFilePath )
        {
            mDatabase->get_background_mesh()->save_to_vtk( aFilePath );
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_bsplines_to_vtk( const std::string & aFilePath )
        {
            for( uint k=0; k<mDatabase->get_number_of_lagrange_meshes(); ++k  )
            {
                // pick mesh
                if( mDatabase->get_lagrange_mesh_by_index( k )->get_activation_pattern()
                        == mParameters->get_lagrange_output_pattern() )
                {
                    // dump mesh
                    mDatabase->get_lagrange_mesh_by_index( k )
                            ->get_bspline_mesh( mDatabase->get_lagrange_mesh_by_index( k )->get_order() )->save_to_vtk( aFilePath );
                    break;
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_faces_to_vtk( const std::string & aFilePath )
        {
            for( uint k=0; k<mDatabase->get_number_of_lagrange_meshes(); ++k  )
            {
                // pick mesh
                if( mDatabase->get_lagrange_mesh_by_index( k )->get_activation_pattern()
                        == mParameters->get_lagrange_output_pattern() )
                {
                    // dump mesh
                    mDatabase->get_lagrange_mesh_by_index( k )->save_faces_to_vtk( aFilePath );

                    break;
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_edges_to_vtk( const std::string & aFilePath )
        {
            if( mParameters->get_number_of_dimensions() == 3 )
            {
                for( uint k=0; k<mDatabase->get_number_of_lagrange_meshes(); ++k  )
                {
                    // pick mesh
                    if( mDatabase->get_lagrange_mesh_by_index( k )->get_activation_pattern()
                            == mParameters->get_lagrange_output_pattern() )
                    {
                        // dump mesh
                        mDatabase->get_lagrange_mesh_by_index( k )->save_edges_to_vtk( aFilePath );

                        break;
                    }
                }
            }
            else
            {
                MORIS_ERROR( false, "save_edges_to_vtk() can only be called for 3D meshes");
            }
        }

// ----------------------------------------------------------------------------

        void
        HMR::save_mesh_to_vtk( const std::string & aFilePath )
        {
            for( uint k=0; k<mDatabase->get_number_of_lagrange_meshes(); ++k  )
            {
                // pick mesh
                if( mDatabase->get_lagrange_mesh_by_index( k )->get_activation_pattern()
                        == mParameters->get_lagrange_output_pattern() )
                {
                    // dump mesh
                    mDatabase->get_lagrange_mesh_by_index( k )->save_to_vtk( aFilePath );

                    break;
                }
            }
        }

// ----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::load_field_from_hdf5_file( const std::string & aLabel,
                                        const std::string & aFilePath,
                                        const uint          aLagrangeOrder,
                                        const uint          aBSpineOrder )
        {
            // opens an existing file with read and write access
            hid_t tFileID = open_hdf5_file( aFilePath );

            // error handler
            hid_t tStatus = 0;

            // load order from file
            uint tLagrangeOrder;
            if( aLagrangeOrder == 0 )
            {
                // assume that lagrange order is maximum of input
                tLagrangeOrder = mParameters->get_lagrange_orders().max();
            }
            else
            {
                // take passed parameter
                tLagrangeOrder = aLagrangeOrder;
            }

            uint tFieldIndex = mFields.size();

            // add a new field to the list ( zero will be overwritten )
            mFields.push_back( mInputMeshes( mLagrangeOrderToInputMeshIndexMap( tLagrangeOrder ) )->create_field(
                                                                                         aLabel, tLagrangeOrder ) );

            // get a pointer to this field
            std::shared_ptr< Field > aField = mFields( tFieldIndex );

            load_matrix_from_hdf5_file( tFileID, aLabel, aField->get_coefficients(), tStatus );

            // close hdf5 file
            close_hdf5_file( tFileID );

            uint tBSplineOrder = 0;
            uint tNumberOfCoeffs = aField->get_coefficients().length();

            for( uint k = 0; k < mDatabase->get_number_of_bspline_meshes(); ++k )
            {
                // get pointer to B-Spline mesh
                BSpline_Mesh_Base * tBMesh = mDatabase->get_bspline_mesh_by_index( k );

                // test if mesh is an input mesh
                if( tBMesh->get_activation_pattern() == mParameters->get_bspline_input_pattern() )
                {
                    if( tBMesh->get_number_of_active_basis_on_proc() == tNumberOfCoeffs )
                    {
                        tBSplineOrder = tBMesh->get_order();
                        break;
                    }
                }
            }

            MORIS_ERROR( tBSplineOrder != 0, "Could not determine order of B-Spline coefficients" );
            MORIS_ERROR( tBSplineOrder == aBSpineOrder || aBSpineOrder == 0,
                    "Specified B-Spline order does not match order of coefficients" );

            // set order of B-Splines
            aField->set_bspline_order( tBSplineOrder );

            // get number of nodes from input mesh
            uint tNumberOfNodes = mInputMeshes( mLagrangeOrderToInputMeshIndexMap( tLagrangeOrder ) )->get_num_nodes();

            // allocate field of nodes
            aField->get_node_values().set_size( tNumberOfNodes, 1 );

            // evaluate node values
            aField->evaluate_node_values();

            // return the pointer
            return aField;
        }

// ----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::load_field_from_exo_file( const std::string & aLabel,
                                       const std::string & aFilePath,
                                       const uint          aLagrangeOrder,
                                       const uint          aBSpineOrder )
        {
            // create mesh object
            mtk::Mesh * tMesh = mtk::create_mesh( MeshType::STK, aFilePath, nullptr );

            // load order from file
            uint tLagrangeOrder;

            // try to detect Lagrange order from mesh
            if( aLagrangeOrder == 0 )
            {
                // nodes of first element
                Matrix< IndexMat > tNodes =
                tMesh->get_nodes_connected_to_element_loc_inds( 0 );

                // convert number of nodes to Interpolation rank and to numeric order
                tLagrangeOrder = mtk::interpolation_order_to_uint(
                                     mtk::get_interpolation_order_from_element_nodes_and_dimension(
                                                       tNodes.length(),
                                                       mParameters->get_number_of_dimensions() ) );
            }
            else
            {
                // take passed parameter
                tLagrangeOrder = aLagrangeOrder;
            }

            // if B-Spline order is zero, assume same order as lagrange
            uint tBSplineOrder;
            if( aBSpineOrder == 0 )
            {
                tBSplineOrder = tLagrangeOrder;
            }
            else
            {
                tBSplineOrder = aBSpineOrder;
            }

            uint tFieldIndex = mFields.size();

            // get pointer to HMR Mesh
            auto tHmrMesh = mInputMeshes( mLagrangeOrderToInputMeshIndexMap( tLagrangeOrder ) );

            // add a new field to the list ( zero will be overwritten )
            mFields.push_back( tHmrMesh->create_field( aLabel, tLagrangeOrder ) );

            // get a pointer to this field
            std::shared_ptr< Field > aField = mFields( tFieldIndex );

            // load nodes from mesh
            uint tNumberOfNodes = tMesh->get_num_nodes();

            // make sure that mesh is correct
            MORIS_ASSERT( tNumberOfNodes == tHmrMesh->get_num_nodes(),
                            "Number of Nodes does not match" );

            // create array of indices for MTK interface
            Matrix< IndexMat > tIndices( 1, tNumberOfNodes  );
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                tIndices( k ) = k;
            }

            // get ref to values
            Matrix< DDRMat > & tValues = aField->get_node_values();

            // allocate nodal field
            tValues.set_size( tNumberOfNodes, 1 );

            tValues = tMesh->get_entity_field_value_real_scalar( tIndices,
                                                                 aLabel,
                                                                 EntityRank::NODE );

            // having the values, we must no rearrange them in the order of the HMR mesh.
            // Therefore, we create a map
            map< moris_id, real > tMap;
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // get ID of this node and map it with value
                tMap[ tMesh->get_glb_entity_id_from_entity_loc_index(
                        k,
                        EntityRank::NODE ) ] = tValues( k );
            }
            // make sure that field is a row matrix
            tValues.set_size( tNumberOfNodes, 1 );

            // now, we rearrange the values according to the ID of the Lagrange Mesh
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                tValues( k ) = tMap.find( tHmrMesh->get_mtk_vertex( k ).get_id() );
            }

            // finally, we set the order of the B-Spline coefficients
            aField->set_bspline_order( tBSplineOrder );

            // delete mesh pointer
            delete tMesh;

            // return the pointer
            return aField;
        }

// ----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::load_field_from_file( const std::string & aLabel,
                                   const std::string & aFilePath,
                                   const uint          aLagrangeOrder,
                                   const uint          aBSpineOrder )
        {
            // detect file type
            std::string tType = aFilePath.substr( aFilePath.find_last_of(".")+1, aFilePath.length() );

            if( tType == "hdf5" || tType == "h5" )
            {
                // assume this is a hdf file
                return this->load_field_from_hdf5_file( aLabel, aFilePath, aLagrangeOrder, aBSpineOrder );
            }
            else
            {
                // assume this is an exodus file
                return this->load_field_from_exo_file( aLabel, aFilePath, aLagrangeOrder, aBSpineOrder );
            }
        }

// ----------------------------------------------------------------------------

        void
        HMR::create_input_and_output_meshes()
        {
            // clear memory
            mInputMeshes.clear();
            mOutputMeshes.clear();

            // get orders for Lagrange meshes from patterns
            Matrix< DDUMat > tOrders;
            unique( mParameters->get_lagrange_orders(), tOrders );

            // get number of Lagrange meshes
            uint tNumberOfMeshes = tOrders.length();

            mLagrangeOrderToInputMeshIndexMap.set_size( 4, 1, MORIS_UINT_MAX );

            // create meshes
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                mLagrangeOrderToInputMeshIndexMap( tOrders( k ) ) = k;
                mInputMeshes.push_back(  this->create_mesh( tOrders( k ), mParameters->get_lagrange_input_pattern() ) );
                mOutputMeshes.push_back( this->create_mesh( tOrders( k ), mParameters->get_lagrange_output_pattern() ) );
            }
        }

// ----------------------------------------------------------------------------

        void
        HMR::perform_initial_refinement()
        {
            // get minimum refinement from parameters object
            uint tInitialRefinement = mParameters->get_minimum_initial_refimenent();

            // get pointer to background mesh
            Background_Mesh_Base * tBackMesh =  mDatabase->get_background_mesh();

            // get number of active elements on mesh
            uint tNumberOfElements = tBackMesh->get_number_of_active_elements_on_proc();

            // flag all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to background element
                Background_Element_Base * tElement = tBackMesh->get_element( e );

                // set minumum level for this element
                tElement->set_min_refimenent_level( tInitialRefinement );

                // flag this element
                tElement->put_on_refinement_queue();
            }

            // run the refiner
            this->perform_refinement( gRefinementModeBSpline );
        }

// ----------------------------------------------------------------------------

        void
        HMR::user_defined_flagging(
                bool (*aFunction)(
                        const Element                    * aElement,
                        const Cell< Matrix< DDRMat > >   & aElementLocalValues,
                              ParameterList              & aParameters ),
                        Cell< std::shared_ptr< Field > > & aFields,
                              ParameterList              & aParameters,
                        const bool                         aRefinementMode )
        {

            // remember current active scheme
            uint tActivePattern = mDatabase->get_activation_pattern();

            // define patterns
            uint tInputPattern;

            if( aRefinementMode == gRefinementModeBSpline )
            {
                tInputPattern = mParameters->get_bspline_input_pattern();
            }
            else
            {
                tInputPattern = mParameters->get_bspline_output_pattern();
            }

            // set activation pattern to input
            if( tActivePattern != tInputPattern )
            {
                // set active pattern to input mesh
                mDatabase->set_activation_pattern( tInputPattern );
            }

            // get number of fields
            uint tNumberOfFields = aFields.size();

            // create empty cell of fields
            Matrix< DDRMat> tEmpty;
            Cell< Matrix< DDRMat> > tFields( tNumberOfFields, tEmpty );

            // get number of elements from input mesh
            uint tNumberOfElements
                = mDatabase->get_background_mesh()->get_number_of_active_elements_on_proc();

            // loop over all fields
            for( auto tField: aFields )
            {
                // test if field has node values
                if( tField->get_node_values().length() == 0 )
                {
                    tField->evaluate_node_values();
                }
            }

            // loop over all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // loop over all fields
                for( uint f = 0; f<tNumberOfFields; ++f )
                {
                    // grab nodal values
                    aFields( f )->get_element_local_node_values( e, tFields( f ) );
                }

                // perform flagging test
                if( aFunction(
                        mInputMeshes( 0 )->get_lagrange_mesh()->get_element( e ),
                        tFields,
                        aParameters ) )
                {
                    // flag this element
                    this->flag_element( e );
                }
            }

            // reset activation pattern of database
            if( tActivePattern != tInputPattern )
            {
                mDatabase->set_activation_pattern( tActivePattern );
            }
        }

// ----------------------------------------------------------------------------

        uint
        HMR::get_mesh_index( const uint aOrder, const uint aPattern )
        {
            uint aIndex = MORIS_UINT_MAX;

            // find correct output mesh
            uint tNumberOfMeshes = mDatabase->get_number_of_lagrange_meshes();
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = mDatabase->get_lagrange_mesh_by_index( k );

                if(        tMesh->get_activation_pattern() == aPattern
                        && tMesh->get_order() == aOrder )
                {
                    aIndex = k;
                    break;
                }
            }

            return aIndex;
        }

// ----------------------------------------------------------------------------

        void
        HMR::get_candidates_for_refinement(
                Cell< mtk::Cell* > & aCandidates,
                const bool           aRefinementMode,
                const uint           aMaxLevel )
        {
            // reset candidate list
            aCandidates.clear();

            uint tPattern;

            if( aRefinementMode == gRefinementModeBSpline )
            {

                tPattern =  mParameters->get_bspline_input_pattern();
            }
            else
            {
                tPattern =  mParameters->get_lagrange_input_pattern();
            }

            // make sure that input pattern is active
            mDatabase->set_activation_pattern( tPattern );

            // get pointer to background mesh
            Background_Mesh_Base * tBackgroundMesh = mDatabase->get_background_mesh();

            // pick first Lagrange mesh on input pattern
            // fixme: add option to pick another one
            Lagrange_Mesh_Base * tMesh = mDatabase->get_lagrange_mesh_by_index( 0 );

            // get max level of this mesh
            //uint tMaxLevel = std::min( tBackgroundMesh->get_max_level(), aMaxLevel );

            // counter for elements
            uint tCount = 0;

            // loop over all levels
            for( uint l = 0; l < aMaxLevel; ++l )
            {
                Cell< Background_Element_Base * > tBackgroundElements;

                tBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tBackgroundElements );

                // element must be active or refined
                for( Background_Element_Base * tElement : tBackgroundElements )
                {
                    if( ( tElement->is_active( tPattern ) ||  tElement->is_refined( tPattern ) )
                            && ! tElement->is_padding() )
                    {
                        // increment counter
                        ++tCount;
                    }
                }
            }

            // allocate memory for output
            aCandidates.resize( tCount, nullptr );

            // reset counter
            tCount = 0;
            // loop over all levels
            for( uint l=0; l<aMaxLevel; ++l )
            {
                Cell< Background_Element_Base * > tBackgroundElements;
                tBackgroundMesh->collect_elements_on_level_within_proc_domain( l, tBackgroundElements );

                // element must be active or refined
                for(  Background_Element_Base * tElement : tBackgroundElements )
                {
                    if( ( tElement->is_active( tPattern ) ||  tElement->is_refined( tPattern ) )
                            && ! tElement->is_padding() )
                    {
                        aCandidates( tCount++ ) = tMesh->get_element_by_memory_index( tElement->get_memory_index() );
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        uint
        HMR::flag_volume_and_surface_elements(
                const std::shared_ptr<Field> aScalarField )
        {
            // the funciton returns the number of flagged elements
            uint aElementCounter = 0;

            // create geometry engine
            gen::Geometry_Engine tRefMan;

            // candidates for refinement
            Cell< mtk::Cell*  > tCandidates;

            // elements to be flagged for refinement
            Cell< mtk::Cell* > tRefinementList;

            // get candidates for surface
            this->get_candidates_for_refinement( tCandidates,
                    gRefinementModeBSpline,
                    aScalarField->get_max_surface_level() );

            // call refinement manager and get intersected cells
            tRefMan.find_cells_intersected_by_levelset( tRefinementList,
                    tCandidates,
                    aScalarField->get_node_values() );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in HMR
            this->flag_elements( tRefinementList, aScalarField->get_min_surface_level() );

            // get candidates from volume
            this->get_candidates_for_refinement( tCandidates,
                    aScalarField->get_max_volume_level() );

            // call refinement manager and get volume cells
            tRefMan.find_cells_within_levelset( tRefinementList,
                    tCandidates,
                    aScalarField->get_node_values() );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in database
            this->flag_elements( tRefinementList, aScalarField->get_min_volume_level()  );

            // return number of flagged elements
            return aElementCounter;
        }

// -----------------------------------------------------------------------------

        uint
        HMR::flag_surface_elements( const std::shared_ptr<Field> aScalarField )
        {
            // the funciton returns the number of flagged elements
            uint aElementCounter = 0;

            // create geometry engine
            gen::Geometry_Engine tRefMan;

            // candidates for refinement
            Cell< mtk::Cell* > tCandidates;

            // elements to be flagged for refinement
            Cell< mtk::Cell* > tRefinementList;

            // get candidates for surface
            this->get_candidates_for_refinement(
                    tCandidates,
                    gRefinementModeBSpline,
                    aScalarField->get_max_surface_level() );

            // call refinement manager and get intersected cells
            tRefMan.find_cells_intersected_by_levelset( tRefinementList,
                    tCandidates,
                    aScalarField->get_node_values() );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in HMR
            this->flag_elements( tRefinementList, aScalarField->get_min_surface_level() );

            // return number of flagged elements
            return aElementCounter;
        }

// ----------------------------------------------------------------------------

        /**
         * needed for tutorials
         */
        void
        HMR::perform_refinement_and_map_fields()
        {
            // - - - - - - - - - - - - - - - - - - - - - -
            // step 0: perform simple refinement
            // - - - - - - - - - - - - - - - - - - - - - -

            // in the tutorial, lagrange and B-Spline are the same refinement
            this->perform_refinement( gRefinementModeBSpline );

            // create union of input and output
            mDatabase->create_union_pattern();

            // finalize database
            this->finalize();

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 1: find out which orders are needed
            // - - - - - - - - - - - - - - - - - - - - - -

            // number of input fields
            uint tNumberOfFields = mFields.size();

            // counter
            uint tCount = 0;

            // container for orders of fields
            Matrix< DDUMat > tInputFieldOrders( 2*tNumberOfFields, 1 );

            // loop over all fields
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                tInputFieldOrders( tCount++ ) = mFields( f )->get_bspline_order();
                tInputFieldOrders( tCount++ ) = mFields( f )->get_lagrange_order();
            }

            // chop container
            tInputFieldOrders.resize( tCount, 1 );

            // make orders unique
            Matrix< DDUMat > tMeshOrders;
            unique( tInputFieldOrders, tMeshOrders );

            uint tNumberOfMappers = tMeshOrders.length();

            // create map for mappers
            Matrix< DDUMat > tMapperIndex( gMaxBSplineOrder+1, 1, MORIS_UINT_MAX );

            for( uint k = 0; k<tNumberOfMappers; ++k )
            {
                tMapperIndex( tMeshOrders( k ) ) = k;
            }

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 2: create union meshes and mappers
            // - - - - - - - - - - - - - - - - - - - - - -

            Cell< std::shared_ptr< Mesh > > tUnionMeshes;
            Cell< std::shared_ptr< Mesh > > tInputMeshes;
            Cell< mapper::Mapper * > tMappers( tNumberOfMappers, nullptr );

            for( uint m=0; m<tNumberOfMappers; ++m )
            {
                // get pointer to input mesh
                tInputMeshes.push_back( this->create_mesh(
                        tMeshOrders( m ),
                        mParameters->get_lagrange_input_pattern() ) );

                // create union mesh from HMR object
                tUnionMeshes.push_back( this->create_mesh(
                        tMeshOrders( m ),
                        mParameters->get_union_pattern() ) );

                // create mapper
                tMappers( m ) = new mapper::Mapper( tUnionMeshes( m ) );
            }

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 3: map and project fields
            // - - - - - - - - - - - - - - - - - - - - - -

            for( uint f=0; f<tNumberOfFields; ++f )
            {

                // get pointer to input field
                std::shared_ptr< Field > tInputField = mFields( f );

                // get order
                uint tBSplineOrder = tInputField->get_bspline_order();

                // get index of mapper
                uint m = tMapperIndex( tBSplineOrder );

                // get pointer to field on union mesh
                std::shared_ptr< Field > tUnionField =  tUnionMeshes( m )->create_field(
                        tInputField->get_label(),
                        tBSplineOrder );


                if( tInputField->get_lagrange_order() >= tBSplineOrder )
                {
                    // interpolate field onto union mesh
                    mDatabase->interpolate_field(
                            mParameters->get_lagrange_input_pattern(),
                            tInputField,
                            mParameters->get_union_pattern(),
                            tUnionField );
                }
                else
                {
                    // first, project field on mesh with correct order
                    std::shared_ptr< Field > tTemporaryField =
                            tInputMeshes( m )->create_field(
                                    tInputField->get_label(),
                                    tBSplineOrder );

                    mDatabase->change_field_order(
                            tInputField, tTemporaryField );

                    // now, interpolate this field onto the inion
                    mDatabase->interpolate_field(
                            mParameters->get_lagrange_input_pattern(),
                            tTemporaryField,
                            mParameters->get_union_pattern(),
                            tUnionField );
                }

                // perform mapping
                tMappers( m )->perform_mapping(
                        tInputField->get_label(),
                        EntityRank::NODE,
                        tInputField->get_label(),
                        tUnionField->get_bspline_rank() );

                // a small sanity test
                MORIS_ASSERT(  tUnionField->get_coefficients().length()
                        == tUnionMeshes( m )->get_num_entities(
                                mtk::order_to_entity_rank( tBSplineOrder ) ),
                                "Number of B-Splines does not match" );

                // get pointer to output mesh
                std::shared_ptr< Mesh >  tOutputMesh = this->create_mesh(
                        tInputField->get_lagrange_order(),
                        mParameters->get_lagrange_output_pattern() );

                // create output field
                std::shared_ptr< Field >  tOutputField =
                        tOutputMesh->create_field(
                                tInputField->get_label(),
                                tBSplineOrder );

                // move coefficients to output field
                // fixme: to be tested with Eigen also
                tOutputField->get_coefficients() = std::move( tUnionField->get_coefficients() );

                // allocate nodes for output
                tOutputField->get_node_values().set_size( tOutputMesh->get_num_nodes(), 1 );

                // evaluate nodes
                tOutputField->evaluate_node_values();

                // make this field point to the output mesh
                tInputField->change_mesh(
                        tOutputField->get_mesh(),
                        tOutputField->get_field_index() );
            }

            // delete mappers
            for( mapper::Mapper * tMapper : tMappers )
            {
                delete tMapper;
            }
        }
// ----------------------------------------------------------------------------

        void
        HMR::map_field_to_output( std::shared_ptr< Field > aField )
        {
            // grab orders of meshes
            uint tLagrangeOrder = aField->get_lagrange_order();
            uint tBSplineOrder =  aField->get_bspline_order();

            // get order of Union Mesh
            uint tOrder = std::max( tLagrangeOrder, tBSplineOrder );

            // create union pattern
            mDatabase->create_union_pattern();

            // create union mesh
            std::shared_ptr< Mesh > tUnionMesh =  this->create_mesh(
                    tOrder, mParameters->get_union_pattern() );

            // create union field
            std::shared_ptr< Field > tUnionField =  tUnionMesh->create_field(
                    aField->get_label(),
                    aField->get_bspline_order() );

            if( tLagrangeOrder >=  tBSplineOrder )
            {
                // interpolate field onto union mesh
                mDatabase->interpolate_field(
                        mParameters->get_lagrange_input_pattern(),
                        aField,
                        mParameters->get_union_pattern(),
                        tUnionField );
            }
            else
            {
                // mesh the input field is based on
                std::shared_ptr< Mesh > tInputMesh = this->create_mesh(
                        tOrder, mParameters->get_lagrange_input_pattern() );

                // first, project field on mesh with correct order
                std::shared_ptr< Field > tTemporaryField =
                        tInputMesh->create_field(
                                aField->get_label(),
                                tBSplineOrder );

                mDatabase->change_field_order( aField, tTemporaryField );

                // now, interpolate this field onto the union
                mDatabase->interpolate_field(
                        mParameters->get_lagrange_input_pattern(),
                        tTemporaryField,
                        mParameters->get_union_pattern(),
                        tUnionField );
            }

            // create mapper
            mapper::Mapper tMapper( tUnionMesh );

            // project field to union
            tMapper.perform_mapping(
                                    aField->get_label(),
                                    EntityRank::NODE,
                                    aField->get_label(),
                                    tUnionField->get_bspline_rank() );

            // a small sanity test
            MORIS_ASSERT(  tUnionField->get_coefficients().length()
                    == tUnionMesh->get_num_entities(
                            mtk::order_to_entity_rank( tBSplineOrder ) ),
                            "Number of B-Splines does not match" );

            // get pointer to output mesh
            std::shared_ptr< Mesh >  tOutputMesh = this->create_mesh(
                    tLagrangeOrder,
                    mParameters->get_lagrange_output_pattern() );

            // create output field
            std::shared_ptr< Field >  tOutputField =
                    tOutputMesh->create_field(
                            aField->get_label(),
                            tBSplineOrder );

            // move coefficients to output field
            tOutputField->get_coefficients() = std::move( tUnionField->get_coefficients() );

            // allocate nodes for output
            tOutputField->get_node_values().set_size( tOutputMesh->get_num_nodes(), 1 );

            // evaluate nodes
            tOutputField->evaluate_node_values();

            // make this field point to the output mesh
            aField->change_mesh(
                    tOutputField->get_mesh(),
                    tOutputField->get_field_index() );
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
