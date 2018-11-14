

/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

// see http://pubs.opengroup.org/onlinepubs/7908799/xsh/dlfcn.h.html
#include "dlfcn.h"

#include "typedefs.hpp"
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

            mDatabase ->calculate_t_matrices_for_input();
        }

// -----------------------------------------------------------------------------

        void
        HMR::finalize()

        {
            // if mesh has not been refined, copy input to output before finalizing
            if( ! mDatabase->have_refined_at_least_one_element() )
            {
                // select output pattern
                mDatabase->set_activation_pattern(
                        mParameters->get_output_pattern()  );

                // copy input to output
                mDatabase->get_background_mesh()->copy_pattern(
                        mParameters->get_input_pattern(),
                        mParameters->get_output_pattern() );

                // update database
                mDatabase->update_meshes();
            }

            // finish database
            mDatabase->finalize();

        }
// -----------------------------------------------------------------------------

        void
        HMR::load_output_pattern_from_path( const std::string & aPath )
        {
            mDatabase->load_pattern_from_hdf5_file(
                    mParameters->get_output_pattern(),
                    aPath );
        }

// -----------------------------------------------------------------------------
        void
        HMR::save_to_exodus( const std::string & aPath, const double aTimeStep )
        {
            /*if( ! mPerformRefinementCalled )
            {
                this->save_to_exodus(
                        mParameters->get_input_pattern(),
                        aPath,
                        aTimeStep );
            }
            else
            {*/
            this->save_to_exodus(
                    mParameters->get_output_pattern(),
                    aPath,
                    aTimeStep );
            //}

        }

// -----------------------------------------------------------------------------

        void
        HMR::save_last_step_to_exodus(
                            const std::string & aPath,
                            const double aTimeStep )
        {
            MORIS_ERROR( ! mUpdateRefinementCalled,
                    "HMR does not feel comfortable with you calling save_last_step_to_exodus() after you have overwritten the input pattern using update_refinement_pattern()");

            this->save_to_exodus( mParameters->get_input_pattern(),
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
            tHDF5.save_refinement_pattern(
                    tBackgroundMesh,
                    mParameters->get_output_pattern() );

            if( tActivePattern != tBackgroundMesh->get_activation_pattern() )
            {
                tBackgroundMesh->set_activation_pattern( tActivePattern );
            }

            // close hdf5 file
            tHDF5.close();
        }

// -----------------------------------------------------------------------------

        void
        HMR::save_coeffs_to_binary_files(
                const std::string & aFilePath )
        {
            // get number of meshes
            uint tNumberOfLagrangeMeshes = mDatabase->get_number_of_lagrange_meshes();

            // loop over all meshes
            for( uint m=0; m<tNumberOfLagrangeMeshes; ++m )
            {
                // get pointer to mesh
                Lagrange_Mesh_Base * tMesh = mDatabase->get_lagrange_mesh_by_index( m );

                // test if mesh links to output pattern
                if( tMesh->get_activation_pattern() == mParameters->get_output_pattern() )
                {
                    uint tNumberOfBSplineMeshes = tMesh->get_number_of_bspline_meshes();

                    // loop over all B-Spline meshes
                    for( uint n=0; n<tNumberOfBSplineMeshes; ++n )
                    {

                        // get pointer to bspoine mesh
                        BSpline_Mesh_Base * tBMesh = tMesh->get_bspline_mesh( n );

                        if( tBMesh != NULL )
                        {
                            // get order of B-Spline mesh
                            uint tOrder = tBMesh->get_order();

                            // calculate file path
                            std::string tFilePath =
                                    aFilePath.substr(0,aFilePath.find_last_of(".")) // base path

                                    // get order of lagrange mesh
                                    + "_" + std::to_string( tMesh->get_order() )

                            // get order of bspline mesh
                            + "_" + std::to_string( tOrder )

                            // finish path
                            +  aFilePath.substr( aFilePath.find_last_of("."), aFilePath.length() );

                            tMesh->save_coeffs_to_binary_file( tOrder, tFilePath );
                        }
                    }
                }
            }
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
                if( tMesh->get_activation_pattern() == mParameters->get_output_pattern() )
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
                const uint               aMinRefinementLevel )
        {

            // get  working pattern
            uint tWorkingPattern = mParameters->get_working_pattern();

            // get pointer to background mesh
            Background_Mesh_Base * tBackgroundMesh = mDatabase->get_background_mesh();

            // use output pattern
            tBackgroundMesh->set_activation_pattern( mParameters->get_input_pattern() );

            // pick any lagrange mesh ( it really doesn't matter which one )
            Lagrange_Mesh_Base * tLagrangeMesh = mDatabase->get_lagrange_mesh_by_index( 0 );

            // loop over all active elements
            for( mtk::Cell* tCell :  aElements )
            {
                // get pointer to Background Element
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
        HMR::perform_refinement()
        {
            // refine database and remember flag
            mDatabase->perform_refinement( ! mPerformRefinementCalled );

            // remember that refinement has been called
            mPerformRefinementCalled = true;
        }

// -----------------------------------------------------------------------------

        void
        HMR::update_refinement_pattern()
        {
            mDatabase->copy_pattern(
                    mParameters->get_output_pattern(),
                    mParameters->get_input_pattern() );

            // get number of bspline meshes
            uint tNumberOfBsplineMeshes = mDatabase->get_number_of_bspline_meshes();

            // update bspline meshes
            for( uint k=0; k<tNumberOfBsplineMeshes; ++k )
            {
                // get pointer to bspline mesh
                BSpline_Mesh_Base * tMesh = mDatabase->get_bspline_mesh_by_index( k );

                if( tMesh->get_activation_pattern() == mParameters->get_input_pattern() )
                {
                    tMesh->update_mesh();
                }
            }

            // get number of bspline meshes
            uint tNumberOfLagrangeMeshes = mDatabase->get_number_of_lagrange_meshes();

            // update lagrange meshes
            for( uint k=0; k<tNumberOfLagrangeMeshes; ++k )
            {
                // get pointer to bspline mesh
                Lagrange_Mesh_Base * tMesh = mDatabase->get_lagrange_mesh_by_index( k );

                if( tMesh->get_activation_pattern() == mParameters->get_input_pattern() )
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
                                             mParameters->get_output_pattern() );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Mesh >
        HMR::create_mesh( const uint & aLagrangeOrder )
        {
            return std::make_shared< Mesh >( mDatabase,
                                             aLagrangeOrder,
                                             mParameters->get_output_pattern() );
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
        HMR::create_field(
                const std::string & aLabel,
                const uint        & aLagrangeOrder,
                const uint        & aBSplineOrder )
        {
            //return mInputMesh->create_field( aLabel );
            uint tFieldIndex = mFields.size();

            // add a new field to the list
            mFields.push_back( mInputMeshes(
                    mLagrangeOrderToInputMeshIndexMap( aLagrangeOrder ) )
                    ->create_field( aLabel, aBSplineOrder ) );

            // return the pointer
            return mFields( tFieldIndex );
        }

// -----------------------------------------------------------------------------

        std::shared_ptr< Field>
        HMR::map_field_on_mesh(
                std::shared_ptr< Field > aField,
                std::shared_ptr< Mesh >  aMesh )
        {
            tic tTimer;

            //mDatabase->get_background_mesh()->save_to_vtk("BG.vtk");
            //mDatabase->get_bspline_mesh_by_index(1)->save_to_vtk("BS.vtk");

            // create pointer to output field
            auto aOutField = aMesh->create_field( aField->get_label(), aField->get_bspline_order() );

            // create a temporary union mesh
            Mesh * tUnionMesh = new Mesh( mDatabase,
                    mtk::interpolation_order_to_uint( aField->get_interpolation_order() ),
                    mParameters->get_union_pattern() );

            // calculate T-Matrix for union mesh
            tUnionMesh->get_lagrange_mesh()->calculate_t_matrix( aField->get_bspline_order() );

            // create a temporary union field
            auto tUnionField = tUnionMesh->create_field( aField->get_label(), aField->get_bspline_order() );

            // interpolate input field to union
            mDatabase->interpolate_field(
                    mParameters->get_input_pattern(),
                    aField,
                    mParameters->get_union_pattern(),
                    tUnionField );

            // create IWG object
            moris::fem::IWG_L2 tIWG;

            // fixme: #ADOFORDERHACK
            moris::MSI::gAdofOrderHack = aField->get_bspline_order();

            // create model
            mdl::Model tModel(
                     tUnionMesh,
                     & tIWG );

            // set order of dofs
            //tModel.set_dof_order( aField->get_bspline_order() );

            // set weak bcs
            tModel.set_weak_bcs( tUnionField->get_node_values() );

            // solve problem
            tModel.solve( aOutField->get_coefficients() );

            // delete the pointer to the union mesh
            delete tUnionMesh;

            // finally, the node values on the output mesh are evaluated
            aOutField->evaluate_node_values();

            // print output if verbose level is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                std::fprintf( stdout,"%s Performed L2 projection.\n               Operation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( double ) tElapsedTime / 1000 );
            }

            // return output mesh
            return aOutField;
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
            this->get_candidates_for_refinement(
                    tCandidates,
                    aScalarField->get_max_surface_level() );

            // call refinement manager and get intersected cells
            tRefMan.find_cells_intersected_by_levelset(
                    tRefinementList,
                    tCandidates,
                    aScalarField->get_node_values() );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in HMR
            this->flag_elements( tRefinementList, aScalarField->get_min_surface_level() );

            // get candidates from volume
            this->get_candidates_for_refinement(
                    tCandidates,
                    aScalarField->get_max_volume_level() );


            // call refinement manager and get volume cells
            tRefMan.find_cells_within_levelset(
                    tRefinementList,
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
        HMR::flag_surface_elements(
                const std::shared_ptr<Field> aScalarField )
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
            this->get_candidates_for_refinement( tCandidates, aScalarField->get_max_surface_level() );


            // call refinement manager and get intersected cells
            tRefMan.find_cells_intersected_by_levelset(
                    tRefinementList,
                    tCandidates,
                    aScalarField->get_node_values() );

            // add length of list to counter
            aElementCounter += tRefinementList.size();

            // flag elements in HMR
            this->flag_elements( tRefinementList, aScalarField->get_min_surface_level() );

            // return number of flagged elements
            return aElementCounter;
        }

// -----------------------------------------------------------------------------

        void
        HMR::get_candidates_for_refinement(
                Cell< mtk::Cell* > & aCandidates,
                const uint           aMaxLevel )
        {
            // reset candidate list
            aCandidates.clear();

            // make sure that input pattern is active
            mDatabase->set_activation_pattern( mParameters->get_input_pattern() );

            // get pointer to background mesh
            Background_Mesh_Base * tBackgroundMesh = mDatabase->get_background_mesh();

            // pick first Lagrange mesh on input pattern
            // fixme: add option to pick another one
            Lagrange_Mesh_Base * tMesh = mDatabase->get_lagrange_mesh_by_index( 0 );

            auto tPattern = mParameters->get_input_pattern();

            // make sure that this mesh uses correct pattern
            MORIS_ASSERT( tMesh->get_activation_pattern() ==  tPattern,
                    "wrong pattern picked for get_candidates_for_refinement()");

            // get max level of this mesh
            //uint tMaxLevel = std::min( tBackgroundMesh->get_max_level(), aMaxLevel );

            // counter for elements
            uint tCount = 0;

            // loop over all levels
            for( uint l=0; l<aMaxLevel; ++l )
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
                        aCandidates( tCount++ )
                                = tMesh->get_element_by_memory_index(
                                        tElement->get_memory_index() );
                    }
                }
            }


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
                        == mParameters->get_output_pattern() )
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
                        == mParameters->get_output_pattern() )
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
                            == mParameters->get_output_pattern() )
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
                        == mParameters->get_output_pattern() )
                {
                    // dump mesh
                    mDatabase->get_lagrange_mesh_by_index( k )->save_to_vtk( aFilePath );

                    break;
                }
            }
        }

// ----------------------------------------------------------------------------

        void
        HMR::perform_refinement_and_map_fields()
        {

            // perform refinement
            this->perform_refinement();

            // calculate T-Matrices, edges and faces
            this->finalize();

            // get number of fields
            uint tNumberOfFields = mFields.size();

            // loop over all fields
            for( uint k=0; k<tNumberOfFields; ++k )
            {
                // get index of output mesh
                uint tOutIndex = mLagrangeOrderToInputMeshIndexMap( mFields( k )->get_lagrange_order() );

                /*
                std::shared_ptr< Mesh > tMesh = mOutputMeshes( tOutIndex );

                Lagrange_Mesh_Base* tLMesh = tMesh->get_lagrange_mesh();

                uint tNumNodes = tLMesh->get_number_of_nodes_on_proc();

                for( uint i=0; i<tNumNodes; ++i )
                {
                    Basis* tNode = tLMesh->get_node_by_index( i );



                 print( tNode->get_interpolation( 1 )->get_ids(), "IDs" );

                }

                std::cout << "Exiting MORIS" << std::endl;
                exit( 0 ); */

                // perform projection
                auto tOutField = this->map_field_on_mesh(
                        mFields( k ), mOutputMeshes( tOutIndex ) );

                // point mesh in database to new data
                mFields( k )->change_mesh(
                        mOutputMeshes( tOutIndex )->get_lagrange_mesh(),
                        tOutField->get_field_index() );
            }

        }

// ----------------------------------------------------------------------------

        std::shared_ptr< Field >
        HMR::load_field_from_hdf5_file(
                const std::string & aFilePath,
                const uint          aLagrangeOrder,
                const uint          aBSpineOrder )
        {
            // figure out Order from field

            // make path parallel
            std::string tFilePath = parallelize_path( aFilePath );

            // opens an existing file with read and write access
            hid_t tFileID = H5Fopen(
                    tFilePath.c_str(),
                    H5F_ACC_RDWR,
                    H5P_DEFAULT);

            // error handler
            hid_t tStatus;


            // load order from file
            uint tLagrangeOrder;
            if( aLagrangeOrder == 0 )
            {
                // try to load value from HDF5 file
                load_scalar_from_hdf5_file(
                        tFileID,
                        "LagrangeOrder",
                        tLagrangeOrder,
                        tStatus );
            }
            else
            {
                // take passed parameter
                tLagrangeOrder = aLagrangeOrder;
            }

            // close file
            tStatus = H5Fclose( tFileID );

            uint tFieldIndex = mFields.size();

            // add a new field to the list ( zero will be overwritten )
            mFields.push_back(
                    mInputMeshes( mLagrangeOrderToInputMeshIndexMap( tLagrangeOrder ) )->create_field( "", tLagrangeOrder ) );

            // get a pointer to this field
            std::shared_ptr< Field > aField = mFields( tFieldIndex );

            // load data
            aField->load_field_from_hdf5( aFilePath, aBSpineOrder );

            // return the pointer
            return aField;

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
                mInputMeshes.push_back(  this->create_mesh( tOrders( k ), mParameters->get_input_pattern() ) );
                mOutputMeshes.push_back( this->create_mesh( tOrders( k ), mParameters->get_output_pattern() ) );
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
            this->perform_refinement();
        }

// ----------------------------------------------------------------------------

        void
        HMR::user_defined_flagging(
                bool (*aFunction)(
                        const Element                    * aElement,
                        const Cell< Matrix< DDRMat > >   & aElementLocalValues,
                              ParameterList              & aParameters ),
                        Cell< std::shared_ptr< Field > > & aFields,
                              ParameterList              & aParameters  )
        {
            // remember current active scheme
            uint tActivePattern = mDatabase->get_activation_pattern();

            // set activation pattern to input
            if( tActivePattern != mParameters->get_input_pattern() )
            {
                // set active pattern to input mesh
                mDatabase->set_activation_pattern( mParameters->get_input_pattern() );
            }

            // get number of fields
            uint tNumberOfFields = aFields.size();

            // create empty cell of fields
            Matrix< DDRMat> tEmpty;
            Cell< Matrix< DDRMat> > tFields( tNumberOfFields, tEmpty );

            // get number of elements from input mesh
            uint tNumberOfElements
                = mDatabase->get_background_mesh()->get_number_of_active_elements_on_proc();

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
                    std::cout << "Flag element " << e << std::endl;
                }
            }

            // reset activation pattern of database
            if( tActivePattern != mParameters->get_input_pattern() )
            {
                mDatabase->set_activation_pattern( tActivePattern );
            }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
