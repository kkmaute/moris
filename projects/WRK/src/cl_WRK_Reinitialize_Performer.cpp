/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Reinitialize_Performer.cpp
 *
 */

#include "cl_WRK_Reinitialize_Performer.hpp"

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_File.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "HMR_Globals.hpp"

#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_sort.hpp"

namespace moris
{
    namespace wrk
    {
        //------------------------------------------------------------------------------

        Reinitialize_Performer::Reinitialize_Performer( std::shared_ptr< Library_IO > aLibrary )
                : mLibrary( aLibrary )
        {
            // get the parameter lists
            ModuleParameterList tMORISParameterList = aLibrary->get_parameters_for_module( Parameter_List_Type::MORISGENERAL );
            ModuleParameterList tMSIParameterList = aLibrary->get_parameters_for_module( Parameter_List_Type::MSI );

            mAdofMeshIndex = tMSIParameterList( 0 )( 0 ).get< moris::sint >( tMORISParameterList( 2 )( 0 ).get< std::string >( "dof_type" ) );

            // get the adv field name that wll be reinitialized
            mADVFiledName = tMORISParameterList( 2 )( 0 ).get< std::string >( "adv_field" );

            // get msi string to dof type map
            moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                    moris::MSI::get_msi_dof_type_map();

            // get the quantity dof type from parameter list
            string_to_cell(
                    tMORISParameterList( 2 )( 0 ).get< std::string >( "dof_type" ),
                    mDofTypes,
                    tMSIDofTypeMap );

            mReinitializationFrequency = tMORISParameterList( 2 )( 0 ).get< sint >( "reinitialization_frequency" );

            // get the mesh output info
            mOutputMeshFile = tMORISParameterList( 2 )( 0 ).get< std::string >( "output_mesh_file" );
            mTimeOffset     = tMORISParameterList( 2 )( 0 ).get< real >( "time_offset" );
        }

        //------------------------------------------------------------------------------

        void
        Reinitialize_Performer::perform(
                moris::Cell< std::shared_ptr< hmr::HMR > >&            aHMRPerformers,
                moris::Cell< std::shared_ptr< ge::Geometry_Engine > >& aGENPerformer,
                moris::Cell< std::shared_ptr< mtk::Mesh_Manager > >&   aMTKPerformer,
                moris::Cell< std::shared_ptr< mdl::Model > >           aMDLPerformer )
        {
            // Tracer to trace the time
            Tracer tTracer( "WRK", "Reinitialize ADVs", "Perform Reinitialize" );

            // initialize and populate the fields
            moris::Cell< std::shared_ptr< mtk::Field > > tGENFields;
            tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

            // find the index of the desired adv field that will be reinitialized
            auto itr = std::find_if( tGENFields.begin(), tGENFields.end(), [ & ]( std::shared_ptr< mtk::Field > const & aFiled )    //
                    { return aFiled->get_label() == mADVFiledName; } );

            // find the index of the adv field
            moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );

            // get the the adv discretization mesh index
            uint tDiscretizationMeshIndex = ( *itr )->get_discretization_mesh_index();

            // get interpolation mesh from mesh pair
            moris::mtk::Mesh* tTargetMesh = ( *itr )->get_mesh_pair().get_interpolation_mesh();

            // get the solution field and get a matrix of the solutions
            // generate a cell containing the indices of the bspline coefficients
            // since indices are consecutive and they start from 0
            moris::Cell< moris_index > tLocalCoeffIndices( tTargetMesh->get_num_entities( EntityRank::BSPLINE ) );
            std::iota( tLocalCoeffIndices.begin(), tLocalCoeffIndices.end(), 0 );

            moris::sol::Dist_Vector* tPartialSolutionVector = aMDLPerformer( 0 )->get_solver_interface()->get_solution_vector( mDofTypes, tLocalCoeffIndices );
            tPartialSolutionVector->extract_copy( mCoefficients );

            // delete the pointer as it is not needed anymore
            delete tPartialSolutionVector;

            // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
            std::shared_ptr< mtk::Field_Discrete > tFieldSource = std::make_shared< mtk::Field_Discrete >( aMTKPerformer( 0 )->get_mesh_pair( 0 ), mAdofMeshIndex );

            // unlock fields and set the coeff
            tFieldSource->unlock_field();
            tFieldSource->set_coefficients( mCoefficients );

            // compute the nodal values based on the coeff
            tFieldSource->compute_nodal_values();

            // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
            std::shared_ptr< mtk::Field_Discrete > tFieldTarget = std::make_shared< mtk::Field_Discrete >( ( *itr )->get_mesh_pair(), tDiscretizationMeshIndex );
            tFieldTarget->set_label( mADVFiledName );

            // set the nodal values
            tFieldTarget->unlock_field();
            tFieldTarget->set_values( tFieldSource->get_values() );

            // invoke the mapper and map to the target field
            mtk::Mapper tMapper;
            tFieldTarget->unlock_field();
            tMapper.map_input_field_to_output_field_2( tFieldTarget.get() );

            // compute the nodal value
            tFieldTarget->compute_nodal_values();

            // get the coefficents and store them
            mCoefficients = tFieldTarget->get_coefficients();

            // clip the values and
            this->impose_upper_lower_bound( aGENPerformer, tFieldTarget.get() );

            // replace the newly constructed field
            tGENFields( tADVFieldIndex ) = tFieldTarget;

            // store the fields
            mMTKFields = tGENFields;

            // output the fields if asked
            if ( mOutputMeshFile != "" )
            {
                this->output_fields( tFieldTarget.get(), tFieldSource.get(), mOutputMeshFile );
            }
        }

        //------------------------------------------------------------------------------

        moris::sint
        Reinitialize_Performer::get_reinitialization_frequency() const
        {
            return mReinitializationFrequency;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > const &
        Reinitialize_Performer::get_coefficients() const
        {
            return mCoefficients;
        }
        //------------------------------------------------------------------------------
        void
        Reinitialize_Performer::impose_upper_lower_bound( moris::Cell< std::shared_ptr< ge::Geometry_Engine > >& aGENPerformer, mtk::Field* aField )
        {
            // lower bound and upper bound are defined on proc 0 and they need to be communicated to other
            // Note:  we make an assumption that all the lower bounds and upper bounds are equal

            // initialize  the upper and lower abound
            moris::real tLowerBound;
            moris::real tUpperBound;

            // assign the values on processor 0
            if ( par_rank() == 0 )
            {
                tLowerBound = aGENPerformer( 0 )->get_lower_bounds()( 0 );
                tUpperBound = aGENPerformer( 0 )->get_upper_bounds()( 0 );
            }    // Bcast the values to other processeors
            MPI_Bcast( &tLowerBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
            MPI_Bcast( &tUpperBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

            // clip the values of the adv
            for ( uint iADV = 0; iADV < mCoefficients.numel(); iADV++ )
            {
                mCoefficients( iADV ) = std::max( tLowerBound, std::min( mCoefficients( iADV ), tUpperBound ) );
            }

            // update the field based on the newly clipped coeff
            aField->unlock_field();
            aField->set_coefficients( mCoefficients );
            aField->compute_nodal_values();
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< mtk::Field > >
        Reinitialize_Performer::get_mtk_fields() const
        {
            return mMTKFields;
        }

        //------------------------------------------------------------------------------

        void
        Reinitialize_Performer::output_fields( mtk::Field* aTarget, mtk::Field* aSource, std::string aExoFileName ) const
        {
            //
            Tracer tTracer( "WRK", "Reinitialize ADVs", "Outputting Fields" );
            // time shift
            real tTimeShift = 0.0;

            if ( mTimeOffset > 0 )
            {
                // get optimization iteration
                uint tOptIter = gLogger.get_opt_iteration();

                // set name
                std::string tOptIterStrg = std::to_string( tOptIter );
                aExoFileName += ".e-s." + std::string( 4 - tOptIterStrg.length(), '0' ) + tOptIterStrg;

                // determine time shift
                tTimeShift = tOptIter * mTimeOffset;
            }
            // call the lagrange mesh
            moris::mtk::Mesh* tMesh = aTarget->get_mesh_pair().get_interpolation_mesh();
            // Write mesh
            mtk::Writer_Exodus tWriter( tMesh );
            tWriter.write_mesh( "./", aExoFileName, "./", "gen_temp.exo" );

            // write time to file
            tWriter.set_time( tTimeShift );

            // Set nodal fields based on field names
            moris::Cell< std::string > tNodalFieldNames = { "Mapped_Field", "Original_Field" };
            tWriter.set_nodal_fields( tNodalFieldNames );

            // Create field on mesh
            tWriter.write_nodal_field( tNodalFieldNames( 0 ), aTarget->get_values() );
            tWriter.write_nodal_field( tNodalFieldNames( 1 ), aSource->get_values() );

            // Finalize
            tWriter.close_file( true );
        }
    }    // namespace wrk
}    // namespace moris
