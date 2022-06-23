#include "cl_WRK_Reinitialize_Performer.hpp"

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"
//#include "cl_MTK_Mesh_Pair.hpp"

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

namespace moris
{
    namespace wrk
    {
        // Parameter function
        typedef void ( *Parameter_Function )( moris::Cell< moris::Cell< moris::ParameterList > >& aParameterList );


        Reinitialize_Performer::Reinitialize_Performer( std::shared_ptr< Library_IO > aLibrary )
                : mLibrary( aLibrary )
        {

            // call the work flow parameter list and see if it is empty
            std::string        tWRKString            = "WRKParameterList";
            Parameter_Function tWRKParameterListFunc = aLibrary->load_function< Parameter_Function >( tWRKString, true );

            // extract the PRM list
            moris::Cell< moris::Cell< ParameterList > > tWRKParameterList;
            tWRKParameterListFunc( tWRKParameterList );

            // call the msi paramater list;
            std::string        tMSIString = "MSIParameterList";
            Parameter_Function tMSIParameterListFunc =
                    aLibrary->load_function< Parameter_Function >( tMSIString, true );

            // extract the PRM list
            moris::Cell< moris::Cell< ParameterList > > tMSIParameterList;
            tMSIParameterListFunc( tMSIParameterList );

            mAdofMeshIndex = tMSIParameterList( 0 )( 0 ).get< moris::sint >( tWRKParameterList( 0 )( 0 ).get< std::string >( "dof_type" ) );

            // get the adv field name that wll be reinitialized
            mADVFiledName = tWRKParameterList( 0 )( 0 ).get< std::string >( "adv_field" );

            // get msi string to dof type map
            moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                    moris::MSI::get_msi_dof_type_map();

            // get the quantity dof type from parameter list
            string_to_cell(
                    tWRKParameterList( 0 )( 0 ).get< std::string >( "dof_type" ),
                    mDofTypes,
                    tMSIDofTypeMap );


            mReinitializationFrequency = tWRKParameterList( 0 )( 0 ).get< sint >( "reinitialization_frequency" );

            // get the mesh output info
            mOutputMeshFile = tWRKParameterList( 0 )( 0 ).get< std::string >( "output_mesh_file" );
            mTimeOffset     = tWRKParameterList( 0 )( 0 ).get< real >( "time_offset" );
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
            Tracer tTracer( "WRK", "Reinitialize ADVs", "Perform reinitialize" );

            // initialize and populate the fields
            moris::Cell< std::shared_ptr< mtk::Field > > tGENFields;
            tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

            // find the index of the desired adv field that will be reinitialized
            auto itr = std::find_if( tGENFields.begin(), tGENFields.end(), [ & ]( std::shared_ptr< mtk::Field > const & aFiled )    //
                    { return aFiled->get_label() == mADVFiledName; } );

            // find the index of the adv field
            moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );

            // get the necessary data to build the target HMR mesh
            uint tTargetLagrangeOrder     = ( *itr )->get_lagrange_order();
            uint tDiscretizationOrder     = ( *itr )->get_discretization_order();
            uint tDiscretizationMeshIndex = ( *itr )->get_discretization_mesh_index();

            // get interpolation mesh from mesh pair
            moris::mtk::Mesh* tTargetMesh = ( *itr )->get_mesh_pair().get_interpolation_mesh();

            uint tTargetLagrangePattern = tTargetMesh->get_HMR_lagrange_mesh()->get_activation_pattern();
            uint tTargetBSplinePattern  = tTargetMesh->get_HMR_lagrange_mesh()->get_bspline_pattern( tDiscretizationMeshIndex );

            // build the target HMR mesh based on the adv field
            hmr::Interpolation_Mesh_HMR* tInterpolationMeshTarget = new hmr::Interpolation_Mesh_HMR(
                    aHMRPerformers( 0 )->get_database(),
                    tTargetLagrangeOrder,
                    tTargetLagrangePattern,
                    tDiscretizationOrder,
                    tTargetBSplinePattern );    // order, Lagrange pattern, bspline pattern

            // Create  mesh pair for the discrete field
            mtk::Mesh_Pair tMeshPairTarget( tInterpolationMeshTarget, nullptr, true );

            // get the source mesh
            moris::mtk::Mesh* tSourceMesh = aMTKPerformer( 0 )->get_mesh_pair( 0 ).get_interpolation_mesh();

            // get the source mesh bspline info
            uint tSourceBSplinePattern = tSourceMesh->get_HMR_lagrange_mesh()->get_bspline_pattern( mAdofMeshIndex );
            uint tSourceBSplineOrder   = tSourceMesh->get_HMR_lagrange_mesh()->get_bspline_order( mAdofMeshIndex );

            // build the hmr mesh which is a source
            hmr::Interpolation_Mesh_HMR* tInterpolationMeshSource = new hmr::Interpolation_Mesh_HMR(
                    aHMRPerformers( 0 )->get_database(),
                    tTargetLagrangeOrder,
                    tTargetLagrangePattern,
                    tSourceBSplineOrder,
                    tSourceBSplinePattern );    // order, Lagrange pattern, bspline pattern

            // Create  mesh pair
            mtk::Mesh_Pair tMeshPairSource( tInterpolationMeshSource, nullptr, true );

            // get the solution field and get a matrix of the solutions
            sol::Dist_Vector* tPartialSolutionVector = aMDLPerformer( 0 )->get_solver_interface()->get_solution_vector( mDofTypes );
            tPartialSolutionVector->extract_copy( mCoefficients );

            // delete the pointer as it is not needed anymore
            delete tPartialSolutionVector;

            // create field object for this mesh. with the source
            std::shared_ptr< mtk::Field_Discrete > tFieldSource = std::make_shared< mtk::Field_Discrete >( tMeshPairSource, 0 );

            // unlock fields and set the coeff
            tFieldSource->unlock_field();
            tFieldSource->set_coefficients( mCoefficients );

            // compute the nodal values based on the coeff
            tFieldSource->compute_nodal_values();

            // get the solution vector of the dof type that will be mapped to the adv field
            std::shared_ptr< mtk::Field_Discrete > tFieldTarget = std::make_shared< mtk::Field_Discrete >( tMeshPairTarget, 0 );
            tFieldTarget->set_label( mADVFiledName );

            //set the nodal values
            tFieldTarget->unlock_field();
            tFieldTarget->set_values( tFieldSource->get_values() );

            //invoke the mapper and map to the target field
            mtk::Mapper tMapper;
            tFieldTarget->unlock_field();
            tMapper.map_input_field_to_output_field_2( tFieldTarget.get() );

            //compute the nodal value
            tFieldTarget->compute_nodal_values();

            //get the coefficents and store them
            mCoefficients = tFieldTarget->get_coefficients();

            //clip the values
            this->impose_upper_lower_bound( aGENPerformer );

            //replace the newly constructed field
            tGENFields( tADVFieldIndex ) = tFieldTarget;

            //store the fields
            mMTKFields = tGENFields;

            //output the fields if asked
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
        Reinitialize_Performer::impose_upper_lower_bound( moris::Cell< std::shared_ptr< ge::Geometry_Engine > >& aGENPerformer )
        {
            Matrix< DDRMat > tLowerBounds = aGENPerformer( 0 )->get_lower_bounds();
            Matrix< DDRMat > tUpperBounds = aGENPerformer( 0 )->get_upper_bounds();

            // clip the values of the adv before putting in GEN
            for ( uint iADV = 0; iADV < mCoefficients.numel(); iADV++ )
            {
                mCoefficients( iADV ) = std::max( tLowerBounds( iADV ), std::min( mCoefficients( iADV ), tUpperBounds( iADV ) ) );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< mtk::Field > >
        Reinitialize_Performer::get_mtk_fields() const
        {
            return mMTKFields;
        }

        //------------------------------------------------------------------------------

        void
        Reinitialize_Performer::output_fields( mtk::Field* aTarget, mtk::Field* aSource , std::string aExoFileName) const
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
            tWriter.write_nodal_field( tNodalFieldNames(0), aTarget->get_values() );
            tWriter.write_nodal_field( tNodalFieldNames(1), aSource->get_values() );

            // Finalize
            tWriter.close_file( true );
        }
    }    // namespace wrk

}    // namespace moris
