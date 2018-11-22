//------------------------------------------------------------------------------
#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "HDF5_Tools.hpp"

#define protected public
#define private   public
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR.hpp" //HMR/src
#undef protected
#undef private

// select namespaces
using namespace moris;
using namespace hmr;


//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
//------------------------------------------------------------------------------

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );
//------------------------------------------------------------------------------

    std::string tMorisRoot = std::getenv("MORISROOT");

                // do this for first and second dimension
                for( uint tDimension=2; tDimension<=3; ++tDimension )
                {

                    // do this for first second and third order
                    for( uint tOrder=1; tOrder<3; ++tOrder )
                    {

                        for( uint tMultigrid = 0; tMultigrid<2; ++tMultigrid )
                        {

                            // The parameter object controls the behavior of HMR.
                            moris::hmr::Parameters tParameters;

                            tParameters.set_multigrid(  tMultigrid == 1  );

                            moris::Matrix< moris::DDLUMat > tNumberOfElements(  tDimension, 1,  2*tOrder );
                            tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

                            moris::Matrix< moris::DDRMat > tDomainOffset( tDimension, 1, 0.0 );
                            tParameters.set_domain_offset( tDomainOffset );

                            // set order of B-Splines
                            tParameters.set_mesh_orders_simple( tOrder );

                            // set buffer
                            tParameters.set_buffer_size( tOrder );

                            tParameters.set_verbose( false );

                            // create HMR Object
                            moris::hmr::HMR tHMR( tParameters );

                            // refine the mesh three times
                            for( uint k=0; k<3; ++k )
                            {
                                tHMR.flag_element( 0 );
                                tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
                            }

                            // finish mesh
                            tHMR.finalize();

                            // create matrix with IDs
                            moris::Matrix< moris::IdMat > tActiveBasis;
                            moris::Matrix< moris::IdMat > tRefinedBasis;

                            // reset counter
                            moris::uint tCount = 0;

                            moris::hmr::BSpline_Mesh_Base * tMesh = tHMR.mDatabase->mBSplineMeshes( 1 );

                            moris::Cell< moris::hmr::Basis *  > & mActiveBasisOnProc = tMesh->mActiveBasisOnProc;

                            tActiveBasis.set_size( mActiveBasisOnProc.size(), 1 );

                            // loop over all active basis
                            for( moris::hmr::Basis * tBasis :  mActiveBasisOnProc )
                            {
                                tActiveBasis( tCount++ ) = tBasis->get_id();
                            }

                            if( tParameters.use_multigrid() )
                            {
                                moris::Cell< moris::hmr::Basis *  > & mRefinedBasisOnProc = tMesh->mRefinedBasisOnProc;

                                tCount = 0;
                                tRefinedBasis.set_size( mRefinedBasisOnProc.size(), 1 );

                                for( moris::hmr::Basis * tBasis :  mRefinedBasisOnProc )
                                {
                                    tRefinedBasis( tCount++ ) = tBasis->get_id();
                                }
                            }

                            std::string tFilePath =  tMorisRoot + "/projects/HMR/test/data/hmr_bspline_ids_"
                                    + std::to_string( tDimension )
                            + std::to_string( tOrder )
                            + std::to_string( tMultigrid )
                            + ".hdf5";

                            /*
                            hid_t tFile = create_hdf5_file( tFilePath  );
                            herr_t tStatus = 0;
                            save_matrix_to_hdf5_file( tFile, "ActiveBasis", tActiveBasis, tStatus );

                            if( tParameters.use_multigrid() )
                            {
                                save_matrix_to_hdf5_file( tFile, "RefinedBasis", tRefinedBasis, tStatus );
                            }
                            close_hdf5_file( tFile ); */

                            moris::Matrix< moris::IdMat > tActiveBasisExpect;
                            moris::Matrix< moris::IdMat > tRefinedBasisExpect;

                            hid_t  tFile = moris::open_hdf5_file( tFilePath  );
                            herr_t tStatus = 0;
                            moris::load_matrix_from_hdf5_file( tFile,  "ActiveBasis", tActiveBasisExpect, tStatus );
                            //REQUIRE( moris::all_true( tActiveBasis == tActiveBasisExpect ) );
                            if( tParameters.use_multigrid() )
                            {
                                moris::load_matrix_from_hdf5_file( tFile,  "RefinedBasis", tRefinedBasisExpect, tStatus );
                                //REQUIRE( moris::all_true( tRefinedBasis == tRefinedBasisExpect ) );
                            }

                            moris::close_hdf5_file( tFile );

                        }
                    }
                }

//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;
}
