//------------------------------------------------------------------------------
#include <memory>
#include <string>

// dynamik linker function
#include "dlfcn.h"

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from LINALG
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "HDF5_Tools.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"


//------------------------------------------------------------------------------
// from MTK
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_Mesh_Factory.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>


//------------------------------------------------------------------------------
// HMR
#define protected public
#define private   public
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#undef protected
#undef private
//------------------------------------------------------------------------------


namespace moris
{

    namespace hmr
    {
        /*
         * Interface for user defined function
         */
        typedef bool ( *MORIS_HMR_USER_FUNCTION )
                (
                        const Element                  * aElement,
                        const Cell< Matrix< DDRMat > > & aElementLocalValues,
                        ParameterList                  & aParameters
                );

    }
}

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
                    tHMR.perform_refinement();
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

                hid_t tFile = create_hdf5_file( tFilePath  );
                herr_t tStatus = 0;
                save_matrix_to_hdf5_file( tFile, "ActiveBasis", tActiveBasis, tStatus );

                if( tParameters.use_multigrid() )
                {
                    save_matrix_to_hdf5_file( tFile, "RefinedBasis", tRefinedBasis, tStatus );
                }
                close_hdf5_file( tFile );


          }
     }
  }
//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;
}
