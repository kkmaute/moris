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

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"


#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"


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

    /*bool
        user_defined_refinement(
                const Element                  * aElement,
                const Cell< Matrix< DDRMat > > & aElementLocalValues,
                      ParameterList            & aParameters )
        {
            // max level
            uint tMaxLevel   = aParameters.get< sint >("max_level");

            // minimal value
            real tLowerBound =  aParameters.get< real >("lower_bound");

            return  aElementLocalValues( 0 ).max() >= tLowerBound && aElement->get_level() < tMaxLevel;
        } */
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


    // determine path for object file
    std::string tHdf5FilePath = "AbsDesVariables0100.hdf5" ;
    std::string tDatabasePath = "hmr_data.hdf5";
    std::string tMeshPath     = "mbeam.e-s.0100";
    std::string tUserFun      = "/home/messe/build/userdef/hmr_userdef.so";

//------------------------------------------------------------------------------
// Open user defined function
//------------------------------------------------------------------------------

    // http://man7.org/linux/man-pages/man0/dlfcn.h.0p.html
    void * tFunctionHandle = dlopen( tUserFun.c_str(), RTLD_NOW );
    if ( !tFunctionHandle )
    {
        std::string tError = dlerror();
        std::cout << "Could not open the library" << std::endl;
        std::cout <<  tError << std::endl;
        return 1;
    }

    MORIS_HMR_USER_FUNCTION user_defined_refinement
        = reinterpret_cast<MORIS_HMR_USER_FUNCTION>
        ( dlsym( tFunctionHandle, "user_defined_refinement" ) );

    if (!user_defined_refinement) {
        std::cout << "Could not find symbol user_defined_refinement" << std::endl;
        dlclose( tFunctionHandle );
        return 1;
    }

//------------------------------------------------------------------------------

      HMR tHMR( tDatabasePath );

//------------------------------------------------------------------------------
// Mapping routine
//------------------------------------------------------------------------------
      std::string tLabel = "AbsDesVariables";
      uint tLagrangeOrder = tHMR.get_parameters()->get_lagrange_orders().max();

      // prepare parameters
      ParameterList tParameters;
      tParameters.insert( "max_level",   ( sint ) 4 );
      tParameters.insert( "lower_bound", ( real ) 0.25 );


      // create input mesh
      std::shared_ptr< Mesh > tInputMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_input_pattern() );

      // - - - - - - - - - - - - - - - - - - - -
      // step 1: load coefficients from file
      // - - - - - - - - - - - - - - - - - - - -

      std::shared_ptr< Field > tField = std::make_shared< Field >(
              tLabel,
              tInputMesh,
              tHdf5FilePath );


      // prapare cell with all fields
      Cell< std::shared_ptr< Field > > tFields;
      tFields.push_back( tField );

      tHMR.user_defined_flagging( user_defined_refinement, tFields, tParameters );

      // refine all flagged elements
      tHMR.perform_refinement();

      // finish mesh
      tHMR.finalize();


      // - - - - - - - - - - - - - - - - - - - -
      // step 2: map coefficients on union mesh
      // - - - - - - - - - - - - - - - - - - - -

      // create union
      auto tUnionMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_union_pattern() );

      // create union field object
      auto tUnionField = tUnionMesh->create_field(
              tLabel,
              tField->get_bspline_order() );

      // interpolate input field to union
      tHMR.get_database()->interpolate_field(
              tHMR.get_parameters()->get_input_pattern(),
              tField,
              tHMR.get_parameters()->get_union_pattern(),
              tUnionField );

      // - - - - - - - - - - - - - - - - - - - -
      // step 3: Create B-Spline coefficients on union mesh
      // - - - - - - - - - - - - - - - - - - - -

      // create mapper for union mesh
      mapper::Mapper tUnionMapper( tUnionMesh.get() );

      // map interpolated values to B-Spline Coefficients
      tUnionMapper.perform_mapping(
              tLabel,
              EntityRank::NODE,
              tLabel,
              tUnionField->get_bspline_rank() );


      // - - - - - - - - - - - - - - - - - - - -
      // step 4: Node values on output mesh
      // - - - - - - - - - - - - - - - - - - - -

      // get pointer to output mesh
      auto tOutputMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_output_pattern() );

      // create output field object
      std::shared_ptr< Field > tOutputField =
              tOutputMesh->create_field(
                      tLabel,
                      tField->get_bspline_order() );

      // allocate node values
      tOutputField->get_node_values().set_size( tOutputMesh->get_num_nodes(), 1 );

      // copy coefficients into output field ( @note may use move here in the future )
      tOutputField->get_coefficients() = tUnionField->get_coefficients();

      // evaluate field
      tOutputField->evaluate_node_values();

      auto tFilterField =  tOutputMesh->create_field(
              "FilteredValues",
              tField->get_bspline_order() );

      // create mapper for output
      mapper::Mapper tOutputMapper( tOutputMesh.get() );

      tOutputMapper.perform_filter( tLabel, 0.15, tFilterField->get_node_values() );

      // save exo files
      tHMR.save_to_exodus( 0, "LastStep.exo" );
      tHMR.save_to_exodus( 2, "UnionMesh.exo" );
      tHMR.save_to_exodus( 1, "Mesh.exo" );


      // close user defined function
      dlclose( tFunctionHandle );
//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
