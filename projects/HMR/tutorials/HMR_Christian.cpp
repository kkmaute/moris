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
// Parameters, to be passed from XML input file
//------------------------------------------------------------------------------

    // determine path for object file
    std::string tHdf5FilePath = "AbsDesVariables0100.hdf5" ;

    // Name of Field we want to read from this file
    std::string tLabel = "AbsDesVariables";

    // the database for this HMR Object
    std::string tDatabasePath = "hmr_data.hdf5";

    // from this file, we will later read the fields for the refinement
    std::string tMeshPath     = "mbeam.e-s.0100";

    // path to the user defined shared object
    std::string tUserFun      = "/home/messe/build/userdef/hmr_userdef.so";

    // prepare parameters for user defined refinement Since the parameters
    // are very specific to the user defined function, they shold also
    // be generated in the shared object file
    ParameterList tParameters;
    tParameters.insert( "max_level",   ( sint ) 4 );
    tParameters.insert( "lower_bound", ( real ) 0.25 );

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
// Create the database
//------------------------------------------------------------------------------

      // we create a new HMR Object by passing the path to the Database
      HMR tHMR( tDatabasePath );

//------------------------------------------------------------------------------
// Load the field and perform a refinement
//------------------------------------------------------------------------------

      // - - - - - - - - - - - - - - - - - - - -
      // step 1: load coefficients from file
      // - - - - - - - - - - - - - - - - - - - -

      /** HMR Entertains the following refinement patterns
       *
       * 0 : input configuration   ( before refinement )
       * 1 : output configuration ( after refinement )
       * 2 : union configuration
       * 3 : special pattern if Lagrange refinement is independent from B-Spline
       * 4 : reserved for swapping
       */

      // at the moment, I only have one mesh per pattern, but I can create
      // linear, quadratic, and cubic meshes for each patterna at the same time.
      // the current function assumes that I am only interested in the mesh
      // with the maximum order. This needs to be changed
      uint tLagrangeOrder = tHMR.get_parameters()->get_lagrange_orders().max();

      // HMR has an internal parameter object which is created from the ParameterList.
      // But it also sets some internal things automatically and makes sure that the passed
      // parameter list makes sense. In the long run, the user is not meant to see that this
      // parameter object exists.

      // before I can crate a Field object, I need to create a mesh.
      // the moris::hmr::Mesh is a derivative of moris::mtk::Mesh, and can be also
      // used for any module that understands MTK, such as SDF or the Mapper

      // the following lines create a shared pointer to a Lagrange mesh that
      // is linked to the input configurations
      std::shared_ptr< Mesh > tInputMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_input_pattern() );


      // with the Mesh defined, I can now create a Field object.
      // I give the field a label, and look for the corresponding entry in the
      // HDF file. HMR assumes that the data refer to B-Spline values.

      // This function will be extended with an Entity Rank that defines if this
      // is a Nodal, Elemental or B-Spline field. ( see projects/MTK/cl_Mesh_Enums.hpp )
      // currently, I assume that I read in B-Spline coefficients.

      std::shared_ptr< Field > tField = std::make_shared< Field >(
              tLabel,
              tInputMesh,
              tHdf5FilePath );

      // In general, we might have an arbitrary number of fields. The API for the
      // user defined function requires a Cell to be passed
      Cell< std::shared_ptr< Field > > tFields;
      tFields.push_back( tField );

      // now we call the flagging. We can call this function an arbitrary number
      // of times, or maybe also run a function from the Geometry engine to flag more elements
      tHMR.user_defined_flagging( user_defined_refinement, tFields, tParameters );


      // After I have flagged all elements, I perform the refinement
      // By default, this function activates new elements on pattern 1.
      // I have to provide a nice API, so that I can switch this to pattern 3 if I only
      // want to perform Lagrange Elements. The call would look something like
      // tHMR.perform_refimenent( tHMR->get_lagrange_pattern() ) or so
      tHMR.perform_refinement();


      // There might be a case Where I want to do the flagging - and refienement routines multiple times,
      // for example if I want to call SDF in core. While doing that, there is some information I do not
      // need, such as the T-Matrices for the output mesh, or the faces and the Edges.
      //
      // Therefore, the finalize function is called only after the last call of perform_refinement
      tHMR.finalize();


      // - - - - - - - - - - - - - - - - - - - -
      // step 2: map coefficients on union mesh
      // - - - - - - - - - - - - - - - - - - - -

      // Now I have created the new mesh, and I could even dump it using tHMR->save_to_exodus( tPath ).
      // for this example, I call the mapper in core. The following lines need to be called
      // in an HMR specific mapper executable, into which I have to pass the database for the old and the new
      // refinement configuration.


      // first, I create a MTK Interface to the union mesh
      std::shared_ptr< Mesh > tUnionMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_union_pattern() );

      // Now I create a field object which will contain the unton field
      std::shared_ptr< Field > tUnionField = tUnionMesh->create_field(
              tLabel,
              tField->get_bspline_order() );

      // The following command takes the nodal data from the field of the Input mesh,
      // and projects it to the output field
      tHMR.get_database()->interpolate_field(
              tHMR.get_parameters()->get_input_pattern(),  // < -- since the field is linked to the pattern already, this parameter is redundant and will be removed
              tField,
              tHMR.get_parameters()->get_union_pattern(),  // < -- same goes for this one
              tUnionField );

      // - - - - - - - - - - - - - - - - - - - -
      // step 3: Create B-Spline coefficients on union mesh
      // - - - - - - - - - - - - - - - - - - - -

      // With the union mesh being created and populated with a field, I can now create an MTK mapper

      mapper::Mapper tUnionMapper( tUnionMesh );

      // Note that from within MTK, the field object does not exist.
      // Data is accesed via Label and Rank
      // map interpolated values to B-Spline Coefficients
      tUnionMapper.perform_mapping(
              tLabel,
              EntityRank::NODE,
              tLabel,
              tUnionField->get_bspline_rank() );


      // - - - - - - - - - - - - - - - - - - - -
      // step 4: Node values on output mesh
      // - - - - - - - - - - - - - - - - - - - -

      // With the B-Spline coefficients being created, we can now create an
      // output mesh and put a field on it
      std::shared_ptr< Mesh >  tOutputMesh = tHMR.create_mesh(
              tLagrangeOrder,
              tHMR.get_parameters()->get_output_pattern() );

      // The field is still called AbsDesVariables, and grabbs the B-Spline order
      // from the field we put on the input mesh
      std::shared_ptr< Field > tOutputField =
              tOutputMesh->create_field(
                      tLabel,
                      tField->get_bspline_order() );

      // Now I copy the coefficients from the union field.
      // Noah is currently working on a move operator for the moris::Matrix, so that
      // we can skip the copying, since we do not need the Coefficients in the union field
      // any more
      tOutputField->get_coefficients() = tUnionField->get_coefficients();

      // The following lines allocate the Matrix for the nodal values and evaluate
      // them according to the T-Matrices of the Mesh.
      // I wrote the evaluate_node_values function before I had a mapper.
      // To remove redundancy, this funciton will soon call the mtk::Mapper internally
      // instead of doing its own thing
      tOutputField->get_node_values().set_size( tOutputMesh->get_num_nodes(), 1 );
      tOutputField->evaluate_node_values();


      // - - - - - - - - - - - - - - - - - - - -
      // step 5: Filter
      // - - - - - - - - - - - - - - - - - - - -

      // this is some Fun I had on Thursday night. First, I create a new field object ...
      std::shared_ptr< Field >  tFilterField =  tOutputMesh->create_field(
              "FilteredValues",
              tField->get_bspline_order() );

      //  ... and create a mapper that is linked to the output mesh
      mapper::Mapper tOutputMapper( tOutputMesh );

      // Now I call the filter from the mapper object. Yesterday night, it made
      // sense to me that the filter is also a responsibility of the mapper.
      // I tried to write it it a way where it is independent of HMR.
      // I see the point of wanting to expoit the Tensor-like structure of a hierarchical
      // mesh. However, I am sceptical that this will be much faster, because the
      // "looking up" algorithm we talked about and implemented in the old HMR
      // contained many modulo operations, which are very expensive.
      //
      // This call just puts the filtered data on a moris::Matrix.
      // The filters stores the Node IDs and the weights internally
      // sence we want to build T-Matrices from this information later.
      //
      // If you call the filter for a second time without changing the radius,
      // it does not need to recalculate the weights
      //tOutputMapper.perform_filter( tLabel, 0.15, tFilterField->get_node_values() );

      // - - - - - - - - - - - - - - - - - - - -
      // step 6: Exodus Files
      // - - - - - - - - - - - - - - - - - - - -

      // Having four possible refinement patterns with three different orders,
      // HMR will create up to twelve different Lagrange Meshes.

      // calling
      // HMR.save_to_exodus( "Mesh.exo" );
      // will by default store  a first or second order mesh linked to pattern
      // 2.

      // By passing the field index, we can chose any of these meshes.
      // If there is only one Lagrange order present, the indices are by default : */

      // the configuration of the input mesh
      tHMR.save_to_exodus( 0, "LastStep.exo" );

      // the configuration of the output mesh
      //tHMR.save_to_exodus( 1, "Mesh.exo" );

      // the configuration of the Union Mesh
      //tHMR.save_to_exodus( 2, "UnionMesh.exo" );

      // - - - - - - - - - - - - - - - - - - - -
      // Apart from these meshes, you can also store debug files which
      // reveal internal information that is not meant to be shown in the Exodus.
      //
      //
      // If you want to try out the debugging files, you can do the Tutorial 3.
      //
      // - - - - - - - - - - - - - - - - - - - -

      // close user defined function
      dlclose( tFunctionHandle );

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
