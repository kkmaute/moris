// standard library
#include <string>
#include <iostream>
#include <sstream>

// basic type definitions
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
#include "cl_Ascii.hpp"
#include "linalg_typedefs.hpp"

// Modules
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_HMR_Element.hpp"

// Standard Parameter Lists
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"

// SOL
#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Matrix Tools
#include "fn_equal_to.hpp"
#include "fn_find.hpp"

// Other Tools
#include "fn_stringify.hpp"
#include "fn_stringify_matrix.hpp"


//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    /* ------------------------------------------------------------------------ */
    // basic setup parameters

    // std::string gGeometryFilePath = "/home/wunsch/work/PotMs/nTop_Sandwich/Geometry/"; // BLUEFLAME
    // std::string gGeometryFilePath = "/home/wunsch/work/XTK_Paper/scaling_study/geometry/"; // MIRAGE
    std::string gGeometryFilePath = "/scratch/alpine/niwu6047/work/XTK_Paper/nTop_Sandwich/geometry/"; // ALPINE

    // file name of the input file and the output files
    std::string gFileName = "nTop_Sandwich";

    // solver setup
    uint gSolverType = 2;
    // 0 - Amesos-Pardiso (Direct Solver)
    // 1 - Amesos-Mumps (Direct Solver)
    // 2 - BELOS Flexible GMRES (Iterative Solver) + Non-Symmetric Smoothed Aggregation, NSSA (Pre-Cond)
    // 3 - AZTEC GMRES (Iterative Solver) + Non-Symmetric Smoothed Aggregation, NSSA (Pre-Cond)
    
    // scale the problem size by 1/gScale in X and Y directions (use {1,2,4,8} to retain pattern of lattices)
    uint gScaleX = i-0000;
    uint gScaleY = i-0001; // current setup has 2x2 unit cells
    // std::string gProcDomainSplit = "1,1,1"; // = 1 proc
    std::string gProcDomainSplit = "i-0002";
    
    // refine the base mesh
    real gRefinementFactor = 1.0;
    
    // interpolation order for both the Lagrange and B-spline meshes
    std::string tPolyOrder = "2";

    /* ------------------------------------------------------------------------ */
    // image parameters

    // dimension of image domain
    real gImgDimX = 7.650000e-02;    // x-direction
    real gImgDimY = 7.650000e-02;    // y-direction
    real gImgDimZ = 4.050000e-02;    // z-direction

    // offset of image domain
    real gImgOffX = -3.815000e-02;    // x-direction
    real gImgOffY = -3.815000e-02;    // y-direction
    real gImgOffZ = -2.015000e-02;    // z-direction

    // image voxel dimension
    real gDeltaX = 3.000000e-04;
    real gDeltaY = 3.000000e-04;
    real gDeltaZ = 3.000000e-04;

    // number of voxels in image
    uint gNumVoxelX = 255;
    uint gNumVoxelY = 255;
    uint gNumVoxelZ = 135;

    // approximate element size
    real gApproxVoxelSize = gImgDimX / 200.0;

    // refinement parameters
    real             gMaxLevelsetValue = gApproxVoxelSize / 1.5;    // 4.384164e-03;
    Matrix< DDRMat > gRefineRegions    = { { gMaxLevelsetValue, 0.5 * gMaxLevelsetValue, 0.25 * gMaxLevelsetValue, 0.125 * gMaxLevelsetValue } };

    // TODO ?
    int gRefinementLimit = 0;

    /* ------------------------------------------------------------------------ */
    // background mesh parameters

    // offset domain by some number of voxels into the interior of the image
    uint gNumDomainElemOffX = 6;
    uint gNumDomainElemOffY = 6;
    uint gNumDomainElemOffZ = 4;
    
    // get the mesh dimensions
    uint gNumElemsX = std::ceil( gImgDimX / gApproxVoxelSize ) - 2 * gNumDomainElemOffX;
    uint gNumElemsY = std::ceil( gImgDimY / gApproxVoxelSize ) - 2 * gNumDomainElemOffY;
    uint gNumElemsZ = std::ceil( gImgDimZ / gApproxVoxelSize ) - 2 * gNumDomainElemOffZ;
    
    // get the actual element size of the non-refined mesh
    real gElemSizeX = gImgDimX / ( (real)gNumElemsX + 2.0 * (real)gNumDomainElemOffX );
    real gElemSizeY = gImgDimY / ( (real)gNumElemsY + 2.0 * (real)gNumDomainElemOffY );
    real gElemSizeZ = gImgDimZ / ( (real)gNumElemsZ + 2.0 * (real)gNumDomainElemOffZ );
    
    // get the background mesh offset
    real gDomainOffX = gImgOffX + (real)gNumDomainElemOffX * gElemSizeX;
    real gDomainOffY = gImgOffY + (real)gNumDomainElemOffY * gElemSizeY;
    real gDomainOffZ = gImgOffZ + (real)gNumDomainElemOffZ * gElemSizeZ;
    
    // domain size
    real gDomainLengthX = gNumElemsX * gElemSizeX;
    real gDomainLengthY = gNumElemsY * gElemSizeY;
    real gDomainLengthZ = gNumElemsZ * gElemSizeZ;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tFiberPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tFillerPhase  = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tAllBulkSets  = tFiberPhase + "," + tFillerPhase;
    
    std::string tTopSide = "SideSet_6_n_p0,SideSet_6_c_p0";
    std::string tBottomSide = "SideSet_5_n_p0,SideSet_5_c_p0";
    std::string tInterfaceFiberSide = "iside_b0_0_b1_1";
    std::string tInterfaceFillerSide = "iside_b0_1_b1_0";
    std::string tInterface = "dbl_iside_p0_0_p1_1";
    std::string tAllInterfaces = tTopSide + "," + tBottomSide + "," + tInterfaceFiberSide + "," + tInterfaceFillerSide + "," + tInterface;

    /* ------------------------------------------------------------------------ */

    int
    Func_Refinement(
            hmr::Element*           aElement,
            const Matrix< DDRMat >& aElementLocalValues )
    {
        // current refinement level of element
        uint tCurrentLevel = aElement->get_level();

        // get a representative distance value of the element
        real tMinDistVal = std::min( std::abs( aElementLocalValues.max() ), std::abs( aElementLocalValues.min() ) );

        // default: do not refine
        int tRefine = 0;

        // if element is still within the refinement region for the next higher level, refine it
        if ( (int)tCurrentLevel < gRefinementLimit )
        {
            if ( tMinDistVal < gRefineRegions( tCurrentLevel ) )
            {
                tRefine = 1;
            }
        }

        // return refinement instruction
        return tRefine;
    }

    /* ------------------------------------------------------------------------ */

    // Phase assignment
    uint
    get_phase_index( const Bitset< 1 >& aGeometrySigns )
    {
        // In Filler
        if ( aGeometrySigns.test( 0 ) )
        {
            return 1;    // Filler
        }

        return 0;    // Fiber
    }

    /* ------------------------------------------------------------------------ */
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 3 );

        tParameterList( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        tParameterList( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).push_back( prm::create_hmr_parameter_list() );

        // define how the domain gets split up across the processors
        tParameterList( 0 )( 0 ).set( "processor_decomposition_method", (int)0 );
        tParameterList( 0 )( 0 ).set( "processor_dimensions", gProcDomainSplit );

        // convert to strings
        std::string tNumElemsPerDim = 
                ios::stringify( (uint)( gRefinementFactor * gNumElemsX / gScaleX ) ) + "," + 
                ios::stringify( (uint)( gRefinementFactor * gNumElemsY / gScaleY ) )+ "," + 
                ios::stringify( (uint)( gRefinementFactor * gNumElemsZ ) );
        std::string tDomainDims = 
                ios::stringify( gDomainLengthX / (real)gScaleX ) + "," + 
                ios::stringify( gDomainLengthY / (real)gScaleY ) + "," + 
                ios::stringify( gDomainLengthZ );
        std::string tDomainOffset = 
                ios::stringify( gDomainOffX ) + "," + 
                ios::stringify( gDomainOffY ) + "," + 
                ios::stringify( gDomainOffZ );
    
        // define base grid
        tParameterList( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterList( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterList( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterList( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4,5,6" );

        // define patterns to be created and how they are refined - all geometry is based off pattern 0
        tParameterList( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
        tParameterList( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterList( 0 )( 0 ).set( "refinement_function_names", "Func_Refinement" );
        tParameterList( 0 )( 0 ).set( "refinement_buffer", 0 );
        tParameterList( 0 )( 0 ).set( "staircase_buffer", 0 );

        // define basis functions to be constructed on grids - the lagrange mesh should use pattern 0 which all geometry is based off
        tParameterList( 0 )( 0 ).set( "lagrange_orders", tPolyOrder );
        tParameterList( 0 )( 0 ).set( "lagrange_pattern", "0" );
        tParameterList( 0 )( 0 ).set( "bspline_orders", tPolyOrder );
        tParameterList( 0 )( 0 ).set( "bspline_pattern", "0" );
        tParameterList( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
        tParameterList( 0 )( 0 ).set( "truncate_bsplines", 1 );

        // misc parameters
        tParameterList( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        // advanced refinement options
        tParameterList( 0 )( 0 ).set( "use_advanced_T_matrix_scheme", 1 );
        tParameterList( 0 )( 0 ).set( "use_refine_low_level_elements", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );

        tParameterList( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterList( 0 )( 0 ).set( "decompose", true );
        tParameterList( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterList( 0 )( 0 ).set( "basis_rank", "bspline" );              // use option "node" for Lagrange discretization
        tParameterList( 0 )( 0 ).set( "enrich", true );
        tParameterList( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterList( 0 )( 0 ).set( "ghost_stab", true );                   // generate ghost side sets
        tParameterList( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );    // generate double side sets both ways

        // output
        tParameterList( 0 )( 0 ).set( "verbose", true );                      // turn on/off verbose console output
        tParameterList( 0 )( 0 ).set( "verbose_level", moris::uint( 1 ) );
        tParameterList( 0 )( 0 ).set( "print_enriched_ig_mesh", false );      // output summary of enr. IG mesh to console
        tParameterList( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );    // create xtk_temp.exo
        tParameterList( 0 )( 0 ).set( "output_cut_ig_mesh", false );          // create cut_ig_mesh.exo
        tParameterList( 0 )( 0 ).set( "output_path", "./XTK/" );              // path to write xtk_temp.exo into

        tParameterList( 0 )( 0 ).set( "only_generate_xtk_temp", true );       // kill workflow after XTK        
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 3 );
        tParameterList( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterList( 0 )( 0 ).set( "number_of_phases", 2 );
        tParameterList( 0 )( 0 ).set( "phase_function_name", "get_phase_index" );
        // tParameterList( 0 )( 0 ).set( "output_mesh_file", "./GEN/" + gFileName + "_GEN.exo" );
        tParameterList( 0 )( 0 ).set( "time_offset", 10.0 );

        //------------------------------------------------------------------------------
        // init geometry counter
        uint tGeoCounter = 0;

        // 3D Image
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SIGNED_DISTANCE_IMAGE ) );
        tParameterList( 1 )( tGeoCounter ).set( "name", "NODLEVELSET" );

        // image parameters
        Vector< real > tImgDims = { gImgDimX, gImgDimY, gImgDimZ };
        Vector< real > tImgOffSet = { gImgOffX, gImgOffY, gImgOffZ };
        tParameterList( 1 )( tGeoCounter ).set( "image_file", gGeometryFilePath + "NTop.hdf5" );
        tParameterList( 1 )( tGeoCounter ).set( "image_dimensions", tImgDims );
        tParameterList( 1 )( tGeoCounter ).set( "image_offset", tImgOffSet );
        tParameterList( 1 )( tGeoCounter ).set( "image_sdf_shift", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "image_sdf_default", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "image_sdf_interpolate", true );
        tParameterList( 1 )( tGeoCounter ).set( "image_sdf_scaling", 1.0 );

        // mesh cutting
        tParameterList( 1 )( tGeoCounter ).set( "isocontour_threshold", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1.0e-12 );
        tParameterList( 1 )( tGeoCounter ).set( "intersection_tolerance", 1.0e-12 );

        // meshing parameters
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", true );
        tParameterList( 1 )( tGeoCounter ).set( "number_of_refinements", Vector< uint >({0}) );
        tParameterList( 1 )( tGeoCounter ).set( "refinement_mesh_index", Vector< uint >({0}) );
        // tParameterList( 1 )( tGeoCounter ).set( "refinement_function_index",  0 );

        // ADV parameters
        tParameterList( 1 )( tGeoCounter ).set( "discretization_mesh_index", -1 );
        tParameterList( 1 )( tGeoCounter ).set( "discretization_lower_bound", -100.0 * gApproxVoxelSize );
        tParameterList( 1 )( tGeoCounter ).set( "discretization_upper_bound", 100.0 * gApproxVoxelSize );
        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 9 );
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------
        // phase info
        uint tPhaseCounter = 0;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseFiber" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseFiller" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseAll" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0,1" );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        //---------------------
        // Structural Material Properties
        
        // Young's modulus, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropYoungsModulusFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 2.0e11 ) );    // ~200 GPa for Carbon fibers
        tPropCounter++;

        // Poisson ratio, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPoissonRatioFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 0.27 ) );    // ~0.27 for Carbon fibers
        tPropCounter++;
        
        // coefficient of thermal expansion, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCteFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( -1.0e-6 ) );    // ~ -1e-6 for Carbon fibers typically found in AS
        tPropCounter++;

        // Young's modulus, filler material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropYoungsModulusFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 3.1e9 ) );    // 2.6~3.6 GPa for cured epoxy resin
        tPropCounter++;

        // Poisson ratio, filler material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPoissonRatioFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 0.32 ) );    // typically 0.29~0.34 for cured epoxy resin
        tPropCounter++;
        
        // coefficient of thermal expansion, filler material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCteFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 7.0e-5 ) );    // 60~80e-6 for cured epoxy (cheap non-AS stuff)
        tPropCounter++;
        
        //---------------------
        // Thermal Material Properties
        
        // thermal conductivity, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivityFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 900.0 ) ); // 900~1000 W/m/K for Carbon fibers
        tPropCounter++;
        
        // heat capacity, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropHeatCapacityFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 1.1e3 ) );    // 1000~1200 J/kg/K for Carbon fibers
        tPropCounter++;
        
        // mass density, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensityFiber" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 1.8e3 ) );    // 1750~2000 kg/m^3 for Carbon fibers
        tPropCounter++;
        
        // thermal conductivity, filler material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivityFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 0.2 ) );    // 0.15~0.25 W/m/K for cured epoxy resin
        tPropCounter++;
        
        // heat capacity, filler material
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropHeatCapacityFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 1.1e3 ) );    // 900~1200 J/kg/K for cured epoxy resin
        tPropCounter++;
        
        // mass density, fiber
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensityFiller" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", ios::stringify( 1.1e3 ) );    // 1000~1200 kg/m^3 for cured epoxy resin
        tPropCounter++;

        //---------------------
        // Structural BC/Loading Properties
        
        // prescribed displacement bottom plate
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPrescDisplBottom" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" ); // clamped
        tPropCounter++;

        // prescribed loading on top plate
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropNeumann" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0;1.0e6" ); // 1 MPa in Z-direction
        tPropCounter++;
        
        //---------------------
        // Thermal BC/Loading Properties
        
        // prescribed temperature bottom plate
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPrescTempBottom" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "273.0" ); // 0°C / 273K at the bottom
        tPropCounter++;

        // prescribed temperature top plate
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPrescTempTop" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "373.0" ); // 100°C / 373K at the top
        tPropCounter++;
        
        // reference temperature in non-loaded state
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropRefTemp" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "293.0" ); // 20°C / 293K room temperature as reference point
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // structural CM for fibers
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CM_Struct_Fiber" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFiber" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "model_type", (uint)fem::Model_Type::PLANE_STRESS );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ;TEMP", "Displacement,Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", 
                "PropYoungsModulusFiber,YoungsModulus;" 
                "PropPoissonRatioFiber,PoissonRatio;"
                "PropCteFiber,CTE;"
                "PropRefTemp,ReferenceTemperature" );
        tCMCounter++;

        // structural CM for filler material
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CM_Struct_Filler" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFiller" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "model_type", (uint)fem::Model_Type::PLANE_STRESS );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ;TEMP", "Displacement,Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", 
                "PropYoungsModulusFiller,YoungsModulus;"
                "PropPoissonRatioFiller,PoissonRatio;"
                "PropCteFiller,CTE;"
                "PropRefTemp,ReferenceTemperature" );
        tCMCounter++;
        
        // thermal CM for fibers
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CM_Thermal_Fiber" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFiber" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", 
                "PropConductivityFiber,Conductivity;"
                "PropDensityFiber,Density;"
                "PropHeatCapacityFiber,HeatCapacity" );
        tCMCounter++;        
        
        // thermal CM for filler material
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CM_Thermal_Filler" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFiller" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", (uint)fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties", 
                "PropConductivityFiller,Conductivity;"
                "PropDensityFiller,Density;"
                "PropHeatCapacityFiller,HeatCapacity" );
        tCMCounter++;   

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        uint tSPCounter = 0;

        //---------------------
        // Nitsche penalty parameters for BCs/ICs

        // displacement BC on fiber
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Displ_BC_Fiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "50.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropYoungsModulusFiber,Material" );
        tSPCounter++;
        
        // tied structural contact between fiber and filler material
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Displ_Interface" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "50.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropYoungsModulusFiber,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", "PropYoungsModulusFiller,Material" );
        tSPCounter++;
        
        // temerature BC on Fiber
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Temp_BC_Fiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "50.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivityFiber,Material" );
        tSPCounter++;
        
        // thermal contact between fiber and filler material
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Temp_Interface" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::NITSCHE_INTERFACE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "50.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivityFiber,Material" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_properties", "PropConductivityFiller,Material" );
        tSPCounter++;

        //---------------------
        // Ghost penalty parameters for stabilization
        
        // displacement field, fiber
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Ghost_Displ_Fiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.001" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropYoungsModulusFiber,Material" );
        tSPCounter++;
        
        // displacement field, filler material
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Ghost_Displ_Filler" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.001" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropYoungsModulusFiller,Material" );
        tSPCounter++;

        // temperature field, fiber
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Ghost_Temp_Fiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiber" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.001" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivityFiber,Material" );
        tSPCounter++;
        
        // temperature field, filler material
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SP_Ghost_Temp_Filler" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "0.001" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivityFiller,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        //---------------------
        // bulk IWGs (structure & thermal)

        // structure, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Struct_Fiber_Bulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,ElastLinIso" );
        tIWGCounter++;
        
        // structure, filler material
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Struct_Filler_Bulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Struct_Filler,ElastLinIso" );
        tIWGCounter++;
        
        // thermal, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Thermal_Fiber_Bulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Thermal_Fiber,Diffusion" );
        tIWGCounter++;
        
        // thermal, filler material
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Thermal_Filler_Bulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Thermal_Filler,Diffusion" );
        tIWGCounter++;

        //---------------------
        // structural BCs/ICs

        // displacement Dirichlet, fiber -> clamped fiber base plate
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Displ_Dirichlet_BC_Bottom_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "5" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescDisplBottom,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,ElastLinIso" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Displ_BC_Fiber,DirichletNitsche" );
        tIWGCounter++;

        // displacement Neumann, fiber -> forcing on top plate
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Struct_Neumnan_BC_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "6" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropNeumann,Traction" );
        tIWGCounter++;
        
        // displacement interface, fiber-filler -> tied contact between the two materials
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Struct_Interface_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,ElastLinIso" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_constitutive_models", "CM_Struct_Filler,ElastLinIso" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Displ_Interface,NitscheInterface" );
        tIWGCounter++;

        //---------------------
        // thermal BCs/ICs
        
        // temperature Dirichlet TOP, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Temp_Dirichlet_BC_Top_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "6" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescTempTop,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Thermal_Fiber,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Temp_BC_Fiber,DirichletNitsche" );
        tIWGCounter++;
        
        // temperature Dirichlet BOTTOM, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Temp_Dirichlet_BC_Bottom_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "5" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescTempBottom,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Thermal_Fiber,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Temp_BC_Fiber,DirichletNitsche" );
        tIWGCounter++;
        
        // temperature interface, fiber-filler
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Thermal_Interface_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CM_Thermal_Fiber,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_constitutive_models", "CM_Thermal_Filler,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Temp_Interface,NitscheInterface" );
        tIWGCounter++;

        //---------------------
        // Ghost IWGs

        // displacement, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Ghost_Displ_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Ghost_Displ_Fiber,GhostSP" );
        tIWGCounter++;
        
        // displacement, filler material
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Ghost_Displ_Filler" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Ghost_Displ_Filler,GhostSP" );
        tIWGCounter++;
        
        // temperature, fiber
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Ghost_Temp_Fiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiber" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Ghost_Temp_Fiber,GhostSP" );
        tIWGCounter++;
        
        // temperature, filler material
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWG_Ghost_Temp_Filler" ); 
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SP_Ghost_Temp_Filler,GhostSP" );
        tIWGCounter++;
        

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        //---------------------
        // structural quantities

        // x displacement
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Disp_X" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // y displacement
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Disp_Y" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tIQICounter++;

        // z displacement
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Disp_Z" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 2 );
        tIQICounter++;
        
        // strain energy - Fiber
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Strain_Energy_Fiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::STRAIN_ENERGY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,Elast" );
        tIQICounter++;
        
        // strain energy - Filler
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Strain_Energy_Filler" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::STRAIN_ENERGY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Filler,Elast" );
        tIQICounter++;

        // von Mises stress - Fiber
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Stress_VM_Fiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::VON_MISES_STRESS );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,ElastLinIso" );
        tIQICounter++;
        
        // von Mises stress - Filler material
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Stress_VM_Filler" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::VON_MISES_STRESS );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Filler,ElastLinIso" );
        tIQICounter++;
        
        // traction on interface for fiber - first component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Fiber_X" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;
        
        // traction on interface for fiber - second component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Fiber_Y" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tIQICounter++;
        
        // traction on interface for fiber - third component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Fiber_Z" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 2 );
        tIQICounter++;
        
        // traction jump on the interface - first component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Jump_X" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::JUMP_TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_constitutive_models", "CM_Struct_Filler,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;
        
        // traction jump on the interface - first component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Jump_Y" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::JUMP_TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_constitutive_models", "CM_Struct_Filler,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tIQICounter++;
        
        // traction jump on the interface - first component
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Traction_Jump_Z" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::JUMP_TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_phase_name", "PhaseFiller" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_dof_dependencies", "UX,UY,UZ;TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Struct_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "follower_constitutive_models", "CM_Struct_Filler,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 2 );
        tIQICounter++;
        
        //---------------------
        // thermal quantities
        
        // temperature
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Temp" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;
        
        // top surface heat flux
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQI_Surface_Heat_Flux" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", (uint)fem::IQI_Type::TRACTION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", (uint)fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFiber" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "side_ordinals", "5,6" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CM_Thermal_Fiber,TractionCM" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        // initialize the SOL parameter lists
        tParameterList.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterList( Ik ).resize( 1 );
        }

        switch ( gSolverType )
        {
            // Amesos-Pardiso Direct Solver
            case 0:
            {
                tParameterList( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
                tParameterList( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
                break;
            }

            // Amesos-Mumps Direct Solver
            case 1:
            {
                tParameterList( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
                tParameterList( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
                tParameterList( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
                break;
            }

            // BELOS Flexible GMRES (Iterative Solver) + NSSA (Pre-Cond)
            case 2:
            {
                tParameterList( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
                tParameterList( 0 )( 0 ).set( "preconditioners", "0" );

                // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES,
                //              Recycling CG, MINRES, LSQR, TFQMR, Pseudoblock TFQMR, Seed GMRES, Seed CG
                tParameterList( 0 )( 0 ).set( "Solver Type", "Flexible GMRES" );

                // Diagnostics
                sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary;
                tParameterList( 0 )( 0 ).set( "Verbosity", tVerbosity );

                // Maximum number of blocks in Krylov factorization
                tParameterList( 0 )( 0 ).set( "Num Blocks", 500 );

                // Block size to be used by iterative solver
                tParameterList( 0 )( 0 ).set( "Block Size", 1 );

                // Allowable Belos solver iterations
                tParameterList( 0 )( 0 ).set( "Maximum Iterations", 1000 );

                // Allowable Belos solver iterations
                tParameterList( 0 )( 0 ).set( "Maximum Restarts", 2 );

                // Convergence criteria
                tParameterList( 0 )( 0 ).set( "Convergence Tolerance", 1e-9 );

                // Left or right preconditioner
                tParameterList( 0 )( 0 ).set( "Left-right Preconditioner", "right" );

                // Preconditioner

                // ifpack - ILU
                // tParameterList( 0 )( 0 ).set( "ifpack_prec_type",  "ILU");
                // tParameterList( 0 )( 0 ).set( "fact: level-of-fill",  5 );

                // ifpack - ILUT
                // tParameterList( 0 )( 0 ).set( "ifpack_prec_type",  "ILUT");
                // tParameterList( 0 )( 0 ).set( "fact: ilut level-of-fill", 5.0 );
                // tParameterList( 0 )( 0 ).set( "fact: drop tolerance", 1e-12 );

                // ifpack with direct solve
                // tParameterList( 0 )( 0 ).set( "ifpack_prec_type",  "amesos");
                // tParameterList( 0 )( 0 ).set( "amesos: solver type", "Amesos_Pardiso");

                // AMG with defaults for non-symmetric system
                tParameterList( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML );
                tParameterList( 7 )( 0 ).set( "ml_prec_type", "NSSA" );
                tParameterList( 7 )( 0 ).set( "PDE equations", 3 );

                break;

            }    // end solver case: Belos Flexible GMRES (Iterative Solver) + NSSA (Pre-Cond)

            // AZTEC GMRES (Iterative Solver) + NSSA (Pre-Cond)
            case 3:
            {
                tParameterList( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
                tParameterList( 0 )( 0 ).set( "preconditioners", "0" );

                // options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
                tParameterList( 0 )( 0 ).set( "AZ_solver", AZ_gmres );

                // Allowable Aztec solver iterations
                tParameterList( 0 )( 0 ).set( "AZ_max_iter", 500 );

                // Allowable Aztec iterative residual
                tParameterList( 0 )( 0 ).set( "rel_residual", 1e-08 );

                // set Az_conv -convergence criteria
                // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
                tParameterList( 0 )( 0 ).set( "AZ_conv", AZ_r0 );

                // set Az_diagnostic parameters
                // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
                tParameterList( 0 )( 0 ).set( "AZ_diagnostics", AZ_all );

                // set AZ_output options
                // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
                tParameterList( 0 )( 0 ).set( "AZ_output", AZ_all );

                // Determines the sub-matrices factored with the domain decomposition algorithms
                // Option to specify with how many rows from other processors each processor\u2019s local sub-matrix is augmented.
                tParameterList( 0 )( 0 ).set( "AZ_overlap", 1 );

                // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
                // Options are AZ_standard, AZ_symmetric
                tParameterList( 0 )( 0 ).set( "AZ_type_overlap", AZ_standard );

                // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
                // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
                tParameterList( 0 )( 0 ).set( "AZ_reorder", 1 );

                // Use preconditioner from a previous Iterate() call
                // Option are AZ_calc, AZ_recalc, AZ_reuse
                tParameterList( 0 )( 0 ).set( "AZ_pre_calc", AZ_calc );

                // Determines  whether  matrix  factorization  information will be kept after this solve
                // for example for preconditioner_recalculation
                tParameterList( 0 )( 0 ).set( "AZ_keep_info", 0 );

                //-------------------------- GMRES specific solver parameters  ------------------------
                // Set AZ_kspace
                // Krylov subspace size for restarted GMRES
                // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
                // For very difficult problems, set it equal to the maximum number of iterations.
                tParameterList( 0 )( 0 ).set( "AZ_kspace", 250 );

                // Set AZ_orthog
                // AZ_classic or AZ_modified
                tParameterList( 0 )( 0 ).set( "AZ_orthog", AZ_classic );

                // Set AZ_rthresh
                // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                tParameterList( 0 )( 0 ).set( "AZ_rthresh", 0.0 );

                // Set AZ_athresh
                // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
                // any of the incomplete factorization preconditioners
                tParameterList( 0 )( 0 ).set( "AZ_athresh", 0.0 );

                //-------------------------- Pre-conditioner specific parameters ---------------------------
                // Determine which preconditioner is used
                // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
                tParameterList( 0 )( 0 ).set( "AZ_precond", AZ_dom_decomp );

                // Set preconditioner subdomain solve - direct solve or incomplete
                // Options are AZ_lu, AZ_ilut, AZ_ilu, AZ_rilu, AZ_bilu, AZ_icc
                tParameterList( 0 )( 0 ).set( "AZ_subdomain_solve", AZ_ilut );

                // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
                tParameterList( 0 )( 0 ).set( "AZ_poly_ord", 3 );

                // Set drop tolerance - for LU, ILUT
                tParameterList( 0 )( 0 ).set( "AZ_drop", 1.0e-12 );

                // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
                tParameterList( 0 )( 0 ).set( "AZ_graph_fill", 3 );

                // Set ilut fill
                tParameterList( 0 )( 0 ).set( "AZ_ilut_fill", 15.0 );

                // Set Damping or relaxation parameter used for RILU
                tParameterList( 0 )( 0 ).set( "AZ_omega", 1.0 );

                // Pre-conditioner using ifpack
                // tParameterList( 0 )( 0 ).set( "ifpack_prec_type"    , "ILU" );
                // tParameterList( 0 )( 0 ).set( "fact: level-of-fill" ,  3 );
                // tParameterList( 0 )( 0 ).set( "fact: drop tolerance",  1.0e-2 );
                // tParameterList( 0 )( 0 ).set( "prec_reuse"          ,  false );

                // AMG with defaults for non-symmetric system                
                tParameterList( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML );
                tParameterList( 7 )( 0 ).set( "ml_prec_type", "NSSA" );
                tParameterList( 7 )( 0 ).set( "PDE equations", 3 );

                break;

            }    // end solver case: AZTEC GMRES (Iterative Solver) + NSSA (Pre-Cond)

        }        // end switch: solver cases

        //------------------------------------------------------------------------------

        tParameterList( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterList( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------
        // Non-linear solver algorithms
        tParameterList( 2 ).resize( 2 );

        // individual solves
        tParameterList( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterList( 2 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterList( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-8 );
        tParameterList( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterList( 2 )( 0 ).set( "NLA_max_iter", 1 );
        
        // overarching block solve
        tParameterList( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterList( 2 )( 1 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterList( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 1e-8 );
        tParameterList( 2 )( 1 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------
        // Non-lienar solvers
        tParameterList( 3 ).resize( 3 );

        // temperature solve
        tParameterList( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterList( 3 )( 0 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" ); // use NEWTON algorithm
        tParameterList( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );
        
        // displacement solve
        tParameterList( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterList( 3 )( 1 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "0" ); // use NEWTON algorithm
        tParameterList( 3 )( 1 ).set( "NLA_DofTypes", "UX,UY,UZ" );
        
        // overarching block solve
        tParameterList( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterList( 3 )( 2 ).set( "NLA_Solver_Implementation", (uint)moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterList( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "1" ); // use NLBGS algorithm
        tParameterList( 3 )( 2 ).set( "NLA_Sub_Nonlinear_Solver", "0,1" );         // set sub nonlinear solvers with index 0 and 1
        tParameterList( 3 )( 2 ).set( "NLA_DofTypes", "TEMP;UX,UY,UZ" );

        // ----------------------------------------------------------

        tParameterList( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterList( 4 )( 0 ).set( "TSA_Nonlinear_Solver", 2 ); // use NLBGS that again uses the individual NEWTON solves

        //------------------------------------------------------------------------------

        tParameterList( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterList( 5 )( 0 ).set( "TSA_DofTypes", "TEMP;UX,UY,UZ" );

        //------------------------------------------------------------------------------

        tParameterList( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        // tParameterList( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "Matrix" );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );
        tParameterList( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterList( 0 )( 0 ).set( "order_adofs_by_host", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );

        tParameterList( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterList( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "VIS/" + gFileName + ".exo" ) );
        tParameterList( 0 )( 0 ).set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );
        tParameterList( 0 )( 0 ).set( "Set_Names", tAllBulkSets + "," + tAllInterfaces );
        tParameterList( 0 )( 0 ).set( "Field_Names", 
                "UX,UY,UZ,TEMP,"
                "STRESS_VM_Fiber,STRESS_VM_Filler,"
                "STRAIN_ENERGY_FIBER,STRAIN_ENERGY_FILLER,"
                "TRACTION_FIBER_X,TRACTION_FIBER_Y,TRACTION_FIBER_Z,"
                "TRACTION_JUMP_X,TRACTION_JUMP_Y,TRACTION_JUMP_Z,"
                "SURFACE_HEAT_FLUX,SURFACE_HEAT_FLUX_GLOBAL" );
        tParameterList( 0 )( 0 ).set( "Field_Type", 
                "NODAL,NODAL,NODAL,NODAL,"
                "ELEMENTAL_AVG,ELEMENTAL_AVG,"
                "GLOBAL,GLOBAL,"
                "FACETED_AVG,FACETED_AVG,FACETED_AVG,"
                "FACETED_AVG,FACETED_AVG,FACETED_AVG,"
                "FACETED_AVG,GLOBAL" );
        tParameterList( 0 )( 0 ).set( "IQI_Names", 
                "IQI_Disp_X,IQI_Disp_Y,IQI_Disp_Z,IQI_Temp,"
                "IQI_Stress_VM_Fiber,IQI_Stress_VM_Filler,"
                "IQI_Strain_Energy_Fiber,IQI_Strain_Energy_Filler,"
                "IQI_Traction_Fiber_X,IQI_Traction_Fiber_Y,IQI_Traction_Fiber_Z,"
                "IQI_Traction_Jump_X,IQI_Traction_Jump_Y,IQI_Traction_Jump_Z,"
                "IQI_Surface_Heat_Flux,IQI_Surface_Heat_Flux" );
        tParameterList( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterList( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    MORISGENERALParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
    }

    /* ------------------------------------------------------------------------ */

}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
