// standard library
#include <string>
#include <iostream>

// Type definitions
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// modules
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"

// parameter lists
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"

// LinAlg Operators
#include "fn_equal_to.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

// Parsing Tools
#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_assert.hpp"
#include "fn_stringify_matrix.hpp"


#ifdef __cplusplus
extern "C" {

#endif
//------------------------------------------------------------------------------
namespace moris
{   

    /* ----------------------- */
    // Geometric & Mesh configuration
    
    // size of the box
    real tLenX = 2.0;
    real tLenY = 1.0;
    real tLenZ = 1.0;
    
    // define the element size
    real tH = 0.1;
    
    // distance of two off-center planes from center plane
    real tR = i-00R * tH;
    
    // offset of center planes from x,y,z=0
    real tX0 = i-00X * tH;
    real tY0 = i-00Y * tH;
    real tZ0 = i-00Z * tH;
    
    // domain of the mesh
    real tDomainX = tLenX + 4.0 * tH;
    real tDomainY = tLenY + 4.0 * tH;
    real tDomainZ = tLenZ + 4.0 * tH;
    std::string tDomainSize = ios::stringify( tDomainX ) + "," + ios::stringify( tDomainY ) + "," + ios::stringify( tDomainZ );
    
    // size of the mesh
    uint tNumElemsX = (uint) std::floor( 1e-12 + tDomainX / tH );
    uint tNumElemsY = (uint) std::floor( 1e-12 + tDomainY / tH );
    uint tNumElemsZ = (uint) std::floor( 1e-12 + tDomainZ / tH );
    std::string tMeshSize = ios::stringify( tNumElemsX ) + "," + ios::stringify( tNumElemsY ) + "," + ios::stringify( tNumElemsZ );
    
    // mesh offset
    real tOffX = tX0 - ( tDomainX / 2.0 );
    real tOffY = tY0 - ( tDomainY / 2.0 );
    real tOffZ = tZ0 - ( tDomainZ / 2.0 );
    std::string tOffSet = ios::stringify( tOffX ) + "," + ios::stringify( tOffY ) + "," + ios::stringify( tOffZ );
    
    // material phases
    // Vector< std::string > tMaterialNames = { "2", "3", "5", "6" };
    // Vector< std::string > tYoungsModuli  = { "1000", "3000", "2000", "5000" };
    // Vector< std::string > tNus           = { "0.26", "0.32", "0.28", "0.34" }; 
    Vector< std::string > tMaterialNames = { "1", "2", "3", "4", "5", "6", "7", "8" };
    Vector< std::string > tYoungsModuli  = { "1000", "3000", "2000", "5000", "4000", "2500", "1500", "3500" };
    Vector< std::string > tNus           = { "0.26", "0.32", "0.28", "0.34", "0.27", "0.33", "0.31", "0.29" }; 
    uint tNumMaterials = tMaterialNames.size();
    std::string tLeftBoundaryName = "10";   // block side ordinal 4
    std::string tRightBoundaryName = "11";  // block side ordinal 2
    // Vector< std::string > tLeftMaterialNames = { "5", "8" };
    // Vector< std::string > tRightMaterialNames = { "4", "7" };
    Vector< std::string > tLeftMaterialNames = { "3", "6" };
    Vector< std::string > tRightMaterialNames = { "2", "5" };
    
    uint
    get_number_interface_pairs()
    {
        return ( tNumMaterials * tNumMaterials - tNumMaterials ) / 2;
    }
    
    Vector< std::string >
    get_interface_set_names()
    {        
        // number of possible combinations
        uint tNumCombinations = get_number_interface_pairs();
        
        // initialize output vector
        Vector< std::string > tInterfaceNames( tNumCombinations );
        
        // initialize cursor
        uint iCombo = 0;
        
        // construct the interface names (low-to-high convention)
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {            
            for ( uint iNeighborMaterial = iMaterial + 1; iNeighborMaterial < tNumMaterials; iNeighborMaterial++ )
            {
                tInterfaceNames( iCombo++ ) = "dbl_iside_p0_" + tMaterialNames( iMaterial ) + "_p1_" + tMaterialNames( iNeighborMaterial );
            }
        }
        
        return tInterfaceNames;
    }
    
    Vector< std::string >
    get_bulk_set_names()
    {
        // initialize output vector        
        Vector< std::string > tBulkNames( tNumMaterials );
        
        // construct the bulk names 
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {            
            tBulkNames( iMaterial ) = "HMR_dummy_c_p" + tMaterialNames( iMaterial ) + ",HMR_dummy_n_p" + tMaterialNames( iMaterial );
        }

        return tBulkNames;
    }
    
    Vector< std::string >
    get_ghost_set_names()
    {
        // initialize output vector
        Vector< std::string > tGhostNames( tNumMaterials );
        
        // construct the bulk names 
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {            
            tGhostNames( iMaterial ) = "ghost_p" + tMaterialNames( iMaterial );
        }
        
        return tGhostNames;
    }
    
    Vector< std::string >
    get_dirichlet_set_names()
    {
        // initialize output vector        
        Vector< std::string > tDirichletNames( tLeftMaterialNames.size() );
        
        // construct the bulk names 
        for ( uint iMaterial = 0; iMaterial < tLeftMaterialNames.size(); iMaterial++ )
        {            
            tDirichletNames( iMaterial ) = "iside_b0_" + tLeftMaterialNames( iMaterial ) + "_b1_" + tLeftBoundaryName;
        }

        return tDirichletNames;
    }
    
    Vector< std::string >
    get_neumann_set_names()
    {
        // initialize output vector        
        Vector< std::string > tNeumannNames( tLeftMaterialNames.size() );
        
        // construct the bulk names 
        for ( uint iMaterial = 0; iMaterial < tRightMaterialNames.size(); iMaterial++ )
        {            
            tNeumannNames( iMaterial ) = "iside_b0_" + tRightMaterialNames( iMaterial ) + "_b1_" + tRightBoundaryName;
        }

        return tNeumannNames;
    }
    
    std::string
    get_all_bulk_set_names_as_one_string()
    {
        Vector< std::string > tBulkSetNames = get_bulk_set_names();
        std::string tAllBulkSets = tBulkSetNames( 0 );
        for ( uint i = 1; i < tBulkSetNames.size(); i++ )
        {
            tAllBulkSets += "," + tBulkSetNames( i );
        }
        return tAllBulkSets;
    }
    
    /* ----------------------- */
    // Other Parameters
    
    // print VIS output
    bool tHaveVisOutput = false;
    
    // define the discretization
    std::string tNumRefinements = "0";
    bool tIsEnriched = true;
    
    // Ghost
    std::string tGhostPenalty = "0.01";
    bool tUseGhost = i-00G;    
    
    // names of the mesh sets used
    // std::string tBulkSets      = "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1,HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p3,HMR_dummy_c_p3,HMR_dummy_n_p4,HMR_dummy_c_p4,HMR_dummy_n_p5,HMR_dummy_c_p5";
    // std::string tGhostSets     = "ghost_p0,ghost_p1,ghost_p2,ghost_p3,ghost_p4,ghost_p5";
    // std::string tDirichletSets = "SideSet_4_n_p0,SideSet_4_c_p0,SideSet_4_n_p3,SideSet_4_c_p3";
    // std::string tNeumannSets   = "SideSet_2_n_p2,SideSet_2_c_p2,SideSet_2_n_p5,SideSet_2_c_p5";
    // std::string tInterfaceSets = "dbl_iside_p0_0_p1_1,dbl_iside_p0_1_p1_2,dbl_iside_p0_4_p1_3,dbl_iside_p0_5_p1_4,dbl_iside_p0_0_p1_3,dbl_iside_p0_0_p1_4,dbl_iside_p0_5_p1_1,dbl_iside_p0_5_p1_2";    

    Vector< std::string > tBulkSets      = get_bulk_set_names();
    Vector< std::string > tGhostSets     = get_ghost_set_names();
    Vector< std::string > tDirichletSets = get_dirichlet_set_names();
    Vector< std::string > tNeumannSets   = get_neumann_set_names();
    Vector< std::string > tInterfaceSets = get_interface_set_names(); 

    std::string tAllBulkSets = get_all_bulk_set_names_as_one_string();


    // polynomial of the background interpolation
    std::string tPolyOrder = "2";
    
    // intersection computation mode
    bool tUseMultiLinear = false;

    // suppress VIS output
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tHaveVisOutput;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );
        tParameterList( 0 )( 0 ) = prm::create_opt_problem_parameter_list();
    }

    //------------------------------------------------------------------------------

    void
    HMRParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );

        tParameterList( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterList( 0 )( 0 ).set( "number_of_elements_per_dimension", tMeshSize );
        tParameterList( 0 )( 0 ).set( "domain_dimensions", tDomainSize );
        tParameterList( 0 )( 0 ).set( "domain_offset", tOffSet );

        tParameterList( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4,5,6" );
        tParameterList( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterList( 0 )( 0 ).set( "lagrange_orders", tPolyOrder );
        tParameterList( 0 )( 0 ).set( "lagrange_pattern", "0" );

        tParameterList( 0 )( 0 ).set( "bspline_orders", tPolyOrder );
        tParameterList( 0 )( 0 ).set( "bspline_pattern", "0" );
        tParameterList( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
        
        tParameterList( 0 )( 0 ).set( "initial_refinement", tNumRefinements );
        tParameterList( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
    }

    //------------------------------------------------------------------------------

    void
    XTKParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );

        tParameterList( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterList( 0 )( 0 ).set( "decompose", true );
        tParameterList( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterList( 0 )( 0 ).set( "enrich", true );
        tParameterList( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterList( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterList( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterList( 0 )( 0 ).set( "visualize_ghost", false );
        tParameterList( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        tParameterList( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true ); //! TURN THIS BACK ON !
        tParameterList( 0 )( 0 ).set( "output_cut_ig_mesh", false );
        
        // for debugging
        tParameterList( 0 )( 0 ).set( "only_generate_xtk_temp", false );
        tParameterList( 0 )( 0 ).set( "triangulate_all", false );
        tParameterList( 0 )( 0 ).set( "verbose", false );
        // tParameterList( 0 )( 0 ).set( "verbose_level", moris::uint( 1 ) );
        tParameterList( 0 )( 0 ).set( "print_enriched_ig_mesh", false );

        if( !tIsEnriched )
        {
            tParameterList( 0 )( 0 ).set( "unenriched_mesh_indices", "0" );
        }
    }

    //------------------------------------------------------------------------------

    void
    GENParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 3 );
        tParameterList( 0 ).push_back( prm::create_gen_parameter_list() );

        // 3-beam version phase table
        Matrix< DDUMat > tPhaseMap( 2048, 1, 0 );
        
        // left side
        tPhaseMap( 128 ) = 10;
        tPhaseMap( 136 ) = 10;
        tPhaseMap( 144 ) = 10;
        tPhaseMap( 152 ) = 10;
        
        // right side
        tPhaseMap( 519 ) = 11;
        tPhaseMap( 527 ) = 11;
        tPhaseMap( 535 ) = 11;
        tPhaseMap( 543 ) = 11;
        
        // bottom bricks
        tPhaseMap(  0 ) = 6;
        tPhaseMap(  4 ) = 7;
        tPhaseMap(  6 ) = 1;
        tPhaseMap(  7 ) = 2;
        tPhaseMap(  8 ) = 6;
        tPhaseMap( 12 ) = 8;
        tPhaseMap( 14 ) = 2;
        tPhaseMap( 15 ) = 2;
        
        // top bricks
        tPhaseMap( 16 ) = 3;
        tPhaseMap( 20 ) = 3;
        tPhaseMap( 22 ) = 1;
        tPhaseMap( 23 ) = 5;
        tPhaseMap( 24 ) = 3;
        tPhaseMap( 28 ) = 8;
        tPhaseMap( 30 ) = 4;
        tPhaseMap( 31 ) = 5;
        
        std::string tPhaseMapString = moris::ios::stringify( tPhaseMap );
        tParameterList( 0 )( 0 ).set( "phase_table", tPhaseMapString );
        tParameterList( 0 )( 0 ).set( "print_phase_table", false );

        // init geometry counter
        uint tGeoCounter = 0;
        
        // Outside Box - Ordinal 1 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", - tLenY / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", - 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Outside Box - Ordinal 2 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", + tLenX / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", + 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Outside Box - Ordinal 3 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", + tLenY / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", + 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Outside Box - Ordinal 4 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", - tLenX / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", - 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Outside Box - Ordinal 5 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", - tLenZ / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", - 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Outside Box - Ordinal 6 (of 6)
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", + tLenZ / 2.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", + 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Middle Plane: X-Y-Plane
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Middle Plane: X-Z-Plane
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Left Plane
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", -tR );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Center Plane
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        
        // Right Plane
        tParameterList( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::PLANE ) );
        tParameterList( 1 )( tGeoCounter ).set( "center_x", tR );
        tParameterList( 1 )( tGeoCounter ).set( "center_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "center_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "normal_z", 0.0 );
        tParameterList( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );
        tGeoCounter++;
        

    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties of bars
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropUnit" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;
        
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0;0.0" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropNeumann" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0;1.0;1.0" );
        tPropCounter++;
        
        // create material properties
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" + tMaterialNames( iMaterial ) );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters", tYoungsModuli( iMaterial ) );
            tPropCounter++;

            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" + tMaterialNames( iMaterial ) );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters", tNus( iMaterial ) );
            tPropCounter++;
        }

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create constitutive models
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso" + tMaterialNames( iMaterial ) );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", (uint) fem::Constitutive_Type::STRUC_LIN_ISO );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties", 
                    "PropYoungs" + tMaterialNames( iMaterial ) + ",YoungsModulus; "
                    "PropPoisson" + tMaterialNames( iMaterial ) + ",PoissonRatio" );
            tCMCounter++;
        }

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // dirichlet nitsche stabilization parameters
        for ( uint iMaterial = 0; iMaterial < tLeftMaterialNames.size(); iMaterial++ )
        {
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Dirichlet_" + tLeftMaterialNames( iMaterial ) );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type", (uint) fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", "10.0" );
            tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs" + tLeftMaterialNames( iMaterial ) + ",Material" );
            tSPCounter++;
        }
        
        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {
            // ghost parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhost" + tMaterialNames( iMaterial ) );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type", (uint) fem::Stabilization_Type::GHOST_DISPL );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
            tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs" + tMaterialNames( iMaterial ) + ",Material" );
            tSPCounter++;
                   
            for ( uint jMaterial = iMaterial + 1; jMaterial < tNumMaterials; jMaterial++ )
            {
                // dirichlet interface parameter
                tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
                tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SP_Nitsche_Interface_p0_" + tMaterialNames( iMaterial ) + "_p1_" + tMaterialNames( jMaterial ) );
                tParameterList( 2 )( tSPCounter ).set( "stabilization_type", (uint) fem::Stabilization_Type::NITSCHE_INTERFACE );
                tParameterList( 2 )( tSPCounter ).set( "function_parameters", "10.0" );
                tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs" + tMaterialNames( iMaterial ) + ",Material" );
                tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropYoungs" + tMaterialNames( jMaterial ) + ",Material" );
                tSPCounter++;
            }
        }

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        for ( uint iMaterial = 0; iMaterial < tNumMaterials; iMaterial++ )
        {
            // BULK IWGs
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG_Bulk_" + tMaterialNames( iMaterial ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_BULK );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso" + tMaterialNames( iMaterial ) + ",ElastLinIso" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulkSets( iMaterial ) );
            tIWGCounter++;
            
            // GHOST IWGs
            if ( tUseGhost )
            {
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG_Ghost_" + tMaterialNames( iMaterial ) );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint) fem::IWG_Type::GHOST_NORMAL_FIELD );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPGhost" + tMaterialNames( iMaterial ) + ",GhostSP" );
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tGhostSets( iMaterial ) );
                tIWGCounter++;
            }
            
            // INTERFACE IWGs
            for ( uint jMaterial = iMaterial + 1; jMaterial < tNumMaterials; jMaterial++ )
            {
                std::string tName = "_p0_" + tMaterialNames( iMaterial ) + "_p1_" + tMaterialNames( jMaterial );
                std::string tSetName = "dbl_iside" + tName;
                
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG_Interface" + tName );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY,UZ" );
                tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso" + tMaterialNames( iMaterial ) + ",ElastLinIso" );
                tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", "CMStrucLinIso" + tMaterialNames( jMaterial ) + ",ElastLinIso" );
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Interface" + tName + ",NitscheInterface" );
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSetName );
                tIWGCounter++;
            }
        }

        for ( uint iMaterial = 0; iMaterial < tLeftMaterialNames.size(); iMaterial++ )
        { 
            // DIRICHLET IWGs
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG_Dirichlet_" + tLeftMaterialNames( iMaterial ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichlet,Dirichlet" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso" + tLeftMaterialNames( iMaterial ) + ",ElastLinIso" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SP_Nitsche_Dirichlet_" + tLeftMaterialNames( iMaterial ) + ",DirichletNitsche" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tDirichletSets( iMaterial ) );
            tIWGCounter++;
        }

        for ( uint iMaterial = 0; iMaterial < tRightMaterialNames.size(); iMaterial++ )
        {
            // NEUMANN IWGs
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWG_Traction_" + tRightMaterialNames( iMaterial ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint) fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropNeumann,Traction" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tNeumannSets( iMaterial ) );
            tIWGCounter++;
        }
       

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Bulk_DISP_X" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tAllBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Bulk_DISP_Y" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tAllBulkSets );
        tIQICounter++;
        
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Bulk_DISP_Z" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::DOF );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tAllBulkSets );
        tIQICounter++;

        // // nodal von-mises stresses
        // tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_VM_Stress" );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::VON_MISES_STRESS );
        // tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        // tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        // tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        // tIQICounter++;
        // 
        // // nodal stresses
        // tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Normal_Stress_X" );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        // tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        // tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        // tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        // tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        // tIQICounter++;
        // 
        // tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Normal_Stress_Y" );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        // tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        // tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        // tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        // tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        // tIQICounter++;
        // 
        // tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQI_Normal_Stress_Z" );
        // tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint) fem::IQI_Type::NORMAL_STRESS );
        // tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY,UZ" );
        // tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        // tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
        // tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        // tIQICounter++;


        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 8 );
        
        /* ----------------------- */
        // Linear Pre-conditioners
        
        // empty pre-conditioner
        tParameterList( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        tParameterList( 7 )( 0 ).set( "PCType", "none" );

        // MUMPS direct solver - used as pre-conditioner (which already fully solves the linear system)
        tParameterList( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        tParameterList( 7 )( 1 ).set( "PCType", "mumps" );
        
        /* ----------------------- */
        // Linear Algorithms

        // TODO: why put a solver on top of an already solved system?
        tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        tParameterList( 0 )( 0 ).set( "KSPType", "fgmres" );
        tParameterList( 0 )( 0 ).set( "preconditioners", "1" );

        // empty linear solver algorithm (pre-conditioner is also empty)
        tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        tParameterList( 0 )( 1 ).set( "KSPType", "preonly" );
        tParameterList( 0 )( 1 ).set( "preconditioners", "0" );

        // linear solver which simply passes on the solution computed by the MUMPS direct solver (which has been re-packaged as a pre-conditioner)
        tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        tParameterList( 0 )( 2 ).set( "KSPType", "preonly" );
        tParameterList( 0 )( 2 ).set( "preconditioners", "1" );

        // find min eigenvalue - needs a solved linear system
        tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        tParameterList( 0 )( 3 ).set( "Eigen_Algorithm", "krylovschur" );
        tParameterList( 0 )( 3 ).set( "Which", "LM" );
        tParameterList( 0 )( 3 ).set( "Num_Eig_Vals", 1 );
        tParameterList( 0 )( 3 ).set( "STType", "shift_invert" ); // search from lowest to highest eigenvalue (H-to-L is the default)
        tParameterList( 0 )( 3 ).set( "sub_linear_solver", "0" );
        tParameterList( 0 )( 3 ).set( "is_symmetric", false );
        tParameterList( 0 )( 3 ).set( "Update_Flag", false ); // whether or not to pass the eigenvector as a solution field to FEM, turn on flag in MSI to use

        // find max eigenvalue - does not need a solved linear system
        tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        tParameterList( 0 )( 4 ).set( "Eigen_Algorithm", "krylovschur" );
        tParameterList( 0 )( 4 ).set( "Which", "LM" );
        tParameterList( 0 )( 4 ).set( "Num_Eig_Vals", 1 );
        tParameterList( 0 )( 4 ).set( "sub_linear_solver", "1" );
        tParameterList( 0 )( 4 ).set( "is_symmetric", false );
        tParameterList( 0 )( 4 ).set( "Update_Flag", false ); // whether or not to pass the eigenvector as a solution field to FEM, turn on flag in MSI to use

        /* ----------------------- */
        // Linear Solvers

        // TODO: What does this solver solve? 
        tParameterList( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
        tParameterList( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        // Solves MIN eigenvalue problem
        tParameterList( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() ); 
        tParameterList( 1 )( 1 ).set( "DLA_Linear_solver_algorithms", "3" );
        tParameterList( 1 )( 1 ).set( "RHS_Matrix_Type", "IdentityMat" );

        // Solves MAX eigenvalue problem
        tParameterList( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
        tParameterList( 1 )( 2 ).set( "DLA_Linear_solver_algorithms", "4" );
        tParameterList( 1 )( 2 ).set( "RHS_Matrix_Type", "IdentityMat" );

        /* ----------------------- */
        // Non-Linear Algorithms

        // TODO: what is this algorithm used for?
        tParameterList( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterList( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-09 );
        tParameterList( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterList( 2 )( 0 ).set( "NLA_max_iter", 1 );
        // tParameterList( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );

        // Solve the MIN eigenvalue problem
        tParameterList( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterList( 2 )( 1 ).set( "NLA_Solver_Implementation", (uint) moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterList( 2 )( 1 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterList( 2 )( 1 ).set( "NLA_max_iter", 1 );
        tParameterList( 2 )( 1 ).set( "NLA_Linear_solver", 1 );
        tParameterList( 2 )( 1 ).set( "NLA_is_eigen_problem", true );

        // Solve the MAX eigenvalue problem
        tParameterList( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterList( 2 )( 2 ).set( "NLA_Solver_Implementation", (uint) moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterList( 2 )( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterList( 2 )( 2 ).set( "NLA_max_iter", 1 );
        tParameterList( 2 )( 2 ).set( "NLA_Linear_solver", 2 );
        tParameterList( 2 )( 2 ).set( "NLA_is_eigen_problem", true );

        /* ----------------------- */
        // Non-Linear Solvers

        // solve the complete eigenvalue problem, i.e. both the MIN and MAX EV problems
        tParameterList( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterList( 3 )( 0 ).set( "NLA_Solver_Implementation", (uint) moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterList( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "2,1,0" );
        tParameterList( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY,UZ" );

        /* ----------------------- */
        // Time Solver Algorithm

        // solve both the structural and the eigenvalue problems here
        tParameterList( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );

        /* ----------------------- */
        // Time Solver

        // solve everything
        tParameterList( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterList( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY,UZ" );
        tParameterList( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        /* ----------------------- */
        // Solver Warehouse

        tParameterList( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
        tParameterList( 6 )( 0 ).set( "SOL_TPL_Type", sol::MapType::Petsc ); // choose data-structure that is compatible with PetsC solvers
        // tParameterList( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "SOE" ); // NOTE: turn on "NLA_combined_res_jac_assembly" in NLA to use this

// /* ----------------------- */
// // debug - check result using standard solver setup
// tParameterList( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
// tParameterList( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );
// tParameterList( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
// tParameterList( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-09 );
// tParameterList( 2 )( 0 ).set( "NLA_max_iter", 1 );
// tParameterList( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );
// tParameterList( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
// tParameterList( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY,UZ" );
// tParameterList( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
// tParameterList( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
// tParameterList( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY,UZ" );
// tParameterList( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
// tParameterList( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
// tParameterList( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "SOE" );
// tParameterList( 7 ).push_back( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );

    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );
        tParameterList( 0 ).resize( 1 );
        tParameterList( 0 )( 0 ) = prm::create_msi_parameter_list();
        // tParameterList( 0 )( 0 ).set( "order_adofs_by_host", true ); // is needed for gamg pre-conditioner to converge
        tParameterList( 0 )( 0 ).set( "number_eigen_vectors", 2 ); // NOTE: this is needed when "Update_Flag" is turned on in linear algorithm
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
        tParameterList.resize( 1 );

        tParameterList( 0 ).push_back( prm::create_vis_parameter_list() );
        tParameterList( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "brick_wall.exo" ) );
        tParameterList( 0 )( 0 ).set( "Mesh_Type", (uint) vis::VIS_Mesh_Type::STANDARD );
        tParameterList( 0 )( 0 ).set( "Set_Names", tAllBulkSets );
        tParameterList( 0 )( 0 ).set( "Field_Names", 
                "UX,UY,UZ" );
        tParameterList( 0 )( 0 ).set( "Field_Type", 
                "NODAL,NODAL,NODAL" );
        tParameterList( 0 )( 0 ).set( "IQI_Names", 
                "IQI_Bulk_DISP_X,IQI_Bulk_DISP_Y,IQI_Bulk_DISP_Z" );
        tParameterList( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    //------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( moris::Vector< moris::Vector< Parameter_List > >& tParameterList )
    {
    }

    //------------------------------------------------------------------------------
    
} // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
