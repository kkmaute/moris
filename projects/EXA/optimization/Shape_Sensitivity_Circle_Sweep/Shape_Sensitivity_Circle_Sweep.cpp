#include <string>
#include <iostream>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_PRM_FEM_Parameters.hpp"
#include "cl_PRM_MSI_Parameters.hpp"
#include "cl_PRM_SOL_Parameters.hpp"
#include "cl_PRM_VIS_Parameters.hpp"
#include "cl_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_PRM_XTK_Parameters.hpp"
#include "cl_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"


#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Constant function for properties
    void Func_Const( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                     moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                     moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDSMat> get_constraint_types()
    {
        Matrix<DDSMat> tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_objectives(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tObjectives( 1, 1 );
        tObjectives(0) = aCriteria( 0 ) / 2.7;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_constraints(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tConstraints( 1, 1 );
        tConstraints(0) = aCriteria( 1 )/2.0 - 0.5;

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dobjective_dadv(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dobjective_dcriteria(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1 / 2.7;
        tDObjectiveDCriteria( 1 ) = 0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dconstraint_dadv(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dconstraint_dcriteria(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0;
        tDConstraintDCriteria( 1 ) = 1/1.0;

        return tDConstraintDCriteria;
    }

    void Func_Neumann_U( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                         moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                         moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        if( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.2 )
        {
            aPropMatrix = {{ 0.0 },
                    { aParameters( 0 )( 0 ) }};
        }
        else
        {
            aPropMatrix = {{ 0.0},
                    { 0.0}};
        }
    }



    //--------------------------------------------------------------------------------------------------------------

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string("9,5"));
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                std::string("2,1"));
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    std::string("0.0,0.0") );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  std::string("1,2,3,4"));
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           std::string("0"));

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ));
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string( "1" ));
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string( "0" ));

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", 0 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 1 );

        tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        std::string("conformal") );
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                std::string("bspline") );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       std::string("0") );
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize(1);
        tParameterlist( 1 ).resize(1);
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set("initial_advs"      , "1.0, 0.5, 0.25");
        tParameterlist( 0 )( 0 ).set("lower_bounds"      , "0.0, 0.0, 0.0");
        tParameterlist( 0 )( 0 ).set("upper_bounds"      , "2.0, 1.0, 0.5");
        tParameterlist( 0 )( 0 ).set("IQI_types"         , "IQIBulkStrainEnergy,IQIBulkVolume");
        tParameterlist( 0 )( 0 ).set("PDV_types"         , "X_COORDINATE,Y_COORDINATE");

        // Geometry parameter lists
        tParameterlist(1)(0) = prm::create_geometry_parameter_list();
        tParameterlist(1)(0).set("type", "circle");
        tParameterlist(1)(0).set("geometry_variable_indices", "0, 1, 2");
        tParameterlist(1)(0).set("adv_indices", "0, 1, 2");

    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 5 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDensity") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropYoungs") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropFlux") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("10.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropDirichletU") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.0;0.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropNeumann") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Neumann_U") );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropPoisson") );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.0") );
        tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMStrucLinIso1") );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        std::string("PropYoungs,YoungsModulus;PropPoisson,PoissonRatio"));
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPNitscheTemp") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropYoungs,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create parameter list for IWG 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGBulkU_1") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMStrucLinIso1,ElastLinIso") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIWGCounter++;

        // create parameter list for IWG 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGDirichletU") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropDirichletU,Dirichlet") );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMStrucLinIso1,ElastLinIso") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheTemp,DirichletNitsche") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("SideSet_4_n_p1") );
        tIWGCounter++;

        // create parameter list for IWG 3
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGNeumannFlux") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("UX,UY") );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropNeumann,Neumann") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("SideSet_2_n_p1") );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMP") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UX ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("UX,UY") );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkStrainEnergy") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("UX,UY") );
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", std::string("CMStrucLinIso1,Elast") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIQICounter++;

        // create parameter list for IQI 4
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkVolume") );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          std::string("PropDensity,Density") );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p1") );
        tIQICounter++;
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        //tParameterlist( 2 )( 0 ).set("NLA_rel_residual", 1e-07 );
        tParameterlist( 2 )( 0 ).set("NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , std::string("UX,UY") );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",     1 );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",         1.0 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           std::string("UX,UY") );
        //tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", std::string("UX,UY,0.0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     std::string("0") );
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",     std::string("Output_Criterion") );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "shape_opt_test.exo" ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p1,HMR_dummy_c_p1" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "VOLUME" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL" ) );
        tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "VOLUME" ) );
    }


    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist(0)(0) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist(0)(0).set("is_optimization_problem", true);
        tParameterlist(0)(0).set("problem", "user_defined");
        tParameterlist(0)(0).set("library", "Shape_Sensitivity_Circle_Sweep.so");

        tParameterlist(2)(0) = moris::prm::create_sweep_parameter_list();
        tParameterlist(2)(0).set("print", true);
        tParameterlist(2)(0).set("hdf5_path", "shape_opt.hdf5");
        tParameterlist(2)(0).set("evaluate_objective_gradients", true);
        tParameterlist(2)(0).set("evaluate_constraint_gradients", true);
        tParameterlist(2)(0).set("num_evaluations_per_adv", "1, 1, 2");
        tParameterlist(2)(0).set("include_bounds", false);

    }

    //------------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
