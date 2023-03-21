/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Cluster_Measure_2Mat_SA.cpp
 *
 */

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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// problem dimension: 2D or 3D
extern uint gDim;

//---------------------------------------------------------------

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // ghost
    bool tUseGhost = true;

    // ghost penalty value
    std::string tGhostPenalty = "0.001";

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName    = "Cluster_Measure_2Mat_SA" + std::to_string(gDim)
                                   + ".exo";
    std::string tGENOutputFileName = "Cluster_Measure_2Mat_SA" + std::to_string(gDim)
                                   + "_GEN.exo";
    std::string tSoFileName        = "Cluster_Measure_2Mat_SA.so";
    std::string tHdf5FileName      = "Cluster_Measure_2Mat_SA" + std::to_string(gDim)
                                   + ".hdf5";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information
    std::string tPhase1        = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tPhase0        = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tTotalDomain   = tPhase0 + "," + tPhase1;

    // FD
    std::string tFDEpsilon = "1e-6";
    std::string tWhichFD   = "all";
    uint        tFEMWhichFDForSA = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );
    uint        tFEMWhichFDForFA = static_cast< uint >( fem::FDScheme_Type::POINT_3_CENTRAL );
    uint        tFEMWhichPerturbation = static_cast< uint >( fem::Perturbation_Type::ABSOLUTE );

    /* ------------------------------------------------------------------------ */
    // boundary condition

    // prescribed temperature
    moris::real sTpre = 0.0;

    // flux at tip fac
    moris::real sP2 = 5.0;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // capacity
    std::string sCap1 = "0.0";

    // density
    std::string tDens1 = "0.0";

    // conductivity
    moris::real sK1 = 1.0;

    // body flux
    moris::real  sQ1 = 1.0;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim     = gDim == 2 ? "1, 1"     : "1, 1, 1";
    std::string tDomainDims         = gDim == 2 ? "1.0, 1.0" : "1.0, 1.0, 1.0";
    std::string tDomainOffset       = gDim == 2 ? "0.0, 0.0" : "0.0, 0.0, 0.0";
    std::string tDomainSidesets     = gDim == 2 ? "1,2,3,4"  : "1,2,3,4,5,6";

    int tRefineBuffer      = 1;

    /* ------------------------------------------------------------------------ */
    // Phase assignment

    std::string tGetPhaseIndex = gDim == 2 ? "get_phase_index_2d" : "get_phase_index_3d";

    uint get_phase_index_2d( const Bitset<1>& aGeometrySigns )
    {
        // by default in void
        uint tPhaseIndex = 0;

        // Phase solid 1
        if(  !aGeometrySigns.test( 0 ) )
        {
            tPhaseIndex = 1;
        }

        return tPhaseIndex;
    }

    uint get_phase_index_3d( const Bitset<1>& aGeometrySigns )
    {
        // by default in void
        uint tPhaseIndex = 0;

        // Phase solid 1
        if(  !aGeometrySigns.test( 0 ) )
        {
            tPhaseIndex = 1;
        }

        return tPhaseIndex;
    }

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int tNLA_max_iter = 2;

    int tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    moris::real tMinLevs = 1.0e-8;

    /* ------------------------------------------------------------------------ */
    // Function for plane

     moris::real Plane(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const moris::Cell< moris::real* > & aGeometryParameters )
    {
        // coordinates in x and y directions
        moris::real tX = aCoordinates(0);
        moris::real tY = aCoordinates(1);

        // get interface position
        moris::real tXInterface = *(aGeometryParameters(0));
        moris::real tYInterface = *(aGeometryParameters(1));

        // compute signed-distance field
        moris::real tVal = ( tX - tXInterface ) + ( tY - tYInterface );

        // clean return value to return non-zero value
        return std::abs(tVal) < tMinLevs ? tMinLevs : tVal;
    }

    /* ------------------------------------------------------------------------ */
     // Derivative function for plane

    void PlaneGrad(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const moris::Cell< moris::real* >& aGeometryParameters,
            moris::Matrix< DDRMat >&           aSensitivities)
    {
        // create copy of geometry parameter cell
        moris::Cell<moris::real*> tGeometryParameters = aGeometryParameters;

        // define finite difference perturbation size
        const real tPerturbation = 1.0e-6;

        // set size of sensitivity vector
        aSensitivities.set_size(1, 2);

        // initialize adv counter
        uint iv=0;

        // loop over all advs
        for (uint ip=0;ip<2;++ip)
        {
            // extract nominal value used for perturbations
            real tNominalValue = *(tGeometryParameters(iv));

            // set pointer to pertubed value
            tGeometryParameters(iv) = &tNominalValue;

            // positive perturbation
            tNominalValue += tPerturbation;
            real tPosValue = Plane(aCoordinates,tGeometryParameters);

            // positive perturbation
            tNominalValue -= 2.0*tPerturbation;
            real tNegValue = Plane(aCoordinates,tGeometryParameters);

            // restore nominal value
            tGeometryParameters(iv) = aGeometryParameters(iv);

            // compute sensitivities
            aSensitivities(iv) = (tPosValue - tNegValue)/(2.0*tPerturbation);

            // increase counter
            iv++;
        }
    }

    /* ------------------------------------------------------------------------ */

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

   /* ------------------------------------------------------------------------ */

   uint tNumConstraints = 8;
    Matrix<DDSMat> get_constraint_types()
    {
        Matrix<DDSMat> tConstraintTypes( 1, tNumConstraints, 1 );

        return tConstraintTypes;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_objectives(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 );

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_constraints(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tConstraints( 1, tNumConstraints );
        tConstraints( 0 ) = aCriteria( 1 );
        tConstraints( 1 ) = aCriteria( 2 );
        tConstraints( 2 ) = aCriteria( 3 );
        tConstraints( 3 ) = aCriteria( 4 );
        tConstraints( 4 ) = aCriteria( 5 );
        tConstraints( 5 ) = aCriteria( 6 );
        tConstraints( 6 ) = aCriteria( 7 );
        tConstraints( 7 ) = aCriteria( 8 );

        return tConstraints;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_dobjective_dadv(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_dobjective_dcriteria(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );
        tDObjectiveDCriteria( 0 ) = 1;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_dconstraint_dadv(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDADV( tNumConstraints, aADVs.numel(), 0.0 );
        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix<DDRMat> compute_dconstraint_dcriteria(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDCriteria( tNumConstraints, aCriteria.numel(), 0.0 );
        tDConstraintDCriteria( 0, 1 ) = 1.0;
        tDConstraintDCriteria( 1, 2 ) = 1.0;
        tDConstraintDCriteria( 2, 3 ) = 1.0;
        tDConstraintDCriteria( 3, 4 ) = 1.0;
        tDConstraintDCriteria( 4, 5 ) = 1.0;
        tDConstraintDCriteria( 5, 6 ) = 1.0;
        tDConstraintDCriteria( 6, 7 ) = 1.0;
        tDConstraintDCriteria( 7, 8 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 )( 0 ).set("is_optimization_problem", true);
        tParameterlist( 0 )( 0 ).set("problem", "user_defined");
        tParameterlist( 0 )( 0 ).set("library", tSoFileName );

        tParameterlist( 1 ).resize( 0 );

        tParameterlist( 2 ).push_back( moris::prm::create_sweep_parameter_list() );
        tParameterlist( 2 )( 0 ).set("hdf5_path", tHdf5FileName );
        tParameterlist( 2 )( 0 ).set("evaluate_objective_gradients", true);
        tParameterlist( 2 )( 0 ).set("evaluate_constraint_gradients", true);
        tParameterlist( 2 )( 0 ).set("num_evaluations_per_adv", "1");
        tParameterlist( 2 )( 0 ).set("include_bounds", false);
        tParameterlist( 2 )( 0 ).set("finite_difference_type", tWhichFD);
        tParameterlist( 2 )( 0 ).set("finite_difference_epsilons",tFDEpsilon);
        tParameterlist( 2 )( 0 ).set("print",true);
    }

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  tDomainSidesets);
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::to_string(gInterpolationOrder) );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern",  "0" )  ;
        tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::to_string(gInterpolationOrder) );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",   "0" )  ;

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0") ;

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        "conformal") ;
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                "bspline") ;
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       "0") ;
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        tParameterList.resize( 3 );

        // main gen parameter list
        tParameterList( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterList( 0 )( 0 ).set("isocontour_tolerance", 1e-12 );
        tParameterList( 0 )( 0 ).set("number_of_phases"    , 2 );
        tParameterList( 0 )( 0 ).set("phase_function_name" , tGetPhaseIndex );
        tParameterList( 0 )( 0 ).set("initial_advs", "0.35,1.0" );
        tParameterList( 0 )( 0 ).set("lower_bounds", "0.35,1.0" );
        tParameterList( 0 )( 0 ).set("upper_bounds", "0.35,1.0" );
        tParameterList( 0 )( 0 ).set("IQI_types"   ,
                                     "IQIVolumeInterface1,"
                                     "IQIVolumeMat0,IQIVolumeMat1,"
                                     "IQISidesetMeasure,"
                                     "IQISidesetLengthMeasure,"
                                     "IQIBulkMeasure,"
                                     "IQIBulkLengthMeasure,"
                                     "IQIBulkStrainEnergy0,IQIBulkStrainEnergy1" );
        //tParameterlist( 0 )( 0 ).set("output_mesh_file", tGENOutputFile );

        // init geometry counter
        uint tGeoCounter = 0;

        // interface plane
        tParameterList( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterList( 1 )( tGeoCounter ).set( "field_function_name",       "Plane" );
        tParameterList( 1 )( tGeoCounter ).set( "sensitivity_function_name", "PlaneGrad");
        tParameterList( 1 )( tGeoCounter ).set( "field_variable_indices",    "0,1");
        tParameterList( 1 )( tGeoCounter ).set( "adv_indices",               "0,1");
        tGeoCounter++;
    }

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
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
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name",       "PhaseMaterial0" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices",    "0"  );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name",       "PhaseMaterial1" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices",    "1"  );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name",       "PhaseAll" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices",    "0,1"  );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // properties for material 0
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropDensity0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      tDens1 );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropCapacity0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      sCap1 );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropConductivity0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(sK1) );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropHeatLoad0") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(sQ1) );
        tPropCounter++;

        // properties for material 1
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropDensity1") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      tDens1 );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropCapacity1") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      sCap1 );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropConductivity1") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(sK1) );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropHeatLoad1") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(sQ1) );
        tPropCounter++;

        // surface flux
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropSurfaceFlux") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      std::to_string(sP2) );
        tPropCounter++;

        // temperature at back surface
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name",            "PropImposedTemperature") ;
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters",      "2.0") ;
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create diffusion CM for material 0
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion0") ;
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name",        "PhaseMaterial0") ;
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity0 , Conductivity;"
                "PropDensity0      , Density;"
                "PropCapacity0     , HeatCapacity") ;
        tCMCounter++;

        // create diffusion CM for material 1
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMDiffusion1") ;
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name",        "PhaseMaterial1") ;
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity") ;
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // Nitsche for Dirichlet on material 0
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPNitscheTemp0") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",   "PropConductivity0,Material") ;
        tSPCounter++;

        // Nitsche for Dirichlet on material 1
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPNitscheTemp1") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",   "PropConductivity1,Material") ;
        tSPCounter++;

        // Nitsche stabilization parameter for mat0 - mat1
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPInterfaceNitsche01") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "slave_phase_name",    "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", "100.0") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",   "PropConductivity0,Material") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "slave_properties",    "PropConductivity1,Material") ;
        tSPCounter++;

        // cluster measure for sideset
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPSidesetMeasure") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MEASURE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "cluster_measures",    std::pair< std::string, std::string >( "CELL_SIDE_MEASURE,PRIMARY,MASTER", "ElementSize" ) );
        tSPCounter++;

        // cluster measure for sideset length
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPSidesetLengthMeasure") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MEASURE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "cluster_measures",    std::pair< std::string, std::string >( "CELL_LENGTH_MEASURE,PRIMARY,MASTER", "ElementSize" ) );
        tSPCounter++;

        // cluster measure for bulk
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPBulkMeasure") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MEASURE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "cluster_measures",    std::pair< std::string, std::string >( "CELL_MEASURE,PRIMARY,MASTER", "ElementSize" ) );
        tSPCounter++;

        // cluster measure for bulk length
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPBulkLengthMeasure") ;
        tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MEASURE ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "cluster_measures",    std::pair< std::string, std::string >( "CELL_LENGTH_MEASURE,PRIMARY,MASTER", "ElementSize" ) );
        tSPCounter++;

        if (tUseGhost)
        {
            // ghost penalty for material 0
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPGPTemp0") ;
            tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial0" );
            tParameterList( tSPIndex )( tSPCounter ).set( "slave_phase_name",   "PhaseMaterial0" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tGhostPenalty ) ;
            tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",   "PropConductivity0,Material") ;
            tParameterList( tSPIndex )( tSPCounter ).set( "slave_properties",   "PropConductivity0,Material") ;
            tSPCounter++;

            // ghost penalty for material 1
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name",  "SPGPTemp1") ;
            tParameterList( tSPIndex )( tSPCounter ).set( "master_phase_name",   "PhaseMaterial1" );
            tParameterList( tSPIndex )( tSPCounter ).set( "slave_phase_name",   "PhaseMaterial1" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tGhostPenalty ) ;
            tParameterList( tSPIndex )( tSPCounter ).set( "master_properties",   "PropConductivity1,Material") ;
            tParameterList( tSPIndex )( tSPCounter ).set( "slave_properties",   "PropConductivity1,Material") ;
            tSPCounter++;
        }

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG for material 0 - bulk diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusion0Bulk") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial0" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion0,Diffusion") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropHeatLoad0,Load") ;
        tIWGCounter++;

        // create IWG for material 1 - bulk diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGDiffusion1Bulk") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropHeatLoad1,Load") ;
        tIWGCounter++;

        // create IWG for interface conditions mat0 - mat1
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGInterface01") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              static_cast< uint >( fem::Element_Type::DOUBLE_SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP");
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial0" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_phase_name",           "PhaseMaterial1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion0,Diffusion");
        tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_constitutive_models",  "CMDiffusion1,Diffusion");
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "SPInterfaceNitsche01,NitscheInterface");
        tIWGCounter++;

        // create IWG for Neumann boundary conditions
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGInletFlux") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropSurfaceFlux,Neumann") ;
        tIWGCounter++;

        // create IWG for Dirichlet boundary conditions on mat 0
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGOutletTemp0") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",
                                                        static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial0" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropImposedTemperature,Dirichlet") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion0,Diffusion") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheTemp0,DirichletNitsche") ;
        tIWGCounter++;

        // create IWG for Dirichlet boundary conditions on mat 1
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGOutletTemp1") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",
                                                        static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals",              "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_properties",          "PropImposedTemperature,Dirichlet") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusion1,Diffusion") ;
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheTemp1,DirichletNitsche") ;
        tIWGCounter++;

        if (tUseGhost)
        {
            // create IWG for ghost on material 0
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGGP0Temp") ;
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",
                                                            static_cast< uint >( fem::Element_Type::DOUBLE_SIDESET ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",
                                                            static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial0" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_phase_name",           "PhaseMaterial0" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp0,GhostSP") ;
            tIWGCounter++;
            // create IWG for ghost on material 1
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name",                   "IWGGP1Temp") ;
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type",
                                                            static_cast< uint >( fem::Element_Type::DOUBLE_SIDESET ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type",
                                                            static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "master_phase_name",          "PhaseMaterial1" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "slave_phase_name",           "PhaseMaterial1" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual",               "TEMP") ;
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters",   "SPGPTemp1,GhostSP") ;
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        // Nodal Temperature IQI
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkTEMP") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseAll" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",               "TEMP");
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      0 );
        tIQICounter++;

        // Volume mat0
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIVolumeMat0") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial0" );
        tIQICounter++;

        // Volume mat1
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIVolumeMat1") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tIQICounter++;

        // Interface length
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIVolumeInterface1") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases",            "PhaseMaterial0" );
        tIQICounter++;

        // strain energy mat0
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkStrainEnergy0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                    static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion0,Elast");
        tIQICounter++;

        // strain energy mat1
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkStrainEnergy1");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                    static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_constitutive_models", "CMDiffusion1,Elast");
        tIQICounter++;

        // IQI with sideset measure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQISidesetMeasure");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STABILIZATION ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases",            "PhaseMaterial0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPSidesetMeasure,Stabilization");
        tIQICounter++;

        // IQI with sideset length measure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQISidesetLengthMeasure");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STABILIZATION ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases",            "PhaseMaterial0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPSidesetLengthMeasure,Stabilization");
        tIQICounter++;

        // IQI with bulk measure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkMeasure");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STABILIZATION ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPBulkMeasure,Stabilization");
        tIQICounter++;

         // IQI with bulk length measure
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkLengthMeasure");
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STABILIZATION ) );
        tParameterList( tIQIIndex )( tIQICounter ).set( "master_phase_name",          "PhaseMaterial1" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPBulkLengthMeasure,Stabilization");
        tIQICounter++;

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
        tParameterList( tFEMIndex )( 0 ).set( "print_physics_model",                 false );
        // sensitivity analysis
        tParameterList( tFEMIndex )( 0 ).set( "is_analytical_sensitivity",           false );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme",            tFEMWhichFDForSA );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size", stod(tFDEpsilon) );
        // forward analysis
        tParameterList( tFEMIndex )( 0 ).set( "is_analytical_forward",                       true );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme_forward",            tFEMWhichFDForFA );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size_forward", stod(tFDEpsilon) );
        // perturbation strategy
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_strategy",     tFEMWhichPerturbation );
    }

    void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", tNLA_relaxation_parameter  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter",             tNLA_max_iter );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , "TEMP") ;

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",     tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           "TEMP") ;
        tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec", "TEMP,0.0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     "0") ;
        tParameterlist( 5 )( 0 ).set("TSA_Output_Criteria",     "Output_Criterion") ;

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set("order_adofs_by_host",false);
    }

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , tTotalDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names",
                                      "TEMP,VOLUME_MAT1,VOLUME_MAT0,VOLUME_INT1") ;
        tParameterlist( 0 )( 0 ).set( "Field_Type" , "NODAL,GLOBAL,GLOBAL,GLOBAL") ;
        tParameterlist( 0 )( 0 ).set( "IQI_Names"  ,
                                      "IQIBulkTEMP,IQIVolumeMat1,IQIVolumeMat0,IQIVolumeInterface1") ;
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

    }

    /* ------------------------------------------------------------------------ */
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

