
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
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"


#ifdef  __cplusplus
extern "C"
{
#endif

//------------------------------------------------------------------------------

namespace moris
{
    //------------------------------------------------------------------------------

    // Main problem parameters

    std::string tName   = "Levelset_Boxbeam_Adaptive_Refinement";

    bool tIs3D     = false;
    bool tIsOpt    = true;
    bool tUseGhost = false;

    // wall thickness
    real tWallThickness = 0.05;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "45,15,15"    : "60,20";
    std::string tDimensions        = tIs3D ? "3,1,1"       : "3,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";
    std::string tSideSets          = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";
    
    int tDispOrder           = 1;

    // Hole Seeding parameters
    sint tNumHolesY = 3;
    sint tNumHolesZ = 1*tNumHolesY;
    sint tNumHolesX = 2*tNumHolesY;

    real tHoleRadius   = tIs3D ? 0.2471/tNumHolesY : 0.2471/tNumHolesY;
    real tHoleExponent = 6.0;
    real tHoleScaling  = 1.0;

    real tHoleXdim  =   3.0;
    real tHoleYdim  =   1.0;
    real tHoleZdim  =   1.0;

    real tHoleXOrg  = 0.0;
    real tHoleYOrg  = 0.0;
    real tHoleZOrg  = 0.0;

    bool tHoleOffsetRow = false;

    // optimization parameters
    real tInitialStrainEnergy = tIs3D ? 1.49721 + 1.47892 : 4.306294e+00 + 4.817279e+00;
    real tInitialPerimeter    = tIs3D ? 27.3111           : 41.5976;
    real tPerimeterPenatly    = 0.2;

    real tMaxMass             = tIs3D ? 1.5 : 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.05;
    int  tMMAMaxIter  = 14;

    real tBsplineLimit  = tHoleRadius;

    // other mesh depedendent parameters
    real tElementEdgeLength = 1.0/tNumHolesX/pow(2,3);
    real tLoadLimitY        = std::floor(0.2/tElementEdgeLength)* tElementEdgeLength;

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = tName+ ".so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets        = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets     = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets  = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost       = "ghost_p2";
    std::string tInteriorGhost    = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets    =
            std::string("SideSet_1_n_p2,SideSet_1_c_p2") +
            std::string("SideSet_3_n_p2,SideSet_3_c_p2");

    std::string tInterfaceVoidSSets = "iside_b0_2_b1_0,iside_b0_1_b1_0";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";


    //------------------------------------------------------------------------------

    // Hole pattern
    real Box_2D3D(
            const Matrix< DDRMat >  & aCoordinates,
            const Cell< real* >     & aGeometryParameters )
    {
        real tBoxExponent = 24.0;

        Matrix< DDRMat > tCenter   = { {1.5, 0.5, 0.5} };
        Matrix< DDRMat > tDimOUter = { {1.5-tWallThickness, 0.5-tWallThickness, 0.5-tWallThickness} };

        real XContrib = std::pow( (tCenter(0) - aCoordinates(0))/tDimOUter(0), tBoxExponent);
        real YContrib = std::pow( (tCenter(1) - aCoordinates(1))/tDimOUter(1), tBoxExponent);

        real tLSval;

        if (tIs3D)
        {
            real tBoxScaling = pow(tDimOUter(0)*tDimOUter(1)*tDimOUter(2),1.0/3.0);

            real ZContrib = std::pow( (tCenter(2) - aCoordinates(2))/tDimOUter(2), tHoleExponent);

            tLSval = tBoxScaling * ( 1.0 - std::pow (XContrib + YContrib + ZContrib, 1.0/tBoxExponent) );
        }
        else
        {
            real tBoxScaling = pow(tDimOUter(0)*tDimOUter(1),1.0/2.0);

            tLSval = tBoxScaling * ( 1.0 - std::pow (XContrib + YContrib, 1.0/tBoxExponent) );
        }

        // clean return value to return non-zero value
        return std::abs(tLSval) < 1e-8 ? 1e-8 : -tLSval;
    }

    //------------------------------------------------------------------------------

    // Hole pattern
    real Hole_Pattern_2D3D(
            const Matrix< DDRMat >  & aCoordinates,
            const Cell< real* >     & aGeometryParameters )
    {
        Matrix< DDRMat > tDelta  = { {tHoleXdim/tNumHolesX}   ,{tHoleYdim/tNumHolesY}     ,{tHoleZdim/tNumHolesZ} };
        Matrix< DDRMat > tOrigin = { {tHoleXOrg+tDelta(0)/2.0},{tHoleYOrg+tDelta(1,0)/2.0},{tHoleZOrg+tDelta(2,0)/2.0} } ;
        Matrix< DDRMat > tOffset = { {tDelta(0,0)/2.0}        ,{0.0}                      ,{tDelta(2,0)/2.0} };

        bool tOddRow    = true;

        real tLSval = -1e20;
        real dist   = 0.0;

        for (int iz=-2;iz<tNumHolesZ+2;++iz)
        {
            for (int iy=-2;iy<tNumHolesY+2;++iy)
            {
                for (int ix=-2;ix<tNumHolesX+2;++ix)
                {

                    real XCenter = tOrigin(0) + (real)ix * tDelta(0) + tOddRow * tOffset(0);
                    real YCenter = tOrigin(1) + (real)iy * tDelta(1) + tOddRow * tOffset(1);

                    real XContrib = std::pow( XCenter - aCoordinates(0), tHoleExponent);
                    real YContrib = std::pow( YCenter - aCoordinates(1), tHoleExponent);

                    if (tIs3D)
                    {
                        real ZCenter  = tOrigin(2) + (real)iz * tDelta(2) + tOddRow * tOffset(2);
                        real ZContrib = std::pow( ZCenter - aCoordinates(2), tHoleExponent);

                        dist = tHoleScaling * ( tHoleRadius - std::pow (XContrib + YContrib + ZContrib, 1.0/tHoleExponent) );
                    }
                    else
                    {
                        dist = tHoleScaling * ( tHoleRadius - std::pow (XContrib + YContrib, 1.0/tHoleExponent) );
                    }
                    tLSval = std::max(tLSval,dist);
                }

                tOddRow = (! tOddRow && tHoleOffsetRow );
            }
        }

        // clean return value to return non-zero value
        return std::abs(tLSval) < 1e-8 ? 1e-8 : -tLSval;
    }

   //--------------------------------------------------------------------------------------------------------------

    void Func_Neumann_U(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        if( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.2 )
        {
            if ( tIs3D )
            {
                aPropMatrix = {{ 0.0 },{ aParameters( 0 )( 0 ) },{ 0.0 }};
            }
            else
            {
                aPropMatrix = {{ 0.0 },{ aParameters( 0 )( 0 ) }};
            }
        }
        else
        {
            if ( tIs3D )
            {
                aPropMatrix = {{ 0.0 },{ 0.0 },{ 0.0 }};
            }
            else
            {
                aPropMatrix =  {{ 0.0 },{ 0.0 }};
            }
        }
    }

    //------------------------------------------------------------------------------

    // Constant function for properties
    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

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

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj3 = tPerimeterPenatly * aCriteria( 3 ) / tInitialPerimeter;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 )      << " ( " <<  obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 )      << " ( " <<  obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 )      << " ( " <<  obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        std::cout << "min ADV                  = " << aADVs.min()         << " \n";
        std::cout << "max ADV                  = " << aADVs.max()         << " \n" << std::flush;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_constraints(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tConstraints( 1, 1 );
        tConstraints(0) = aCriteria( 2 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 2 )    << " \n";
        std::cout << "Constraint = " << tConstraints( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n" << std::flush;

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dobjective_dadv(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dobjective_dcriteria(
            Matrix<DDRMat> aADVs,
            Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenatly / tInitialPerimeter;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dconstraint_dadv(
            Matrix<DDRMat> aADVs,
            Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix<DDRMat> compute_dconstraint_dcriteria(
            Matrix<DDRMat> aADVs,
            Matrix<DDRMat> aCriteria)
    {
        Matrix<DDRMat> tDConstraintDCriteria( 1, aCriteria.numel(), 0.0 );

        tDConstraintDCriteria( 2 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }


    //--------------------------------------------------------------------------------------------------------------

    void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist(0)(0) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist(0)(0).set("is_optimization_problem", tIsOpt);
        tParameterlist(0)(0).set("problem", "user_defined");
        tParameterlist(0)(0).set( "library", tLibraryName );
        tParameterlist(0)(0).set( "restart_file", "" );

        tParameterlist(2)(0) = moris::prm::create_gcmma_parameter_list();
        tParameterlist(2)(0).set("step_size", tMMAStepSize);
        tParameterlist(2)(0).set("penalty"  , tMMAPenalty );
        tParameterlist(2)(0).set("max_its"  , tMMAMaxIter );   // Maximum number of iterations
        tParameterlist(2)(0).set("restart_index", 0);
    }

    //--------------------------------------------------------------------------------------------------------------

    void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElementsPerDir);
        tParameterlist( 0 )( 0 ).set( "domain_dimensions",                tDimensions    );
        tParameterlist( 0 )( 0 ).set( "domain_offset",                    tOffSet);
        tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  tSideSets);
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 )( 0 ).set( "lagrange_orders",  "1"   );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "bspline_orders",   "1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern",  "0" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement",         "0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0");

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer",  1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer",   1 );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1);

        tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

    }

    //--------------------------------------------------------------------------------------------------------------

    void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose",                 true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type",        "conformal");
        tParameterlist( 0 )( 0 ).set( "enrich",                    true );
        tParameterlist( 0 )( 0 ).set( "basis_rank",                "bspline");
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       "0");
        tParameterlist( 0 )( 0 ).set( "ghost_stab",                tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 )( 0 ).set( "verbose",                   true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {

        tParameterlist.resize(3);
        tParameterlist(0).resize(1);
        tParameterlist(0)(0) = moris::prm::create_gen_parameter_list();
        tParameterlist(0)(0).set("IQI_types"       , "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid");
        tParameterlist(0)(0).set("output_mesh_file", tGENOutputFile );
        tParameterlist(0)(0).set("time_offset"     , 10.0 );

        Matrix<DDUMat> tPhaseMap(4,1,0);
        tPhaseMap(1) = 1;
        tPhaseMap(2) = 2;
        tPhaseMap(3) = 2;
        tParameterlist( 0 )( 0 ).set("phase_table",moris::ios::stringify(tPhaseMap));

        tParameterlist( 0 )( 0 ).set("print_phase_table", true);

       // init geometry counter
        uint tGeoCounter = 0;

        // outer frame
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Box_2D3D" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Box" );

        tGeoCounter++;

        // initialize fins as swiss cheese geometry
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Hole_Pattern_2D3D" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Level_Set_Field" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "1" );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );

        if (tIsOpt)
        {
            tParameterlist( 1 )( tGeoCounter ).set("discretization_mesh_index",   0);
            tParameterlist( 1 )( tGeoCounter ).set("discretization_lower_bound", -tBsplineLimit);
            tParameterlist( 1 )( tGeoCounter ).set("discretization_upper_bound",  tBsplineLimit);
        }
        tGeoCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropYoungs");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropBedding");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0e-6");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropFlux");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "10.0");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDirichletU");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      tDirichletStr);
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropTraction");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Neumann_U");
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropPoisson");
        tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0");
        tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const");
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Frame");
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio");
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Interior");
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",        "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio");
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheDirichletBC");
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "100.0");
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropYoungs,Material");
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPNitscheFrameInteriorInterface") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropYoungs,Material") );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropYoungs,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGhost_Frame") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.005") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropYoungs,Material") );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGhost_Interior") );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.005") );
        tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropYoungs,Material") );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGBulkU_Frame");
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso");
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropBedding,Bedding");
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrameSets);
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGBulkU_Frame");
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso");
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropBedding,Bedding");
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInteriorSets);
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGDirichletU");
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropDirichletU,Dirichlet");
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso");
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheDirichletBC,DirichletNitsche");
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrameSupportSSets);
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGTraction");
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropTraction,Traction");
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrameLoadSSets);
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGFrameInteriorInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMStrucLinIso_Frame,ElastLinIso") );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",  std::string("CMStrucLinIso_Interior,ElastLinIso") );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheFrameInteriorInterface,NitscheInterface") );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrameInteriorDSets );
        tIWGCounter++;

        if (tUseGhost)
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGhostFrame") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGhost_Frame,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order",                (uint) tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tFrameGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGhostInterior") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGhost_Interior,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order",                (uint) tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             tInteriorGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkUX");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );      
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",       0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomainSets);
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkUY");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",       1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomainSets);
        tIQICounter++;

        if (tIs3D)
        {
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkUZ");
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               tDofStrg);
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",       2 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tTotalDomainSets);
            tIQICounter++;
        }

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkStrainEnergy_Frame");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,Elast");
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tFrameSets);
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkStrainEnergy_Interior");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", "CMStrucLinIso_Interior,Elast");
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tInteriorSets);
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVolume_Frame");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ));
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          "PropDensity,Density");
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tFrameSets);
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVolume_Interior");
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ));
        tParameterList( 4 )( tIQICounter ).set( "master_properties",          "PropDensity,Density");
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tInteriorSets);
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIPerimeter_InterfaceVoid" ) ;
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    tDofStrg);
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             tInterfaceVoidSSets);
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    //--------------------------------------------------------------------------------------------------------------

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
        tParameterlist( 2 )( 0 ).set("NLA_combined_res_jac_assembly", true );
        tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1.00 );
        tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 1.00  );
        tParameterlist( 2 )( 0 ).set("NLA_max_iter", 1 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , tDofStrg);

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set("TSA_DofTypes",           tDofStrg);
        tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",     "0");
        tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",     "Output_Criterion");

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    //--------------------------------------------------------------------------------------------------------------

    void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set("UX" , 0);
        tParameterlist( 0 )( 0 ).set("UY" , 0);
        if (tIs3D)
        {
            tParameterlist( 0 )( 0 ).set("UZ" , 0);
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names"  , tTotalDomainSets );

        if (tIs3D)
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string("UX,UY,UZ,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid") );
            tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string("NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL") );
            tParameterlist( 0 )( 0 ).set( "IQI_Names"  , std::string("IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid") );
        }
        else
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string("UX,UY,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid") );
            tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string("NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL") );
            tParameterlist( 0 )( 0 ).set( "IQI_Names"  , std::string("IQIBulkUX,IQIBulkUY,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid") );
        }

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset"   , 10.0 );

    }

    //--------------------------------------------------------------------------------------------------------------

    void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );


         prm::create_remeshing_parameterlist( tParameterlist( 0 )( 0 ) );
         tParameterlist( 0 )( 0 ).set( "mode" , "ab_initio" );
         tParameterlist( 0 )( 0 ).set( "remeshing_field_names" , "Box,Level_Set_Field" );
         tParameterlist( 0 )( 0 ).set( "remeshing_levels_of_refinement" , "1" );
         tParameterlist( 0 )( 0 ).set( "remeshing_refinement_pattern" , "0" );
    }

    //--------------------------------------------------------------------------------------------------------------
}


//--------------------------------------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
