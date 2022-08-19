
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


#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    //------------------------------------------------------------------------------

    // Main problem parameters

    std::string tName = "Levelset_Boxbeam_Restart";

    bool tIs3D     = false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    bool tUseAbsoulteValue = true;

    // wall thickness
    real tWallThickness = 0.05;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "45,15,15" : "30,10";
    std::string tDimensions        = tIs3D ? "3,1,1" : "3,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";
    std::string tSideSets          = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";

    int tDispOrder = 1;

    // Hole Seeding parameters
    sint tNumHolesY = 3;
    sint tNumHolesZ = 1 * tNumHolesY;
    sint tNumHolesX = 2 * tNumHolesY + 1;

    real tHoleRadius   = tIs3D ? 0.2771 / tNumHolesY : 0.2771 / tNumHolesY;
    real tHoleExponent = 6.0;
    real tHoleScaling  = 1.0;

    real tHoleXdim = 3.0;
    real tHoleYdim = 1.0;
    real tHoleZdim = 1.0;

    real tHoleXOrg = 0.0;
    real tHoleYOrg = 0.0;
    real tHoleZOrg = 0.0;

    bool tHoleOffsetRow = false;

    // Nitsche penalty
    std::string tNitschePenalty = "100.0";

    // optimization parameters
    real tInitialStrainEnergy = tIs3D ? 1.49721 + 1.47892 : 2.04222e+00 + 2.36298e+00;
    real tInitialPerimeter    = tIs3D ? 27.3111 : 14.0771;

    real tInitialRegularization = 0.243945;
    real tPerimeterPenatly      = 0.05;
    real tRegularizationPenatly = 0.01;

    // prescribed gradient when using non-pdv level set field
    std::string tLevelSetGradxConstant = tIs3D ? "1.0;1.0;1.0" : "1.0;1.0";

    real tMaxMass = tIs3D ? 1.5 : 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.05;
    int  tMMAMaxIter  = 14;

    real tBsplineLimit = tHoleRadius;

    // other mesh depedendent parameters
    real tElementEdgeLength = 1.0 / tNumHolesX / pow( 2, 3 );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = "Levelset_Boxbeam_Restart.so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets    = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tVoidSets     = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tTotalDomainSets = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost    = "ghost_p2";
    std::string tInteriorGhost = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets =
            std::string( "SideSet_1_n_p2,SideSet_1_c_p2" ) + std::string( "SideSet_3_n_p2,SideSet_3_c_p2" );

    std::string tInterfaceVoidSSets  = "iside_b0_2_b1_0,iside_b0_1_b1_0";
    std::string tInterfaceVoidSSets2 = "iside_b0_0_b1_1,iside_b0_0_b1_2";

    std::string tInterfaceVoidSSets3 = "iside_b0_1_b1_0";
    std::string tInterfaceVoidSSets4 = "iside_b0_0_b1_1";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tDofStrgAll = tIs3D ? "UX,UY,UZ;THETA;PHID" : "UX,UY;THETA;PHID";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";

    std::string tTotalDomain = tInteriorSets + "," + tVoidSets;

    std::string tTotalDomain1 = tFrameSets + "," + tInteriorSets + "," + tVoidSets;


    /* ------------------------------------------------------------------------ */
    // material parameters

    // conductivity
    std::string tConductivity = "1.0";

    // capacity
    std::string tCapacityTheta = "0.1";
    std::string tCapacityPhi   = "0.0";

    // density
    std::string tDensityTheta = "1.0";
    std::string tDensityPhi   = "0.0";

    // prescribed theta on interface
    std::string tPrescTheta = "1.0";

    // prescribed phi on interface
    std::string tPrescPhi = "0.0";

    /* ------------------------------------------------------------------------ */


    /* ------------------------------------------------------------------------ */

    void
    tLevelSetFuncReal(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val();
    }


    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        // return absolute value of level set function
        aPropMatrix = factor * value;
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of level set function with respect to PDV
    void
    tDerLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->N();
    }

    /* ------------------------------------------------------------------------ */

    // Spatial derivative of level set function defining property in FEM
    void
    tLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->gradx( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of spatial derivative of level set function with respect to PDV
    void
    tDerLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->dnNdxn( 1 );
    }


    //------------------------------------------------------------------------------

    // Hole pattern
    real
    Box_2D3D(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )
    {
        real tBoxExponent = 24.0;

        Matrix< DDRMat > tCenter   = { { 1.5, 0.5, 0.5 } };
        Matrix< DDRMat > tDimOUter = { { 1.5 - tWallThickness, 0.5 - tWallThickness, 0.5 - tWallThickness } };

        real XContrib = std::pow( ( tCenter( 0 ) - aCoordinates( 0 ) ) / tDimOUter( 0 ), tBoxExponent );
        real YContrib = std::pow( ( tCenter( 1 ) - aCoordinates( 1 ) ) / tDimOUter( 1 ), tBoxExponent );

        real tLSval;

        if ( tIs3D )
        {
            real tBoxScaling = pow( tDimOUter( 0 ) * tDimOUter( 1 ) * tDimOUter( 2 ), 1.0 / 3.0 );

            real ZContrib = std::pow( ( tCenter( 2 ) - aCoordinates( 2 ) ) / tDimOUter( 2 ), tHoleExponent );

            tLSval = tBoxScaling * ( 1.0 - std::pow( XContrib + YContrib + ZContrib, 1.0 / tBoxExponent ) );
        }
        else
        {
            real tBoxScaling = pow( tDimOUter( 0 ) * tDimOUter( 1 ), 1.0 / 2.0 );

            tLSval = tBoxScaling * ( 1.0 - std::pow( XContrib + YContrib, 1.0 / tBoxExponent ) );
        }

        // clean return value to return non-zero value
        return -tLSval;
    }

    //------------------------------------------------------------------------------

    // Hole pattern
    real
    Hole_Pattern_2D3D(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )
    {
        Matrix< DDRMat > tDelta  = { { tHoleXdim / tNumHolesX }, { tHoleYdim / tNumHolesY }, { tHoleZdim / tNumHolesZ } };
        Matrix< DDRMat > tOrigin = { { tHoleXOrg + tDelta( 0 ) / 2.0 }, { tHoleYOrg + tDelta( 1, 0 ) / 2.0 }, { tHoleZOrg + tDelta( 2, 0 ) / 2.0 } };
        Matrix< DDRMat > tOffset = { { tDelta( 0, 0 ) / 2.0 }, { 0.0 }, { tDelta( 2, 0 ) / 2.0 } };

        bool tOddRow = true;

        real tLSval = -1e20;
        real dist   = 0.0;

        for ( int iz = -2; iz < tNumHolesZ + 2; ++iz )
        {
            for ( int iy = -2; iy < tNumHolesY + 2; ++iy )
            {
                for ( int ix = -2; ix < tNumHolesX + 2; ++ix )
                {

                    real XCenter = tOrigin( 0 ) + (real)ix * tDelta( 0 ) + tOddRow * tOffset( 0 );
                    real YCenter = tOrigin( 1 ) + (real)iy * tDelta( 1 ) + tOddRow * tOffset( 1 );

                    real XContrib = std::pow( XCenter - aCoordinates( 0 ), tHoleExponent );
                    real YContrib = std::pow( YCenter - aCoordinates( 1 ), tHoleExponent );

                    if ( tIs3D )
                    {
                        real ZCenter  = tOrigin( 2 ) + (real)iz * tDelta( 2 ) + tOddRow * tOffset( 2 );
                        real ZContrib = std::pow( ZCenter - aCoordinates( 2 ), tHoleExponent );

                        dist = tHoleScaling * ( tHoleRadius - std::pow( XContrib + YContrib + ZContrib, 1.0 / tHoleExponent ) );
                    }
                    else
                    {
                        dist = tHoleScaling * ( tHoleRadius - std::pow( XContrib + YContrib, 1.0 / tHoleExponent ) );
                    }
                    tLSval = std::max( tLSval, dist );
                }

                tOddRow = ( !tOddRow && tHoleOffsetRow );
            }
        }

        // clean return value to return non-zero value
        return -tLSval;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.2 )
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
            }
        }
        else
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { 0.0 }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { 0.0 } };
            }
        }
    }

    //------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj3 = tPerimeterPenatly * aCriteria( 3 ) / tInitialPerimeter;
        real obj4 = tRegularizationPenatly * aCriteria( 4 ) / tInitialRegularization;
        // real obj5 = tRegularizationPenatly* aCriteria( 5 ) ;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3 + obj4;    // + obj5;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "H1Error Const            = " << aCriteria( 4 ) << " ( " << aCriteria( 0 ) / tObjectives( 0, 0 ) << " )\n";
        // std::cout << "H1Error PDV              = " << aCriteria( 5 )      << " ( " <<  aCriteria( 1 ) / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        std::cout << "min ADV                  = " << aADVs.min() << " \n";
        std::cout << "max ADV                  = " << aADVs.max() << " \n"
                  << std::flush;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 2 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 2 ) << " \n";
        std::cout << "Constraint = " << tConstraints( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.numel(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.numel(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenatly / tInitialPerimeter;
        tDObjectiveDCriteria( 4 ) = tRegularizationPenatly / tInitialRegularization;
        // tDObjectiveDCriteria( 5 ) = tRegularizationPenatly* 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.numel(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            Matrix< DDRMat > aADVs,
            Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.numel(), 0.0 );

        tDConstraintDCriteria( 2 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }


    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOpt );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
        tParameterlist( 0 )( 0 ).set( "restart_file", "ADV_Alg_0_Iter_11.hdf5" );
        tParameterlist( 0 )( 0 ).set( "reinitialize_interface_iter", 50 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_gcmma_parameter_list();
        tParameterlist( 2 )( 0 ).set( "step_size", tMMAStepSize );
        tParameterlist( 2 )( 0 ).set( "penalty", tMMAPenalty );
        tParameterlist( 2 )( 0 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        tParameterlist( 2 )( 0 ).set( "restart_index", 10 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElementsPerDir );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDimensions );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tOffSet );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tSideSets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "1,1" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1,1" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0,0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1;-1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        // tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        tParameterlist( 0 )( 0 ).set( "use_refine_low_level_elements", true );

        tParameterlist( 0 )( 0 ).set( "restart_refinement_pattern_file", "HMR_Background_Refinement_Iter_10.hdf5" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {

        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 0 )( 0 ) = moris::prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid,IQIHeatMethodPenalty" );
        // tParameterlist(0)(0).set("IQI_types"       , "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid,IQIH1ErrorConst,IQIH1Error");
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );
        tParameterlist( 0 )( 0 ).set( "isocontour_tolerance", 10e-14 );


        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        tParameterlist( 0 )( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        tParameterlist( 0 )( 0 ).set( "print_phase_table", true );

        // init geometry counter
        uint tGeoCounter = 0;

        // outer frame
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Box_2D3D" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "1,0" );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Box" );

        tGeoCounter++;

        // initialize fins as swiss cheese geometry
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Hole_Pattern_2D3D" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Level_Set_Field" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "1,0" );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );

        if ( tIsOpt )
        {
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 1 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -tBsplineLimit );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", tBsplineLimit );
        }
        tGeoCounter++;

        tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list() );
        uint tParamCounter = 0;

        tParameterlist( 2 )( tParamCounter ).set( "name", "LevelsetField" );
        tParameterlist( 2 )( tParamCounter ).set( "dependencies", "Level_Set_Field" );
        tParameterlist( 2 )( tParamCounter ).set( "type", "scaled_field" );
        tParameterlist( 2 )( tParamCounter ).set( "constant_parameters", "1.0" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomain1 );
        tParamCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0e-6" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 5
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropFlux" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDirichletStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropTraction" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Neumann_U" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // common properties for theta and phi problems

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tConductivity );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for Theta

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPrescTheta" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPrescTheta );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties for phi problem

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensityPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCapacityPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tCapacityPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPrescPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPrescPhi );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropEigenStrainPhi" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetConst" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetGradxConst" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tLevelSetGradxConstant );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetReal" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetFuncReal" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;


        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSetGradx" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetGradxFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // time continuity weights
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "10.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // initial condition
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Frame" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Interior" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        // create parameter list for constitutive model - Theta problem
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionTheta" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model - Phi problem
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMDiffusionPhi" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropYoungs,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost_Frame" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost_Interior" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for theta and phi problems
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.01" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropConductivity,Material" );
        tParameterList( 2 )( tSPCounter ).set( "slave_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for DBC on interface for theta problem
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", tNitschePenalty );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG  - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionThetaBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSSets3 );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSSets4 );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        // create IWG - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSSets3 );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSSets4 );
        tIWGCounter++;


        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Frame" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Frame" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameSupportSSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropTraction,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameLoadSSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies", tDofStrgAll );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string( "CMStrucLinIso_Frame,ElastLinIso" ) );
        tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models", std::string( "CMStrucLinIso_Interior,ElastLinIso" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameInteriorDSets );
        tIWGCounter++;

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostFrame" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies", tDofStrgAll );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_Frame,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostInterior" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", tDofStrgAll );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies", tDofStrgAll );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_Interior,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        if ( tIs3D )
        {
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUZ" );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
            tIQICounter++;
        }

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "master_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInteriorSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Frame" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Interior" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInteriorSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInterfaceVoidSSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTHETA" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;

        // Nodal PHID IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkPHID" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;


        // H1 Error if reference is constant
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1ErrorConst" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;

        // H1 Error if reference is design dependent
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;


        // H1 Error if reference is design dependent
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIHeatMethodPenalty" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::HEAT_METHOD_PENALTY ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.6 / 1.0 / 36.0 / 0.2 / 0.2 / 0.333 / 0.333 " );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;

        // Design level set function
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQILevelSet" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropLevelSetReal,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain1 );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 7 );
        for ( uint Ik = 0; Ik < 7; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();


        tParameterlist( 2 ).resize( 4 );
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 2 );

        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 1 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 2 );


        tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 1
        tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 2 );

        tParameterlist( 2 )( 3 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 3 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        tParameterlist( 2 )( 3 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 3 ).set( "NLA_max_iter", 2 );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 4 );
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "THETA" );
        tParameterlist( 3 )( 0 ).set( "NLA_Secundary_DofTypes", "" );

        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 1
        tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "1" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "PHID" );
        tParameterlist( 3 )( 1 ).set( "NLA_Secundary_DofTypes", "" );


        tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 2
        tParameterlist( 3 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        tParameterlist( 3 )( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2" );       // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "THETA;PHID;UX,UY" );


        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 3 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "THETA;PHID;UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Crteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "UX", 0 );
        tParameterlist( 0 )( 0 ).set( "UY", 0 );
        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "UZ", 0 );
        }

        tParameterlist( 0 )( 0 ).set( "THETA", 1 );
        tParameterlist( 0 )( 0 ).set( "PHID", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        // tParameterlist( 0 )( 0 ).set( "Set_Names"  , tTotalDomain );

        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                                                                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" ) );
        }
        else
        {

            tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain1 );
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,PHID,THETA,LEVELSET" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkPHID,IQIBulkTHETA,IQILevelSet" ) );


            // tParameterlist( 0 )( 0 ).set( "Field_Names", std::string("THETA,PHID") );
            // tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string("NODAL,NODAL") );
            // tParameterlist( 0 )( 0 ).set( "IQI_Names"  , std::string("IQIBulkTHETA,IQIBulkPHID") );
        }

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        prm::create_remeshing_parameterlist( tParameterlist( 0 )( 0 ) );
        tParameterlist( 0 )( 0 ).set( "mode", "ab_initio" );
        tParameterlist( 0 )( 0 ).set( "remeshing_field_names", "Box,Level_Set_Field" );
        tParameterlist( 0 )( 0 ).set( "remeshing_levels_of_refinement", "2,1;2,1" );
        tParameterlist( 0 )( 0 ).set( "remeshing_refinement_pattern", "0,1;0,1" );

        tParameterlist( 0 )( 0 ).set( "output_meshes", true );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris


//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
