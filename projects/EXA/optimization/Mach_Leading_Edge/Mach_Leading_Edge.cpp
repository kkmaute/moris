/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Mach_Leading_Edge.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
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
    /* ------------------------------------------------------------------------ */
    // For Parameter sweep

    // file name
    std::string tName = "Mach_Leading_Edge";

    // Hole Seeding
    moris::sint tNumSeedHolesX    = 3;
    moris::sint tNumSeedHolesY    = 2;
    moris::real tHoleSeedDensityX = 1.0;
    moris::real tHoleSeedDensityY = 0.7;

    // Interpolation order
    std::string tOrder = "1";

    moris::uint tFigureCounter = 0;

    /* ------------------------------------------------------------------------ */
    // For Optimization

    int tNumMaxGcmmaIts = 1;

    // relative step size for MMA
    moris::real tStepSize = 0.1;

    // for volume constraint
    moris::real tInitialFinVolume    = 8.69684e-06;
    moris::real tMaxAllowedFinVolume = 1.0 * tInitialFinVolume;
    moris::real tVolumeScaling       = 1.0;

    // for Perimeter penalty
    moris::real tInitialFinPerimeter   = 0.0288672;
    moris::real tPerimeterPenaltyParam = 1.0;

    // for max temperature objective
    moris::real tRefMaxTemp     = 0.0458574;
    moris::real tMaxTempScaling = 1.0;

    // for strain energy objective
    moris::real tInitialStrainEnergy = 1.74024e-05;
    moris::real tStrainEnergyScaling = 5.0;

    // max dof IQI
    std::string tIQIRefTemp = "350.0";
    std::string tExponent   = "10.0";
    std::string tShift      = "0.0";

    /* ------------------------------------------------------------------------ */
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // WITH PHASE TABLE

    ////Bulk Phases
    std::string tSkin            = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tFins            = "HMR_dummy_n_p3,HMR_dummy_c_p3";
    std::string tShell           = tSkin + "," + tFins;
    std::string tPCM             = "HMR_dummy_n_p4,HMR_dummy_c_p4";
    std::string tBSplineGeometry = tFins;

    // interfaces
    std::string tSkinFinsInterface   = "dbl_iside_p0_2_p1_3";
    std::string tShellShellInterface = tSkinFinsInterface;
    std::string tSkinPCMInterface    = "dbl_iside_p0_2_p1_4";
    std::string tFinPCMInterface     = "dbl_iside_p0_3_p1_4";
    std::string tShellPCMInterface   = tSkinPCMInterface + "," + tFinPCMInterface;

    // boundaries
    std::string tOuterShellSurface = "iside_b0_2_b1_0";
    std::string tSkinBackWall      = "iside_b0_2_b1_1";

    // ghost
    std::string tSkinGhost  = "ghost_p2";
    std::string tFinsGhost  = "ghost_p3";
    std::string tShellGhost = tSkinGhost + "," + tFinsGhost;
    std::string tPCMGhost   = "ghost_p4";

    std::string tTotalDomain       = tShell + "," + tPCM;
    std::string tBSplinesPerimeter = "iside_b0_4_b1_3,iside_b0_4_b1_2,iside_b0_4_b1_1";

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    // all measurements in millimeters
    moris::real mm = 1.0e-03;

    // Shell
    moris::real tLength        = 12.0 * mm;
    moris::real tWallThickness = 2.0 * mm;
    moris::real tInnerRadius   = 1.0 * mm;
    moris::real tEdgeAngle     = 5.0;    // in degrees
    moris::real tCenterX       = 0.0 * mm;
    moris::real tCenterY       = 0.0 * mm;
    moris::real tApproxHeight  = tLength * std::tan( tEdgeAngle / 180.0 * M_PI );

    // Initialize Fins
    moris::sint tNumSeedFinsX = tNumSeedHolesX;
    moris::sint tNumSeedFinsY = tNumSeedHolesY;
    // moris::real tHoleSeedDensity = 1.5
    moris::real tHoleWidth   = tLength / ( tHoleSeedDensityX * tNumSeedFinsX );
    moris::real tHoleHeight  = ( 2.0 * tApproxHeight ) / ( tHoleSeedDensityY * tNumSeedFinsY );
    moris::real tFinExponent = 6.0;

    moris::real tXCenterMin = ( tLength * 1.0 ) / ( (real)tNumSeedFinsX ) - 3.0 * mm;         // + 1.0 );
    moris::real tXCenterMax = ( tLength * (real)tNumSeedFinsX ) / ( (real)tNumSeedFinsX );    // + 1.0 );
    moris::real tYCenterMax = 0.9 * tApproxHeight;
    moris::real tYCenterMin = -1.0 * tYCenterMax;

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // Shell (skin & fins)
    std::string tDensityShell          = "1.0e-2";    // 10,000 kg/m^3
    std::string tHeatCapacityShell     = "200.0";     // 0.2 kJ/kg*K
    std::string tConductivityShell     = "1.5e-4";    // 150 W/m*K
    std::string tYoungsModulusShell    = "2.1e5";     // N/m^2
    std::string tPoissonRatioShell     = "0.31";      //
    std::string tThermalExpansionShell = "1.3e-5";    // 1/K

    // PCM
    std::string tDensityPCM          = "6.0e-3";      // 6,000 kg/m^3
    std::string tHeatCapacityPCM     = "500.0";       // 0.5 kJ/kg*K
    std::string tConductivityPCM     = "1.5e-5";      // 15 W/m*K
    std::string tLatentHeatPCM       = "220000.0";    // 220 kJ/kg
    std::string tPCTempPCM           = "750.0";       // deg C
    std::string tPCConstPCM          = "30.0";        // K
    std::string tYoungsModulusPCM    = "1.0e3";       // N/m^2
    std::string tPoissonRatioPCM     = "0.3";         //
    std::string tThermalExpansionPCM = "0.0";         // 1/K

    /* ------------------------------------------------------------------------ */
    // boundary conditions

    // turn radiation on or off
    bool tHaveRadiation = true;

    // Pressure
    moris::real tAppliedPressure = 200.0;    // in ATMs, modify to match objective gradients
    std::string tPressureDelta   = std::to_string( tAppliedPressure * 1.01325e5 * 1.0e-6 );
    // N/m^2; positive as pulling in direction of normal

    // bedding to suppress RBM
    std::string tBedding = std::to_string( 2.1e5 * 1.0e-5 );

    // Initial Temperature
    moris::real tInitialTemp   = 700.00;    // in deg C
    std::string tReferenceTemp = "20.0";

    // Heat flux
    std::string tHeatLoad = "3.75";    // 375 W/cm^2

    // For Radiation
    std::string tEmissivity   = "5.0e-7";     // 0.5
    std::string tAmbientTemp  = "0.0";        // in deg C
    std::string tAbsoluteZero = "-273.15";    // in deg C

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    // std::string tOrder = "1";

    std::string tNumElemsPerDim = "36, 18";
    std::string tDomainDims     = "0.016, 0.008";
    std::string tDomainOffset   = "-0.004,-0.004";
    std::string tDomainSidesets = "1,2,3,4";

    int         tRefineBuffer         = 2;
    int         tAdaptiveRefineBuffer = 3;
    std::string tInitialRefinement    = "0";
    std::string tAdaptiveRefinements  = "0";

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-05;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 10;

    int         tTSA_Num_Time_Steps = 3;
    moris::real tTSA_Time_Frame     = 0.2;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName  = tName + ".exo";
    std::string tLibraryName     = tName + ".so";
    std::string tHDF5Path        = tName + ".hdf5";
    std::string tGENOutputFile   = "GEN_" + tName + ".exo";
    bool        tOutputCriterion = true;

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // GEOMETRY (LEVEL-SET) FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Outer Wedge
    moris::real Outer_Wedge(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        // get angle in rads
        moris::real tAlpha = tEdgeAngle / 180.0 * M_PI;

        // get outer radius
        moris::real tRadius = tInnerRadius + tWallThickness;

        // translate coordinate system
        moris::real tXi             = aCoordinates( 0 ) - tCenterX;
        moris::real tEta            = aCoordinates( 1 ) - tCenterY;
        moris::real tXiIntersection = tRadius / std::sin( tAlpha );

        // check which sector the point lies in
        bool tLeft   = false;
        bool tRight  = false;
        bool tTop    = false;
        bool tBottom = false;

        if ( tEta >= 0.0 )
            tTop = true;
        else
            tBottom = true;

        if ( ( tTop && tEta < -std::tan( 0.5 * M_PI - tAlpha ) * tXi ) || ( tBottom && tEta > std::tan( 0.5 * M_PI - tAlpha ) * tXi ) )
            tLeft = true;
        else
            tRight = true;

        // Compute Signed-Distance field
        moris::real tVal = 0.0;

        // left sector - circle section
        if ( tLeft )
            tVal = tRadius - std::sqrt( std::pow( tXi, 2.0 ) + std::pow( tEta, 2.0 ) );

        // top right sector - straight section
        else if ( tTop && tRight )
        {
            moris::real tDeltaEta = std::tan( tAlpha ) * ( tXi + tXiIntersection ) - tEta;
            tVal                  = std::cos( tAlpha ) * tDeltaEta;
        }

        // bottom right sector - straight section
        else if ( tBottom && tRight )
        {
            moris::real tDeltaEta = std::tan( tAlpha ) * ( tXi + tXiIntersection ) + tEta;
            tVal                  = std::cos( tAlpha ) * tDeltaEta;
        }

        // clean return value to return non-zero value
        return std::abs( tVal ) < 1.0e-8 ? 1.0e-8 : tVal;
    }

    //-----------------------------------------------------------------------------

    // Inner Wall
    moris::real Inner_Wedge(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        // get angle in rads
        moris::real tAlpha = tEdgeAngle / 180.0 * M_PI;

        // get outer radius
        moris::real tRadius = tInnerRadius;

        // translate coordinate system
        moris::real tXi             = aCoordinates( 0 ) - tCenterX;
        moris::real tEta            = aCoordinates( 1 ) - tCenterY;
        moris::real tXiIntersection = tRadius / std::sin( tAlpha );

        // check which sector the point lies in
        bool tLeft   = false;
        bool tRight  = false;
        bool tTop    = false;
        bool tBottom = false;

        if ( tEta >= 0.0 )
            tTop = true;
        else
            tBottom = true;

        if ( ( tTop && tEta < -std::tan( 0.5 * M_PI - tAlpha ) * tXi ) || ( tBottom && tEta > std::tan( 0.5 * M_PI - tAlpha ) * tXi ) )
            tLeft = true;
        else
            tRight = true;

        // Compute Signed-Distance field
        moris::real tVal = 0.0;

        // left sector - circle section
        if ( tLeft )
            tVal = tRadius - std::sqrt( std::pow( tXi, 2.0 ) + std::pow( tEta, 2.0 ) );

        // top right sector - straight section
        else if ( tTop && tRight )
        {
            moris::real tDeltaEta = std::tan( tAlpha ) * ( tXi + tXiIntersection ) - tEta;
            tVal                  = std::cos( tAlpha ) * tDeltaEta;
        }

        // bottom right sector - straight section
        else if ( tBottom && tRight )
        {
            moris::real tDeltaEta = std::tan( tAlpha ) * ( tXi + tXiIntersection ) + tEta;
            tVal                  = std::cos( tAlpha ) * tDeltaEta;
        }

        // clean return value to return non-zero value
        return std::abs( tVal ) < 1.0e-8 ? 1.0e-8 : tVal;
    }

    //-----------------------------------------------------------------------------

    // Back Wall
    moris::real Back_Wall(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        // compute level set value
        moris::real aReturnValue = ( aCoordinates( 0 ) - tLength + tInnerRadius + tWallThickness );

        // clean return value to return non-zero value
        return std::abs( aReturnValue ) < 1e-8 ? 1e-8 : aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // function for concentrating the heat load at the tip
    void
    Func_Heat_Load_Distribution( moris::Matrix<
                                         moris::DDRMat >&  aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get coordinates
        moris::Matrix< DDRMat > tPosition    = aFIManager->get_IG_geometry_interpolator()->valx();
        moris::real             tEtaPosition = tPosition( 1 ) - tCenterY;
        moris::real             tXiPosition  = tPosition( 0 ) - tCenterX;

        // get surface angle of point on surface
        moris::real tTheta = std::atan2( std::abs( tEtaPosition ), -tXiPosition );
        if ( tTheta > ( 90.0 - tEdgeAngle ) * M_PI / 180.0 )
            tTheta = ( 90.0 - tEdgeAngle ) * M_PI / 180.0;

        // compute heat load distribution
        // help parameter c
        moris::real tC = 1.0 / ( 1.4 * std::pow( 7.0, 2.0 ) );

        // help parameter G
        moris::real tG = ( 1.0 - tC ) * ( std::pow( tTheta, 2.0 ) - 0.5 * tTheta * std::sin( 4.0 * tTheta ) + 0.125 * ( 1.0 - std::cos( 4.0 * tTheta ) ) )
                       + 4.0 * tC * ( std::pow( tTheta, 2.0 ) - tTheta * std::sin( 2.0 * tTheta ) + 0.5 * ( 1.0 - std::cos( 2.0 * tTheta ) ) );

        // compute heat flux relative to stagnation point
        moris::real tR = 2.0 * tTheta * std::sin( tTheta ) * ( ( 1.0 - tC ) * std::pow( std::cos( tTheta ), 2.0 ) + tC ) * std::pow( tG, -0.5 );

        // compute final heat flux using stagnation point heat flux
        aPropMatrix = aParameters( 0 ) * tR;
    }

    // initial temperature
    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { tInitialTemp } };
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    // Dummy function for unused sensitivities if needed
    moris::Matrix< DDRMat > Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // FOR OPTIMIZATION
    /* ------------------------------------------------------------------------ */

    moris::Matrix< moris::DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    moris::Matrix< moris::DDRMat >
    compute_objectives(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tObjectives( 1, 1, 0.0 );

        tObjectives( 0, 0 ) = aCriteria( 3 ) * tStrainEnergyScaling / tInitialStrainEnergy + aCriteria( 0 ) * tMaxTempScaling / tRefMaxTemp + aCriteria( 2 ) * tPerimeterPenaltyParam / tInitialFinPerimeter;

        // std::cout << "% --------------------------------- % \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;
        // std::cout << "Max Temp Value = " << aCriteria( 0 ) << " \n" << std::flush;
        // std::cout << "Fin Volume = "     << aCriteria( 1 ) << " \n" << std::flush;
        // std::cout << "Perimeter = "      << aCriteria( 2 ) << " \n" << std::flush;
        ////std::cout << "Total Volume = "   << aCriteria( 3 ) << " \n" << std::flush;
        // std::cout << "Strain Energy = "  << aCriteria( 3 ) << " \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;

        // std::cout << " \n";
        // std::cout << "% --------------------------------- % \n" << std::flush;

        // std::cout << "min ADV       = " << aADVs.min()         << " \n";
        // std::cout << "max ADV       = " << aADVs.max()         << " \n" << std::flush;

        // std::cout << "Objective = " << tObjectives( 0, 0 ) << " \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;

        //  IQIMaxTemp, IQIStrainEnergy, IQIBSplineGeometryVolume, IQIBSplinesPerimeter, IQITotalVolume

        return tObjectives;
    }

    moris::Matrix< moris::DDRMat >
    compute_constraints(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tConstraints( 1, 1, 0.0 );

        tConstraints( 0, 0 ) = aCriteria( 1 ) * tVolumeScaling / tMaxAllowedFinVolume - 1.0 * tVolumeScaling;

        // std::cout << "Constraint = " << tConstraints( 0, 0 ) << " \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;

        return tConstraints;
    }

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );
        return tDObjectiveDADV;
    }

    //  IQIMaxTemp, IQIStrainEnergy, IQIBSplineGeometryVolume, IQIBSplinesPerimeter, IQITotalVolume

    moris::Matrix< moris::DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDObjectiveDCriteria( 1, 4, 0.0 );

        tDObjectiveDCriteria( 0 ) = tMaxTempScaling / tRefMaxTemp;
        tDObjectiveDCriteria( 2 ) = tPerimeterPenaltyParam / tInitialFinPerimeter;
        tDObjectiveDCriteria( 3 ) = tStrainEnergyScaling / tInitialStrainEnergy;

        // std::cout << "% --------------------------------- % \n" << std::flush;
        // std::cout << "dz/d(Max_Temp) = " << tDObjectiveDCriteria( 0 ) << " \n" << std::flush;
        // std::cout << "dz/d(Perimeter) = " << tDObjectiveDCriteria( 2 ) << " \n" << std::flush;
        ////std::cout << "dz/d(Total_Volume) = " << tDObjectiveDCriteria( 3 ) << " \n" << std::flush;
        // std::cout << "dz/d(Strain_Energy) = " << tDObjectiveDCriteria( 3 ) << " \n" << std::flush;

        return tDObjectiveDCriteria;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );
        return tDConstraintDADV;
    }

    moris::Matrix< moris::DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        moris::Matrix< moris::DDRMat > tDConstraintDCriteria( 1, 4, 0.0 );

        tDConstraintDCriteria( 1 ) = tVolumeScaling / tMaxAllowedFinVolume;

        // std::cout << "dg/d(Fin_Volume) = " << tDConstraintDCriteria( 1 ) << " \n" << std::flush;
        // std::cout << "% --------------------------------- % \n" << std::flush;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", true );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
        aParameterLists( 2 ).set( "max_its", tNumMaxGcmmaIts );
        aParameterLists( 2 ).set( "step_size", tStepSize );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists( 0 ).set( "bspline_orders", tOrder );
        aParameterLists( 0 ).set( "bspline_pattern", std::string( "0" ) );

        //        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        //        aParameterLists( 0 ).set( "refinement_buffer",  tRefineBuffer );
        //        aParameterLists( 0 ).set( "staircase_buffer",   tRefineBuffer );
        //        aParameterLists( 0 ).set( "initial_refinement", tInitialRefinement );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        //        aParameterLists( 0 ).set( "adaptive_refinement_level", tAdaptiveRefineBuffer );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterLists( 0 ).set( "output_path", "./" );
        aParameterLists( 0 ).set( "keep_all_opt_iters", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIMaxTemp", "IQIBSplineGeometryVolume", "IQIBSplinesPerimeter", "IQIStrainEnergy" );    // IQIBSplinesPerimeter
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );

        // phase table
        Matrix< DDUMat > tPhaseMap( 16, 1, 0 );
        tPhaseMap( 8 )              = 2;    // Skin
        tPhaseMap( 9 )              = 2;    // Skin
        tPhaseMap( 10 )             = 1;    // Back wall
        tPhaseMap( 11 )             = 1;    // Back wall
        tPhaseMap( 12 )             = 4;    // PCM
        tPhaseMap( 13 )             = 3;    // Fins
        tPhaseMap( 14 )             = 1;    // Back wall
        tPhaseMap( 15 )             = 1;    // Back wall
        std::string tPhaseMapString = moris::ios::stringify( tPhaseMap );
        aParameterLists( 0 ).set( "print_phase_table", true );
        aParameterLists( 0 ).set( "phase_table", tPhaseMapString );

        // Outer Wall
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Outer_Wedge" );
        // aParameterLists( 1 ).set( "number_of_refinements", tAdaptiveRefinements );

        // Inner Wall
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Inner_Wedge" );
        // aParameterLists( 1 ).set( "number_of_refinements", tAdaptiveRefinements );

        // Back Wall
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Back_Wall" );
        // aParameterLists( 1 ).set( "number_of_refinements", tAdaptiveRefinements );

        // initialize fins as swiss cheese geometry
        aParameterLists( 1 ).add_parameter_list( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        aParameterLists( 1 ).set( "semidiameter_x", 0.5 * tHoleWidth );
        aParameterLists( 1 ).set( "semidiameter_y", 0.5 * tHoleHeight );
        aParameterLists( 1 ).set( "exponent", tFinExponent );
        aParameterLists( 1 ).set( "lower_bound_x", tXCenterMin );           // Left-most hole center
        aParameterLists( 1 ).set( "upper_bound_x", tXCenterMax );           // Right-most hole center
        aParameterLists( 1 ).set( "lower_bound_y", tYCenterMin );           // Bottom-most hole center
        aParameterLists( 1 ).set( "upper_bound_y", tYCenterMax );           // Top-most hole center

        aParameterLists( 1 ).set( "number_of_fields_x", tNumSeedFinsX );    // Number of holes in the x direction
        aParameterLists( 1 ).set( "number_of_fields_y", tNumSeedFinsY );    // Number of holes in the y direction

        aParameterLists( 1 ).set( "discretization_mesh_index", 0 );           // Index of B-spline mesh to create level set field on (-1 = none)
        aParameterLists( 1 ).set( "discretization_lower_bound", -0.0025 );    // Lower bound of level set field (if bspline_mesh_index >= 0)
        aParameterLists( 1 ).set( "discretization_upper_bound", 0.0025 );     // Upper bound of level set field (if bspline_mesh_index >= 0)
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tDensityShell );    // divide by 1000
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Heat Capacity Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatCapacity_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tHeatCapacityShell );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Conductivity Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tConductivityShell );    // divide by  1e6, assume ~90W/mK
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungsModulus_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tYoungsModulusShell );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoissonRatio_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tPoissonRatioShell );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // CTE for Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropThermalExpansion_Shell" );
        aParameterLists( 0 ).set( "function_parameters", tThermalExpansionShell );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - PCM (Al-Cu-Si-alloy?)
        //------------------------------------------------------------------------------

        // Density of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tDensityPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Heat Capacity of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatCapacity_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tHeatCapacityPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Conductivity of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tConductivityPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Latent Heat of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLatentHeat_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tLatentHeatPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Phase Change Temperature of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCTemp_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tPCTempPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Phase Change Temperature Range of PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCconst_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tPCConstPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Cubic Phase State Function for phase change model
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPhaseState_PCM" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Youngs Modulus for PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungsModulus_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tYoungsModulusPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Poisson Ratio for PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoissonRatio_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tPoissonRatioPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // CTE for PCM
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropThermalExpansion_PCM" );
        aParameterLists( 0 ).set( "function_parameters", tThermalExpansionPCM );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // properties for bedding (suppression for RBMs, both Shell and PCM)
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", tBedding );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Dummy latent heat for non-pc material
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLatentHeat_Dummy" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCTemp_Dummy" );
        aParameterLists( 0 ).set( "function_parameters", "10000.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // heat flux from outside
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedFlux" );
        aParameterLists( 0 ).set( "function_parameters", tHeatLoad );
        aParameterLists( 0 ).set( "value_function", "Func_Heat_Load_Distribution" );
        // aParameterLists( 0 ).set( "value_function",           "Func_Const" );

        // reference temperature for thermal expansion
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropReferenceTemp" );
        aParameterLists( 0 ).set( "function_parameters", tReferenceTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // pressure load
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPressure" );
        aParameterLists( 0 ).set( "function_parameters", tPressureDelta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Dirichlet structure
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichletStruct" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // emissivity for radiation BC
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropEmissivity" );
        aParameterLists( 0 ).set( "function_parameters", tEmissivity );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // ambient temperature
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropAmbientTemp" );
        aParameterLists( 0 ).set( "function_parameters", tAmbientTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // absolute zero
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropAbsoluteZero" );
        aParameterLists( 0 ).set( "function_parameters", tAbsoluteZero );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Initial Temperature
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "value_function", "Func_Initial_Condition" );

        //------------------------------------------------------------------------------
        // IQIs
        //------------------------------------------------------------------------------

        // Reference Temperature for MAX_DOF - IQI
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropMaxTempReference" );
        aParameterLists( 0 ).set( "function_parameters", tIQIRefTemp );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Exponent for MAX_DOF - IQI
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropMaxTempExponent" );
        aParameterLists( 0 ).set( "function_parameters", tExponent );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Shift for MAX_DOF - IQI
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropMaxTempShift" );
        aParameterLists( 0 ).set( "function_parameters", tShift );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Reference Temperature for MAX_VMStress - IQI
        // aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        // aParameterLists( 0 ).set( "property_name",            "PropMaxStressReference" );
        // aParameterLists( 0 ).set( "function_parameters",      "1.0" );
        // aParameterLists( 0 ).set( "value_function",           "Func_Const" );
        //
        // Exponent for MAX_VMStress - IQI
        // aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        // aParameterLists( 0 ).set( "property_name",            "PropMaxStressExponent" );
        // aParameterLists( 0 ).set( "function_parameters",      "2.0" );
        // aParameterLists( 0 ).set( "value_function",           "Func_Const" );
        //
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - shell - 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_Shell_1" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity" );

        // create parameter list for constitutive model - shell - 2 (for Nitsche interfaces)
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_Shell_2" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity" );

        // diffusion with phase change - PCM - 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusion_PCM" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity_PCM, Conductivity;"
                "PropDensity_PCM     , Density;"
                "PropHeatCapacity_PCM, HeatCapacity;"
                "PropLatentHeat_PCM  , LatentHeat;"
                "PropPCTemp_PCM      , PCTemp;"
                "PropPhaseState_PCM  , PhaseStateFunction;"
                "PropPCconst_PCM     , PhaseChangeConst" );

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - shell - 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Shell_1" );
        aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungsModulus_Shell,    YoungsModulus;"
                "PropPoissonRatio_Shell,     PoissonRatio;"
                "PropThermalExpansion_Shell, CTE;"
                "PropReferenceTemp,          ReferenceTemperature" );

        // linear elasticity - shell - 2
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Shell_2" );
        aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungsModulus_Shell,    YoungsModulus;"
                "PropPoissonRatio_Shell,     PoissonRatio;"
                "PropThermalExpansion_Shell, CTE;"
                "PropReferenceTemp,          ReferenceTemperature" );

        // linear elasticity - PCM - 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_PCM" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropYoungsModulus_PCM,    YoungsModulus;"
                "PropPoissonRatio_PCM,     PoissonRatio;"
                "PropThermalExpansion_PCM, CTE;"
                "PropReferenceTemp,         ReferenceTemperature" );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // GGLS
        //------------------------------------------------------------------------------

        // create parameter list for GGLS stabilization parameter for Skin
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGGLSDiffusion_Shell" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 2 ).set( "leader_properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity;"
                "PropLatentHeat_Dummy   , LatentHeat;"
                "PropPCTemp_Dummy       , PCTemp;"
                "PropPhaseState_PCM     , PhaseStateFunction;"
                "PropPCconst_PCM        , PhaseChangeConst" );

        // create parameter list for GGLS stabilization parameter for PCM
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGGLSDiffusion_PCM" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 2 ).set( "leader_properties",
                "PropConductivity_PCM , Conductivity;"
                "PropDensity_PCM      , Density;"
                "PropHeatCapacity_PCM , HeatCapacity;"
                "PropLatentHeat_PCM   , LatentHeat;"
                "PropPCTemp_PCM       , PCTemp;"
                "PropPhaseState_PCM   , PhaseStateFunction;"
                "PropPCconst_PCM      , PhaseChangeConst" );

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus_Shell,Material" );

        //------------------------------------------------------------------------------
        // NITSCHE INTERFACE
        //------------------------------------------------------------------------------

        // Temperature - Shell - PCM
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceShellPCMNitsche" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_Shell,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivity_PCM,Material" );

        // Displacements - Shell - PCM - not coupled!

        // Temperature - Shell - Shell
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceShellShellNitsche" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_Shell,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivity_Shell,Material" );

        // Displacements - Shell - Shell
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPInterfaceShellShellNitscheStruct" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus_Shell,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropYoungsModulus_Shell,Material" );

        //------------------------------------------------------------------------------
        // GHOST
        //------------------------------------------------------------------------------

        // bulk Ghost - Shell - Temperature
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp_Shell" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_Shell,Material" );

        // bulk Ghost - PCM - Temperature
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp_PCM" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity_PCM,Material" );

        // bulk Ghost - Shell - Displacements
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPStruct_Shell" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus_Shell,Material" );

        // bulk Ghost - PCM - Displacements
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPStruct_PCM" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus_PCM,Material" );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionShellBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGGLSDiffusion_Shell,GGLSParam" );
        aParameterLists( 3 ).set( "mesh_set_names", tShell );

        // linear elasticity - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGStructShellBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tShell );

        // diffusion - PCM
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionPCMBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_PCM,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGGLSDiffusion_PCM,GGLSParam" );
        aParameterLists( 3 ).set( "mesh_set_names", tPCM );

        // linear elasticity - PCM
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGStructPCBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_PCM,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tPCM );

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInletFlux" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedFlux,Neumann" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterShellSurface );

        // pressure pulling on outside of Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannPressure" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropPressure,Pressure" );
        aParameterLists( 3 ).set( "mesh_set_names", tOuterShellSurface );

        if ( tHaveRadiation )
        {
            // radiation heat flux
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGHeatRadiation" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_RADIATION ) ;
            aParameterLists( 3 ).set( "dof_residual", "TEMP" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            aParameterLists( 3 ).set( "leader_properties",
                    "PropEmissivity,Emissivity;"
                    "PropAmbientTemp,AmbientTemperature;"
                    "PropAbsoluteZero,AbsoluteZero" );
            aParameterLists( 3 ).set( "mesh_set_names", tOuterShellSurface );
            }

        //------------------------------------------------------------------------------
        // INTERFACE BCS - IWGs
        //------------------------------------------------------------------------------

        // Temperature - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceShellShellTEMP" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMDiffusion_Shell_2,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceShellShellNitsche     ,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tShellShellInterface );

        // Displacements - Shell - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceShellShellStruct" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMStrucLinIso_Shell_2,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceShellShellNitscheStruct ,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tShellShellInterface );

        // Temperature - Shell - PCM
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGInterfaceShellPCMTEMP" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists( 3 ).set( "follower_constitutive_models", "CMDiffusion_PCM,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPInterfaceShellPCMNitsche,NitscheInterface" );
        aParameterLists( 3 ).set( "mesh_set_names", tShellPCMInterface );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletStructShell" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tSkinBackWall );

        //------------------------------------------------------------------------------
        // IWGs - GHOST
        //------------------------------------------------------------------------------

        // temperature - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPShellTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp_Shell,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", tShellGhost );

        // displacements - Shell
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPShellStruct" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPStruct_Shell,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", tShellGhost );

        // temperature - PCM
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPPCMTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPTemp_PCM,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", tPCMGhost );

        // displacements - PCM
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGGPPCMStrut" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGPStruct_PCM,GhostSP" );
        aParameterLists( 3 ).set( "mesh_set_names", tPCMGhost );

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", tTotalDomain );
        aParameterLists( 3 ).set( "time_continuity", true );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // Nodal Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Max Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxTemp" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "function_parameters", tIQIRefTemp + "/" + tExponent + "/" + tShift );
        aParameterLists( 4 ).set( "mesh_set_names", tSkin );

        // Strain Energy of Structure
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIStrainEnergy" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tShell );
        // aParameterLists( 4 ).set( "normalization",              "design" ); // <-- what does this do?

        // Volume IQI - Fin Volume - temporary to test optimization capabilities
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBSplineGeometryVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tBSplineGeometry );

        // Volume IQI - TotalDomain - use once to find total volume to compute max dof
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQITotalVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Volume IQI - Fin Perimeter
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBSplinesPerimeter" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 4 ).set( "mesh_set_names", tBSplinesPerimeter );

        // X-displacement
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Y-displacement
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // nodal von-mises stresses for shell
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQINodalVMStress" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VON_MISES_STRESS ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        // aParameterLists( 4 ).set( "vectorial_field_index",      0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // create computation parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        aParameterLists( 0 ).set( "Convergence Tolerance", 1e-12 );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        // ----------------------------------------------------------

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

        // ----------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY,TEMP" );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

        // ----------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 2 );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 3 );

        // ----------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        aParameterLists( 7 ).set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        // aParameterLists( 0 ).set( "Mesh_Type"  ,  vis::VIS_Mesh_Type::STANDARD_WITH_OVERLAP ) ;
        aParameterLists( 0 ).set( "Set_Names", tTotalDomain );
        aParameterLists( 0 ).set( "Field_Names",
                "TEMP,UX,UY,STRESS,"
                "MAX_DOF,VOLUME,VOLUME" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkDISPX,IQIBulkDISPY,IQINodalVMStress,"
                "IQIMaxTemp,IQIBSplineGeometryVolume,IQITotalVolume" );
        aParameterLists( 0 ).set( "Save_Frequency", 5 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
