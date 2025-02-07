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

    Vector< uint > tNumElemsPerDim = { 36, 18 };
    Vector< real > tDomainDims     = { 0.016, 0.008 };
    Vector< real > tDomainOffset   = { -0.004, -0.004 };

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
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "max_its", tNumMaxGcmmaIts );
        aParameterLists.set( "step_size", tStepSize );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tOrder );
        aParameterLists.set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists.set( "bspline_orders", tOrder );
        aParameterLists.set( "bspline_pattern", std::string( "0" ) );

        //        aParameterLists.set( "lagrange_to_bspline", "0" );

        //        aParameterLists.set( "refinement_buffer",  tRefineBuffer );
        //        aParameterLists.set( "staircase_buffer",   tRefineBuffer );
        //        aParameterLists.set( "initial_refinement", tInitialRefinement );
        //        aParameterLists.set( "adaptive_refinement_level", tAdaptiveRefineBuffer );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "output_path", "./" );
        aParameterLists.set( "keep_all_opt_iters", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIMaxTemp", "IQIBSplineGeometryVolume", "IQIBSplinesPerimeter", "IQIStrainEnergy" );    // IQIBSplinesPerimeter
        aParameterLists.set( "output_mesh_file", tGENOutputFile );

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
        aParameterLists.set( "print_phase_table", true );
        aParameterLists.set( "phase_table", tPhaseMapString );

        // Outer Wall
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Outer_Wedge" );
        // aParameterLists.set( "number_of_refinements", tAdaptiveRefinements );

        // Inner Wall
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Inner_Wedge" );
        // aParameterLists.set( "number_of_refinements", tAdaptiveRefinements );

        // Back Wall
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Back_Wall" );
        // aParameterLists.set( "number_of_refinements", tAdaptiveRefinements );

        // initialize fins as swiss cheese geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
        aParameterLists.set( "semidiameter_x", 0.5 * tHoleWidth );
        aParameterLists.set( "semidiameter_y", 0.5 * tHoleHeight );
        aParameterLists.set( "exponent", tFinExponent );
        aParameterLists.set( "lower_bound_x", tXCenterMin );           // Left-most hole center
        aParameterLists.set( "upper_bound_x", tXCenterMax );           // Right-most hole center
        aParameterLists.set( "lower_bound_y", tYCenterMin );           // Bottom-most hole center
        aParameterLists.set( "upper_bound_y", tYCenterMax );           // Top-most hole center

        aParameterLists.set( "number_of_fields_x", tNumSeedFinsX );    // Number of holes in the x direction
        aParameterLists.set( "number_of_fields_y", tNumSeedFinsY );    // Number of holes in the y direction

        aParameterLists.set( "discretization_mesh_index", 0 );           // Index of B-spline mesh to create level set field on (-1 = none)
        aParameterLists.set( "discretization_lower_bound", -0.0025 );    // Lower bound of level set field (if bspline_mesh_index >= 0)
        aParameterLists.set( "discretization_upper_bound", 0.0025 );     // Upper bound of level set field (if bspline_mesh_index >= 0)
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
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity_Shell" );
        aParameterLists.set( "function_parameters", tDensityShell );    // divide by 1000
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat Capacity Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatCapacity_Shell" );
        aParameterLists.set( "function_parameters", tHeatCapacityShell );
        aParameterLists.set( "value_function", "Func_Const" );

        // Conductivity Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity_Shell" );
        aParameterLists.set( "function_parameters", tConductivityShell );    // divide by  1e6, assume ~90W/mK
        aParameterLists.set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsModulus_Shell" );
        aParameterLists.set( "function_parameters", tYoungsModulusShell );
        aParameterLists.set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoissonRatio_Shell" );
        aParameterLists.set( "function_parameters", tPoissonRatioShell );
        aParameterLists.set( "value_function", "Func_Const" );

        // CTE for Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropThermalExpansion_Shell" );
        aParameterLists.set( "function_parameters", tThermalExpansionShell );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - PCM (Al-Cu-Si-alloy?)
        //------------------------------------------------------------------------------

        // Density of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity_PCM" );
        aParameterLists.set( "function_parameters", tDensityPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat Capacity of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatCapacity_PCM" );
        aParameterLists.set( "function_parameters", tHeatCapacityPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Conductivity of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity_PCM" );
        aParameterLists.set( "function_parameters", tConductivityPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Latent Heat of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLatentHeat_PCM" );
        aParameterLists.set( "function_parameters", tLatentHeatPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Phase Change Temperature of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCTemp_PCM" );
        aParameterLists.set( "function_parameters", tPCTempPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Phase Change Temperature Range of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCconst_PCM" );
        aParameterLists.set( "function_parameters", tPCConstPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Cubic Phase State Function for phase change model
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPhaseState_PCM" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Youngs Modulus for PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsModulus_PCM" );
        aParameterLists.set( "function_parameters", tYoungsModulusPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // Poisson Ratio for PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoissonRatio_PCM" );
        aParameterLists.set( "function_parameters", tPoissonRatioPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        // CTE for PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropThermalExpansion_PCM" );
        aParameterLists.set( "function_parameters", tThermalExpansionPCM );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // properties for bedding (suppression for RBMs, both Shell and PCM)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", tBedding );
        aParameterLists.set( "value_function", "Func_Const" );

        // Dummy latent heat for non-pc material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLatentHeat_Dummy" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPCTemp_Dummy" );
        aParameterLists.set( "function_parameters", "10000.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // heat flux from outside
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImposedFlux" );
        aParameterLists.set( "function_parameters", tHeatLoad );
        aParameterLists.set( "value_function", "Func_Heat_Load_Distribution" );
        // aParameterLists.set( "value_function",           "Func_Const" );

        // reference temperature for thermal expansion
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropReferenceTemp" );
        aParameterLists.set( "function_parameters", tReferenceTemp );
        aParameterLists.set( "value_function", "Func_Const" );

        // pressure load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPressure" );
        aParameterLists.set( "function_parameters", tPressureDelta );
        aParameterLists.set( "value_function", "Func_Const" );

        // Dirichlet structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletStruct" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // emissivity for radiation BC
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropEmissivity" );
        aParameterLists.set( "function_parameters", tEmissivity );
        aParameterLists.set( "value_function", "Func_Const" );

        // ambient temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropAmbientTemp" );
        aParameterLists.set( "function_parameters", tAmbientTemp );
        aParameterLists.set( "value_function", "Func_Const" );

        // absolute zero
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropAbsoluteZero" );
        aParameterLists.set( "function_parameters", tAbsoluteZero );
        aParameterLists.set( "value_function", "Func_Const" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Initial Temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "value_function", "Func_Initial_Condition" );

        //------------------------------------------------------------------------------
        // IQIs
        //------------------------------------------------------------------------------

        // Reference Temperature for MAX_DOF - IQI
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropMaxTempReference" );
        aParameterLists.set( "function_parameters", tIQIRefTemp );
        aParameterLists.set( "value_function", "Func_Const" );

        // Exponent for MAX_DOF - IQI
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropMaxTempExponent" );
        aParameterLists.set( "function_parameters", tExponent );
        aParameterLists.set( "value_function", "Func_Const" );

        // Shift for MAX_DOF - IQI
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropMaxTempShift" );
        aParameterLists.set( "function_parameters", tShift );
        aParameterLists.set( "value_function", "Func_Const" );

        // Reference Temperature for MAX_VMStress - IQI
        // aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        // aParameterLists.set( "property_name",            "PropMaxStressReference" );
        // aParameterLists.set( "function_parameters",      "1.0" );
        // aParameterLists.set( "value_function",           "Func_Const" );
        //
        // Exponent for MAX_VMStress - IQI
        // aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        // aParameterLists.set( "property_name",            "PropMaxStressExponent" );
        // aParameterLists.set( "function_parameters",      "2.0" );
        // aParameterLists.set( "value_function",           "Func_Const" );
        //
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - shell - 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_Shell_1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity" );

        // create parameter list for constitutive model - shell - 2 (for Nitsche interfaces)
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_Shell_2" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity" );

        // diffusion with phase change - PCM - 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion_PCM" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
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
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Shell_1" );
        aParameterLists.set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists.set( "properties",
                "PropYoungsModulus_Shell,    YoungsModulus;"
                "PropPoissonRatio_Shell,     PoissonRatio;"
                "PropThermalExpansion_Shell, CTE;"
                "PropReferenceTemp,          ReferenceTemperature" );

        // linear elasticity - shell - 2
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Shell_2" );
        aParameterLists.set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists.set( "properties",
                "PropYoungsModulus_Shell,    YoungsModulus;"
                "PropPoissonRatio_Shell,     PoissonRatio;"
                "PropThermalExpansion_Shell, CTE;"
                "PropReferenceTemp,          ReferenceTemperature" );

        // linear elasticity - PCM - 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_PCM" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
        aParameterLists.set( "properties",
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
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGGLSDiffusion_Shell" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "leader_properties",
                "PropConductivity_Shell , Conductivity;"
                "PropDensity_Shell      , Density;"
                "PropHeatCapacity_Shell , HeatCapacity;"
                "PropLatentHeat_Dummy   , LatentHeat;"
                "PropPCTemp_Dummy       , PCTemp;"
                "PropPhaseState_PCM     , PhaseStateFunction;"
                "PropPCconst_PCM        , PhaseChangeConst" );

        // create parameter list for GGLS stabilization parameter for PCM
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGGLSDiffusion_PCM" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "leader_properties",
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
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungsModulus_Shell,Material" );

        //------------------------------------------------------------------------------
        // NITSCHE INTERFACE
        //------------------------------------------------------------------------------

        // Temperature - Shell - PCM
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceShellPCMNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity_Shell,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity_PCM,Material" );

        // Displacements - Shell - PCM - not coupled!

        // Temperature - Shell - Shell
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceShellShellNitsche" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity_Shell,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity_Shell,Material" );

        // Displacements - Shell - Shell
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceShellShellNitscheStruct" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungsModulus_Shell,Material" );
        aParameterLists.set( "follower_properties", "PropYoungsModulus_Shell,Material" );

        //------------------------------------------------------------------------------
        // GHOST
        //------------------------------------------------------------------------------

        // bulk Ghost - Shell - Temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp_Shell" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity_Shell,Material" );

        // bulk Ghost - PCM - Temperature
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp_PCM" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity_PCM,Material" );

        // bulk Ghost - Shell - Displacements
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPStruct_Shell" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungsModulus_Shell,Material" );

        // bulk Ghost - PCM - Displacements
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPStruct_PCM" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungsModulus_PCM,Material" );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionShellBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPGGLSDiffusion_Shell,GGLSParam" );
        aParameterLists.set( "mesh_set_names", tShell );

        // linear elasticity - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGStructShellBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tShell );

        // diffusion - PCM
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionPCMBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_PCM,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPGGLSDiffusion_PCM,GGLSParam" );
        aParameterLists.set( "mesh_set_names", tPCM );

        // linear elasticity - PCM
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGStructPCBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_PCM,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tPCM );

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletFlux" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedFlux,Neumann" );
        aParameterLists.set( "mesh_set_names", tOuterShellSurface );

        // pressure pulling on outside of Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGNeumannPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_properties", "PropPressure,Pressure" );
        aParameterLists.set( "mesh_set_names", tOuterShellSurface );

        if ( tHaveRadiation )
        {
            // radiation heat flux
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGHeatRadiation" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_RADIATION ) ;
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            aParameterLists.set( "leader_properties",
                    "PropEmissivity,Emissivity;"
                    "PropAmbientTemp,AmbientTemperature;"
                    "PropAbsoluteZero,AbsoluteZero" );
            aParameterLists.set( "mesh_set_names", tOuterShellSurface );
            }

        //------------------------------------------------------------------------------
        // INTERFACE BCS - IWGs
        //------------------------------------------------------------------------------

        // Temperature - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceShellShellTEMP" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusion_Shell_2,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceShellShellNitsche     ,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tShellShellInterface );

        // Displacements - Shell - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceShellShellStruct" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIso_Shell_2,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceShellShellNitscheStruct ,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tShellShellInterface );

        // Temperature - Shell - PCM
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterfaceShellPCMTEMP" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion_Shell_1,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusion_PCM,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceShellPCMNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tShellPCMInterface );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletStructShell" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tSkinBackWall );

        //------------------------------------------------------------------------------
        // IWGs - GHOST
        //------------------------------------------------------------------------------

        // temperature - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPShellTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPTemp_Shell,GhostSP" );
        aParameterLists.set( "mesh_set_names", tShellGhost );

        // displacements - Shell
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPShellStruct" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "stabilization_parameters", "SPGPStruct_Shell,GhostSP" );
        aParameterLists.set( "mesh_set_names", tShellGhost );

        // temperature - PCM
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPPCMTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "follower_dof_dependencies", "TEMP" );
        aParameterLists.set( "stabilization_parameters", "SPGPTemp_PCM,GhostSP" );
        aParameterLists.set( "mesh_set_names", tPCMGhost );

        // displacements - PCM
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPPCMStrut" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "stabilization_parameters", "SPGPStruct_PCM,GhostSP" );
        aParameterLists.set( "mesh_set_names", tPCMGhost );

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );
        aParameterLists.set( "time_continuity", true );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Max Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxTemp" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aParameterLists.set( "function_parameters", tIQIRefTemp + "/" + tExponent + "/" + tShift );
        aParameterLists.set( "mesh_set_names", tSkin );

        // Strain Energy of Structure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIStrainEnergy" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,Elast" );
        aParameterLists.set( "mesh_set_names", tShell );
        // aParameterLists.set( "normalization",              "design" ); // <-- what does this do?

        // Volume IQI - Fin Volume - temporary to test optimization capabilities
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBSplineGeometryVolume" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "mesh_set_names", tBSplineGeometry );

        // Volume IQI - TotalDomain - use once to find total volume to compute max dof
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITotalVolume" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Volume IQI - Fin Perimeter
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBSplinesPerimeter" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "mesh_set_names", tBSplinesPerimeter );

        // X-displacement
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Y-displacement
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQINodalVMStress" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VON_MISES_STRESS ) ;
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Shell_1,ElastLinIso" );
        // aParameterLists.set( "vectorial_field_index",      0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create computation parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
        aParameterLists.set( "Convergence Tolerance", 1e-12 );
        aParameterLists.set( "preconditioners", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        // ----------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        // ----------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "1,0" );
        aParameterLists.set( "NLA_DofTypes", "UX,UY;TEMP" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY,TEMP" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );
        aParameterLists.set( "TSA_Nonlinear_Solver", 2 );
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 3 );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY;TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        aParameterLists.set( "ifpack_prec_type", "ILU" );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Time_Offset", 10.0 );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        // aParameterLists.set( "Mesh_Type"  ,  vis::VIS_Mesh_Type::STANDARD_WITH_OVERLAP ) ;
        aParameterLists.set( "Set_Names", tTotalDomain );
        aParameterLists.set( "Field_Names",
                "TEMP,UX,UY,STRESS,"
                "MAX_DOF,VOLUME,VOLUME" );
        aParameterLists.set( "Field_Type",
                "NODAL,NODAL,NODAL,NODAL,"
                "GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTEMP,IQIBulkDISPX,IQIBulkDISPY,IQINodalVMStress,"
                "IQIMaxTemp,IQIBSplineGeometryVolume,IQITotalVolume" );
        aParameterLists.set( "Save_Frequency", 5 );
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
