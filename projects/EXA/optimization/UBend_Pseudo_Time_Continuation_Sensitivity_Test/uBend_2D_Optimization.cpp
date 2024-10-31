#include <string>
#include <iostream>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "GEN_Data_Types.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // interpolation order
    uint gInterpolationOrder = 1;

    std::string
    moris_to_string( real tValue )
    {
        std::ostringstream streamObj;

        // Set precision
        streamObj << std::scientific;
        streamObj << std::setprecision( 15 );

        // Add value to stream
        streamObj << tValue;

        // Get string from output string stream
        return streamObj.str();
    }

    /* ------------------------------------------------------------------------ */
    // set optimization restart iteration (if 0 no restart)
    sint tRestartId = -1;

    bool tIsOpt              = true;
    bool tCheckSensitivities = true;

    uint tNumConstraints = tCheckSensitivities ? 3 : 2;

    /* ------------------------------------------------------------------------ */
    // basic material properties of fluid (air)
    real tDensFld = 1.0;       // kg/m3
    real tKinVisc = 4.0e-5;    // m2/s

    // inlet bc and flux bc
    real tInVel         = 1.0e-2;    // m/s // modified (Re=62.5 -> tInVel=2.5e-2); original (Re=5000 -> tInVel=2.0)
    real tTurbIntensity = -3.0;      // if negative: ratio wrt. tKinVisc

    /* ------------------------------------------------------------------------ */

    // reference values
    real velref = tInVel;      // reference velocity
    real lenref = 0.1;         // reference length
    real rhoref = tDensFld;    // reference density

    std::string tCMTurbFt2      = "0.0";
    std::string tCMTurbAlpha    = "10.0";
    std::string tSupgTurbPower  = "2.0";
    std::string tSupgTurbSource = "1.0";
    std::string tSupgFluidC1    = "36.0";

    uint tNumElements = 2;

    bool tPowerlawProfile = false;    // FIXME to be changed to power law for turbulence profile

    bool tRampViscBC = false;

    /* ------------------------------------------------------------------------ */
    // Solver config
    //
    // The following configurations have been tested
    //
    // Case     tSolverType    tFluidViscSolver  tSubNewtonSolverOption  tPropRampSteps  tUsePseudoTimeStepping
    // 0        NLBGS          Mono              0                       3               No                      ok
    // 1        NLBGS          Mono              5                       3               No                      ok
    // 2        NLBGS          Mono              5                       3               yes                     ok
    // 3        NLBGS          Staggered         5                       3               No                      ok
    // 4        NLBGS          Staggered         5                       3               Yes                     ok  (configuration tested via ctest)
    // 5        Newton         Mono              0                       3               na                      ok

    // Solution strategy for turbulent flow subproblem
    std::string tSolverType            = "NLBGS";        // "Newton" or "NLBGS";
    std::string tFluidViscSolver       = "Staggered";    // "Staggered" or "Mono"; only monolithic scheme possible with "Newton" algorithm
    std::string tSubNewtonSolverOption = "5";            // "0": standard Newton solver  "5": Newton with adaptive relaxation scheme

    bool tUsePseudoTimeStepping = true;    // NLBGS with pseudo time stepping

    // Ramping of flow properties
    sint tPropRampSteps   = 3;        // number of NLBFGS steps for ramping properties (>=1); here dynamic viscosity; nominal value reached in iteration tPropRampSteps+1
    real tPropRampScaling = 100.0;    // factor by which dynamic viscosity is increased

    // NLBGS settings for fluid subproblem
    sint tNLBGS_max_itr  = 35;                                   // max number of step in NLBGS for flow subproblem
    real tNLBGS_rel_res  = 1.0e-9;                               // required convergence on pseudo-dynamic residual
    real tNLBGS_realx    = 1.0;                                  // relaxation parameter
    sint tNLBGS_init_itr = tPropRampSteps;                       // initialization phase during which property has not reached nominal value
    sint tNLBGS_ref_its  = std::max( tPropRampSteps + 1, 2 );    // iteration id when reference residual to be calculated
                                                                 // determines minimum number of steps in NLBFGS; at least two iteration need to be performed

    // NLBGS - Load control
    bool tUseNLBGSLoadControl = false;

    // Newton - Load control
    bool tUseNewtonLoadControl = false;

    // Load control settings used for Newton and NLBGS
    real tRampInitial = 0.01;    // initial load factor
    real tRampRelRes  = 0.01;    // max rel. residual below which load factor is increased
    sint tRampSteps   = 3.;      // number of ramping steps

    // NLBGS - Pseudo time stepping
    bool tUseLoadTimeRamp   = false;
    bool tUseGlobalTimeStep = false;

    real tPseudoTimeIndexCflOld  = -1;
    real tPseudoTimeIndexPropOld = -1;

    sint tMaxNumTimeSteps     = 35;        // max. number of steps in PTC
    real tMaxTimeStepSize     = 1.0e6;     // max. CFL
    real tRampTotalTime       = 1.0;       // total CFL corresponding to end of load ramping (only used if tUseLoadTimeRamp = true)
    real tRelResNormDrop      = 1.0e-9;    // required rel. static res drop for convergence
    real tRelResNormUpdate    = 1.0e1;     // required rel. static res drop for stepping forward in time
    real tRelResInexNewton    = -1.0;      // required rel. static res drop for switching to inexact Newton
    real tSteadyStateRelRes   = -1.0;      // required rel. static res drop for switching to exact Newton
    real tSteadyStateStepSize = 1.0e6;     // CFL in exact Newton
    real tTimeOffSet          = 0.0;       // time offset for writing PTC steps (= 0: no output)

    real tConstantTimeStep = 1.0;
    real tIndexFactor      = 2.0;
    real tIndexExpo        = 1.0;

    real tResFactor = 1.8;     // for Expur: increase factor
    real tResExpo   = 0.33;    // for Expur: decrease factor

    real tComsolParameter1 = 20.0;
    real tComsolParameter2 = 30.0;

    // Newton paramters when using NLBGS
    real tNewton_rel_res  = 2.5e-1;
    real tNewton_relax    = 1.0;
    sint tNewton_max_iter = 15;

    // Newton paramters without NLBGS
    moris::real tNLA_rel_res_norm_drop    = 1.0e-9;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 100;
    sint        tNLA_ref_its              = tSolverType == "Newton" ? std::max( tPropRampSteps + 1, 2 ) : 1;

    // Time solver parameters
    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0;

    std::string tVXInitial = tUseNLBGSLoadControl ? std::to_string( tRampInitial ) : tUseLoadTimeRamp ? std::to_string( 1.0 / tRampTotalTime )
                                                                                                      : "0.001";

    /* ------------------------------------------------------------------------ */

    // derived reference values
    real pressref = rhoref * velref * velref;
    real timeref  = lenref / velref;
    real massref  = rhoref * std::pow( lenref, 3.0 );

    // scaling parameters: input: m -> used for computation: mm
    real tLengthScale = 1.0 / lenref;     // 1 m  = 1.0/lref;
    real tTimeScale   = 1.0 / timeref;    // 1 s  = 1.0 s
    real tMassScale   = 1.0 / massref;    // 1 kg = 1.0 kg

    real tForceScale    = tMassScale * tLengthScale / tTimeScale / tTimeScale;
    real tPressureScale = tForceScale / tLengthScale / tLengthScale;
    real tEnergyScale   = tForceScale * tLengthScale;
    real tDensityScale  = tMassScale / tLengthScale / tLengthScale / tLengthScale;

    /* ------------------------------------------------------------------------ */
    // dimensionality 2D or 3D
    uint gDim = 2;

    // Geometry Parameters (m)
    moris::real sH = 0.1 * tLengthScale;

    // bounding box of computational domain
    moris::real tOffsetInOutlet = 2.0;    // original
    moris::real tOffsetSide     = 1.0;
    moris::real tOffsetCut      = 0.0214;

    moris::real tDimX = ( tOffsetInOutlet + 12.0 + tOffsetSide ) * sH + tOffsetCut;    // x-dimension of computational domain
    moris::real tDimY = ( tOffsetSide + 10.0 + tOffsetSide ) * sH + tOffsetCut;        // y-dimension of computational domain

    moris::real tOffsetX = -tOffsetInOutlet * sH - tOffsetCut;
    moris::real tOffsetY = -tOffsetSide * sH - tOffsetCut;

    moris::real tApproxEleSize = 1.0 / (tNumElements)*tLengthScale;

    /* ------------------------------------------------------------------------ */
    // background mesh
    std::string tNumElemX = moris_to_string( std::ceil( tDimX / tApproxEleSize ) + 1 );
    std::string tNumElemY = moris_to_string( std::ceil( tDimY / tApproxEleSize ) + 1 );

    std::string tDomainDimX = moris_to_string( tDimX );
    std::string tDomainDimY = moris_to_string( tDimY );

    std::string tDomainOffX = moris_to_string( tOffsetX );
    std::string tDomainOffY = moris_to_string( tOffsetY );

    std::string tNumElemsPerDim = tNumElemX + "," + tNumElemY;
    std::string tDomainDims     = tDomainDimX + "," + tDomainDimY;
    std::string tDomainOffset   = tDomainOffX + "," + tDomainOffY;
    std::string tDomainSidesets = "1,2,3,4";

    // Bspline limit
    int tLevelsetOrder = 2;
    int tDispOrder     = 1;
    int tLagMeshOrder  = 1;

    int tLevelsetInitialRef = 0;
    int tDispInitialRef     = 3;
    int tRefineBuffer       = 1;

    // note:  pattern 0 - displacement field   pattern 1 - Levelset field
    std::string tLagrangeOrder   = std::to_string( tLagMeshOrder );
    std::string tBsplineOrder    = std::to_string( tDispOrder ) + "," + std::to_string( tLevelsetOrder );
    std::string tInitialRef      = std::to_string( tDispInitialRef ) + "," + std::to_string( tLevelsetInitialRef );
    std::string tLagrangePattern = tDispInitialRef > tLevelsetInitialRef ? "0" : "1";

    uint tInterfaceRefinementInclusion = 0;

    // Size of FEM element
    moris::real tElementEdgeLength = tApproxEleSize / ( std::pow( 2, tDispInitialRef ) );

    // Bspline limit
    moris::real tBSplineLimit = 2.0 * tApproxEleSize / ( std::pow( 2, tLevelsetInitialRef ) );
    moris::real tPhiBandwidth = 3.0 * tApproxEleSize / ( std::pow( 2, tLevelsetInitialRef ) );
    moris::real tPhiGradient  = tBSplineLimit * std::log( 199.0 ) / ( 2.0 * tPhiBandwidth );
    moris::real tPhiGamma     = 2.0 * std::log( 10.0 ) / std::pow( 2.0 / ( std::exp( -2.0 * tPhiBandwidth * tPhiGradient / tBSplineLimit ) + 1.0 ) - 1.0, 2.0 );

    // heat method
    std::string tWPhi1     = "0.0";
    std::string tWPhi2     = "0.5";
    std::string tWGradPhi1 = "0.0";
    std::string tWGradPhi2 = "0.5";

    /* ------------------------------------------------------------------------ */

    moris::real tMMAPenalty  = 50.0;
    moris::real tMMAStepSize = 0.01;
    int         tMMAMaxIter  = 0;

    moris::real sWeightPowDis = 1.0;    // weight on temperature obj
    moris::real sWeightPerim  = 0.0;    // 0.1;    // weight on perimeter
    moris::real sWeightHM     = 0.0;    // 0.5;    // weight on heat method

    moris::real sFldVolPenalty = 1.0;    // constraint contribution to objective

    // formulation of constraint
    moris::real sWeightFluidV   = 1.0;    // weight on max fluid volume
    moris::real sFractionFluidV = 0.3;    // weight on max fluid volume

    // scaling objective and constraints
    moris::real gRefPowDis = 16.6613;
    moris::real gRefPerim  = 84.4471;
    moris::real gRefHM     = 32.8321;
    moris::real gRefFluidV = sFractionFluidV * ( 10.0 * sH * 10.0 * sH - 4.5 * sH * sH - M_PI * sH * sH / 8.0 ) + 2.0 * ( tOffsetInOutlet + 2.0 ) * sH * 2.0 * sH;

    /* ------------------------------------------------------------------------ */

    // Fluid phase mesh set - 1
    std::string tFluid = "HMR_dummy_c_p1,HMR_dummy_n_p1";

    // Inclusion phase mesh set - 3
    std::string tInclusionSolid = "HMR_dummy_c_p2,HMR_dummy_n_p2";

    // All domain
    std::string tAllDomain = tFluid + "," + tInclusionSolid;

    /* ------------------------------------------------------------------------ */
    // material parameters for wall distance

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
    // FD
    std::string tFDEpsilon = "1e-5";
    std::string tFDSweep   = "1e-5";

    std::string tAdvIndices = "30";

    /* ------------------------------------------------------------------------ */
    // File names
    std::string tName          = "uBend_2D_Optimization";
    std::string tExoFile       = tName + ".exo";
    std::string tSoFile        = tName + ".so";
    std::string tHdf5File      = tName + ".hdf5";
    std::string tGENOutputFile = tName + "_GEN.exo";

    /* ------------------------------------------------------------------------ */
    // material parameters

    // dynamic viscosity
    real tDynVis = tDensFld * tKinVisc;

    // reynolds number based on reference values and fluid properties see below
    real tDimReynolds = lenref * velref / tKinVisc;

    std::string tFluidDensity      = moris_to_string( tDensFld * tDensityScale );                               // kg/m3   air at 300 K and 1 atm
    std::string tFluidDynViscosity = moris_to_string( tDynVis * tPressureScale * tTimeScale );                  // N s/m2  air at 300 K
    std::string tFluidKinViscosity = moris_to_string( tKinVisc * tLengthScale * tLengthScale / tTimeScale );    // m2/s air at 300 K

    std::string tFluidPressureSpring = "1e-12";

    /* ------------------------------------------------------------------------ */
    // boundary condition parameters

    real tInletTurbVisc = tTurbIntensity < 0.0 ? -tTurbIntensity * tKinVisc : std::sqrt( 1.5 ) * velref * lenref * tTurbIntensity;

    std::string tInletKinViscosity = moris_to_string( tInletTurbVisc * tLengthScale * tLengthScale / tTimeScale );    // m2/s air at 800 C
    std::string tInletVelocity     = moris_to_string( tInVel * tLengthScale / tTimeScale );

    /* ------------------------------------------------------------------------ */

    // use ghost
    bool tUseGhost = true;

    std::string alpha_velo = "5.0e-03";    // scaling multiplier for velocity term in nitsche
    std::string alpha_pres = "5.0e-04";    // scaling multiplier for pressure term in nitsche
    std::string alpha_conv = "5.0e-03";    // scaling multiplier for convective term in nitsche
    std::string alpha_temp = "5.0e-03";    // scaling multiplier for incompressible term in nitsche
    std::string alpha_visc = "5.0e-03";
    std::string alpha_heat = "5.0e-03";
    std::string alpha_l2   = "5.0e-03";

    // Nitsche and Ghost parameters for thermal problem
    std::string sIfcNitscheFluid     = "1.0e+1";
    std::string sIfcNitscheTherm     = "1.0e+1";
    std::string sIfcNitscheVisc      = "1.0e+1";
    std::string sIfcNitscheViscosity = "1.0e+1";
    std::string sIfcNitscheHeat      = "1.0e+1";

    /* ------------------------------------------------------------------------ */
    // select VY velocity function
    void
    Func_Select_Inlet_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );

        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if ( tY >= 5.5 * sH && tY <= 7.5 * sH )
        {
            aPropMatrix( 0, 0 ) = 1.0;
            aPropMatrix( 1, 1 ) = 1.0;
        }
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_Inlet_Upwind(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );

        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if ( tY >= 5.5 * sH && tY <= 7.5 * sH )
        {
            aPropMatrix( 0, 0 ) = aParameters( 0 )( 0 );
        }
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_Select_Inlet_V(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );

        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if ( tY >= 5.5 * sH && tY <= 7.5 * sH )
        {
            aPropMatrix( 0, 0 ) = 1.0;
        }
    }

    /* ------------------------------------------------------------------------ */
    // pseudo time stepping function
    void
    Func_Time_Weights(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        if ( gLogger.exists( "MDL", "Model", "Perform Forward Analysis" ) )
        {
            real tWeight = aParameters( 0 )( 0 );

            real tCFL = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "PseudoTimeStep" );

            const Matrix< DDRMat >& tXhat = aFIManager->get_IP_geometry_interpolator()->get_space_coeff();

            real tXmax = tXhat.get_column( 0 ).max();
            real tXmin = tXhat.get_column( 0 ).min();

            real tEleLength = tXmax - tXmin;

            real tVelocNorm = 1.0;

            if ( !tUseGlobalTimeStep )
            {
                moris::fem::Field_Interpolator_Manager* tPrevFIManager = aFIManager->get_field_interpolator_manager_previous();

                const Matrix< DDRMat >& tVelocity = tPrevFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val();

                tVelocNorm = norm( tVelocity );
            }

            aPropMatrix( 0 ) = tWeight * tVelocNorm / std::max( MORIS_REAL_EPS, tCFL * tEleLength );

            real tPseudoTimeStep = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "PseudoTimeStep" );

            if ( tPseudoTimeStep < 1.0 )
            {
                aPropMatrix( 0 ) = 0.0;
            }

            // get pseudo time indes of most recent NLBGS iteration (last parameter set to true)
            uint tPseudoTimeIndex = gLogger.get_iteration( "NonLinearAlgorithm", "NLBGS", "Solve", true );

            if ( std::abs( tPseudoTimeIndex - tPseudoTimeIndexCflOld ) > 1e-8 )
            {
                MORIS_LOG_INFO( "PseudoTimeIndex = %d  in Func_Time_Weights: tPseudoTimeStep = %f  tCFL = %e", tPseudoTimeIndex, tPseudoTimeStep, tCFL );

                tPseudoTimeIndexCflOld = tPseudoTimeIndex;
            }
        }
        else if ( gLogger.exists( "MDL", "Model", "Perform Sensitivity Analysis" ) )
        {
            aPropMatrix( 0 ) = 0.0;
        }
        else
        {
            MORIS_ERROR( false, "Func_Time_Weights - Neither FA, nor SA." );
            aPropMatrix( 0 ) = 0.0;
        }
    }

    /* ------------------------------------------------------------------------ */
    // To turn on/off time continuity residual for inexact Newton iterations
    // in this case time continuitiy residual is omitted but jacobian is still
    // computed
    void
    Func_Time_Weight_Res(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );

        if ( gLogger.exists( "MDL", "Model", "Perform Forward Analysis" ) )
        {
            real tStaticRelRes = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "RelativeStaticResidual" );
            aPropMatrix( 0 )   = tStaticRelRes < tRelResInexNewton ? 0.0 : 1.0;
        }
        else if ( gLogger.exists( "MDL", "Model", "Perform Sensitivity Analysis" ) )
        {
            aPropMatrix( 0 ) = 0.0;
        }
        else
        {
            MORIS_ERROR( false, "Func_Time_Weights - Neither FA, nor SA." );
            aPropMatrix( 0 ) = 0.0;
        }
    }

    /* ------------------------------------------------------------------------ */
    // To decrease viscosity in initialization phase of NLBGS
    void
    Func_Ramp_Property(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tProperty = aParameters( 0 )( 0 );

        if ( gLogger.exists( "MDL", "Model", "Perform Forward Analysis" ) )
        {
            // get pseudo time indes of most recent NLBGS iteration (last parameter set to true)
            uint tPseudoTimeIndex = gLogger.get_iteration( "NonLinearAlgorithm", tSolverType, "Solve", true );

            // return original property value if NLBGS-Solve does not exist (e.g. in visualization)
            if ( tPseudoTimeIndex == 0 )
            {
                aPropMatrix = tProperty;
            }
            else
            {
                const real tMaxProperty      = tPropRampScaling * tProperty;
                const bool tLogInterpolation = false;

                // compute scaling factor; note conversion from uint to real to properly handle negative values
                real tFactor = std::max( 0.0, ( (real)tPropRampSteps - (real)tPseudoTimeIndex + 1.0 ) / ( (real)tPropRampSteps ) );

                if ( tLogInterpolation )
                {
                    aPropMatrix = std::pow( 10.0, tFactor * std::log10( tMaxProperty ) + ( 1.0 - tFactor ) * std::log10( tProperty ) );
                }
                else
                {
                    aPropMatrix = tFactor * tMaxProperty + ( 1.0 - tFactor ) * tProperty;
                }
                if ( std::abs( tPseudoTimeIndex - tPseudoTimeIndexPropOld ) > 1e-8 )
                {
                    MORIS_LOG_INFO( "PseudoTimeIndex = %d  factor = %e  nominal Property = %e  actual Property = %e", tPseudoTimeIndex, tFactor, tProperty, aPropMatrix( 0 ) );

                    tPseudoTimeIndexPropOld = tPseudoTimeIndex;
                }
            }
        }
        else
        {
            aPropMatrix = tProperty;
        }
    }

    /* ------------------------------------------------------------------------ */
    // select integration domain on fluid
    void
    Func_Select_Int_Domain(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // init property container
        aPropMatrix.set_size( 1, 1, 0.0 );

        // grab x-coordinate
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );

        // if within selected domain in x
        if ( tX >= 2.0 && tX <= 12.0 * sH )
        {
            // grab y-coordinate
            real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            // if within selected domain in y
            if ( tY >= 0.0 && tY <= 10.0 * sH )
            {
                // set property to 1
                aPropMatrix( 0 ) = 1.0;
            }
        }
    }

    /* ------------------------------------------------------------------------ */
    // inlet velocity function
    void
    Func_Inlet_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 2, 1, 0.0 );

        moris::real tInVelocity = aParameters( 0 )( 0 );

        real tScalingFactor;
        if ( tUseLoadTimeRamp )
        {
            tScalingFactor = std::min( 1.0, 0.01 + 0.99 * gLogger.get_action_data( "NonLinearAlgorithm", tSolverType, "Solve", "PseudoTotalTime" ) / 50.0 );
        }
        else
        {
            tScalingFactor = std::min( 1.0, gLogger.get_action_data( "NonLinearAlgorithm", tSolverType, "Solve", "LoadFactor" ) );

            real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            if ( tPowerlawProfile )
            {
                real tPower      = 2.0 * std::log10( tDimReynolds / 10.0 );
                aPropMatrix( 0 ) = tScalingFactor * tInVelocity * ( std::pow( 1.0 - ( ( tY - 6.5 * sH ) / std::abs( tY - 6.5 * sH ) ) * ( ( tY - 6.5 * sH ) / sH ), 1.0 / tPower ) );
            }
            else
            {
                aPropMatrix( 0 ) = tScalingFactor * tInVelocity * ( 1.0 - std::pow( ( tY - 6.5 * sH ) / sH, 2.0 ) );
            }
        }
    }

    /* ------------------------------------------------------------------------ */
    // inlet viscosity function
    void
    Func_Inlet_V(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );

        moris::real tInVisc = aParameters( 0 )( 0 );

        real tScalingFactor = tRampViscBC ? std::min( 1.0, gLogger.get_action_data( "NonLinearAlgorithm", tSolverType, "Solve", "LoadFactor" ) ) : 1.0;

        aPropMatrix( 0 ) = tScalingFactor * tInVisc;
    }

    /* ------------------------------------------------------------------------ */
    // wall distance function
    void
    Func_Wall_Distance(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tL2Val = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->val()( 0 );
        aPropMatrix.set_size( 1, 1, std::max( std::sqrt( std::pow( 10.0, -tL2Val ) ), 1e-9 ) );
    }

    /* ------------------------------------------------------------------------ */
    // wall distance derivative function
    void
    Func_Wall_Distance_Der(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tL2Val = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->val()( 0 );

        real tFuncVal = std::sqrt( std::pow( 10.0, -tL2Val ) );

        if ( tFuncVal > 1e-9 )
        {
            aPropMatrix = -log( 10 ) * std::pow( 10.0, -tL2Val ) / 2.0 / tFuncVal * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->N();
        }
        else
        {
            aPropMatrix = 0.0 * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->N();
        }
    }

    /* ------------------------------------------------------------------------ */
    // inverse wall distance squared function
    void
    Func_Wall_InvDistanceSquare(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tPhiDVal = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );

        aPropMatrix.set_size( 1, 1, std::log10( 1.0 / std::pow( std::max( 1e-9, tPhiDVal ), 2.0 ) ) );
    }

    // inverse wall distance squared function derivative
    void
    Func_Wall_InvDistanceSquare_Der(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tPhiDVal = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );
        if ( tPhiDVal > 1e-9 )
        {
            aPropMatrix = -2.0 / tPhiDVal / log( 10.0 ) * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->N();
        }
        else
        {
            aPropMatrix = 0.0 * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->N();
        }
    }

    // prescribed gradient when using non-pdv level set field
    std::string tLevelSetGradxConstant = "1.0;1.0";

    bool tUseAbsoulteValue = false;

    /* ------------------------------------------------------------------------ */

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N();
    }

    /* ------------------------------------------------------------------------ */

    // Spatial derivative of level set function defining property in FEM
    void
    tLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->gradx( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of spatial derivative of level set function with respect to PDV
    void
    tDerLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoulteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->dnNdxn( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Minimum level set value
    moris::real tInclusionRadius = sH / 1.0;          // inclusion radius
    moris::real tLeftBound       = 2.0 * sH + sH;     // 2.0 * tInclusionRadius;     // min x coord for bounding box
    moris::real tRightBound      = 12.0 * sH - sH;    // 2.0 * tInclusionRadius;     // max x coord for bounding box
    moris::real tBottomBound     = 0.0 * sH + sH;     // 2.0 * tInclusionRadius;     // min y coord for bounding box
    moris::real tTopBound        = 10.0 * sH - sH;    // 2.0 * tInclusionRadius;     // max y coord for bounding box

    uint tNumXRods = 3;    // number of inclusions along x
    uint tNumYRods = 3;    // number of inclusions along y

    moris::real tInclusionExponent = 12.0;    // inclusion exponent

    moris::real
    Func_Inclusion(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< moris::real >&   aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        // small rods
        real tReturnValue = 1.0;

        for ( uint iYRod = 0; iYRod < tNumYRods; iYRod++ )
        {
            real tYCenter = tBottomBound + iYRod * ( tTopBound - tBottomBound ) / ( tNumYRods - 1 );

            uint tTempNumXRods = tNumXRods;
            real tXOffset      = 0.0;

            // if ( fmod( iYRod, 2 ) != 0 )
            //{
            //     tXOffset = ( tRightBound - tLeftBound ) / ( 2 * ( tNumXRods - 1 ) );
            //     tTempNumXRods--;
            // }

            for ( uint iXRod = 0; iXRod < tTempNumXRods; iXRod++ )
            // for( uint iXRod = 0; iXRod < tNumXRods; iXRod++ )
            {
                real tXCenter = tLeftBound + tXOffset + iXRod * ( tRightBound - tLeftBound ) / ( tNumXRods - 1 );
                // real tXCenter = tLeftBound + iXRod * ( tRightBound - tLeftBound ) / ( tNumXRods - 1 );

                // skip area near the thin wall
                if ( tY > 6.0 * sH || tY < 4.0 * sH || true )
                {
                    real tCylValue = -tInclusionRadius
                                   + std::pow( std::pow( tX - tXCenter, tInclusionExponent )
                                                       + std::pow( tY - tYCenter, tInclusionExponent ),
                                           1.0 / tInclusionExponent );

                    tReturnValue = std::min( tReturnValue, tCylValue );
                }
                else
                {
                    // skip area near the thin wall
                    if ( tX > 2.5 * sH )    // 8.5 * sH )
                    {
                        real tCylValue = -tInclusionRadius
                                       + std::pow( std::pow( tX - tXCenter, tInclusionExponent )
                                                           + std::pow( tY - tYCenter, tInclusionExponent ),
                                               1.0 / tInclusionExponent );

                        tReturnValue = std::min( tReturnValue, tCylValue );
                    }
                }
            }
        }
        return tReturnValue;
    }

    /* ------------------------------------------------------------------------ */

    moris::real
    Func_Thin_Wall(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< moris::real >&   aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );

        // parameters
        real tXctr = aGeometryParameters( 0 );    // 6.5 * sH;
        real tYctr = aGeometryParameters( 1 );    // 5.0 * sH;
        real tXRad = aGeometryParameters( 2 );    // sH;
        real tYRad = aGeometryParameters( 3 );    // sH;
        real tExpn = aGeometryParameters( 4 );    // 2.0;
        real tBot  = aGeometryParameters( 5 );    // 4.5 * sH;
        real tTop  = aGeometryParameters( 6 );    // 5.5 * sH;

        // init return value
        real tReturnValue = 1.0;

        // compute value of 3 LS (2 planes, 1 circle)
        real tLS0Value = -tY + tBot;
        real tLS1Value = tY - tTop;

        // combine the LS
        tReturnValue = std::max( tLS0Value, tLS1Value );

        if ( tX >= tXctr )
        {
            real tLS2Value = -1.0 + std::pow( std::pow( ( tX - tXctr ) / tXRad, tExpn ) + std::pow( ( tY - tYctr ) / tYRad, tExpn ), 1.0 / tExpn );
            tReturnValue   = std::max( tReturnValue, tLS2Value );
        }

        return tReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // Phase assignement

    std::string tGetPhaseIndex = "get_phase_index_sharp";

    uint get_phase_index_sharp( const Bitset< 9 >& aGeometrySigns )
    {
        // Phase Void
        if ( !aGeometrySigns.test( 0 ) )
        {
            return 0;
        }

        if ( aGeometrySigns.test( 3 ) )
        {
            return 0;
        }

        if ( aGeometrySigns.test( 1 ) )
        {
            return 0;
        }

        if ( !aGeometrySigns.test( 2 ) && aGeometrySigns.test( 5 ) )
        {
            return 0;
        }

        if ( !aGeometrySigns.test( 4 ) )
        {
            return 0;
        }

        if ( !aGeometrySigns.test( 2 ) && aGeometrySigns.test( 5 ) )
        {
            return 0;
        }

        // Phase Inclusion
        if ( aGeometrySigns.test( 0 ) && !aGeometrySigns.test( 1 ) && aGeometrySigns.test( 2 ) && !aGeometrySigns.test( 3 ) && !aGeometrySigns.test( 6 ) )
        {
            return 2;
        }

        // in the fluid
        return 1;
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */
    /*
    0 - IQIVolumePowDisp
    1 - IQIPerimeterItf
    2 - IQIHeatMethodPenaltyFluid
    3 - IQIHeatMethodPenaltyInclusion
    4 - IQIFluidVolume
    */

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, tNumConstraints, 1 );

        return tConstraintTypes;
    }

    Matrix< DDRMat >
    compute_objectives(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        // set init values for scaling
        if ( gLogger.get_opt_iteration() == 1 )
        {
            gRefPowDis = aCriteria( 0 );
            gRefPerim  = aCriteria( 1 );
            gRefHM     = aCriteria( 2 ) + aCriteria( 3 );

            std::cout << "gRefPowDis = " << gRefPowDis << " \n";
            std::cout << "gRefPerim  = " << gRefPerim << " \n";
            std::cout << "gRefHM     = " << gRefHM << " \n";
        }

        moris::real scaledPowDis = sWeightPowDis * aCriteria( 0 ) / gRefPowDis;
        moris::real scaledPerim  = sWeightPerim * aCriteria( 1 ) / gRefPerim;
        moris::real scaledHM     = sWeightHM * ( aCriteria( 2 ) + aCriteria( 3 ) ) / gRefHM;

        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = scaledPowDis + scaledPerim + scaledHM;

        std::cout << "Contributions to objective function:\n";

        std::cout << "Power dissipation   = " << aCriteria( 0 ) << " \n";
        std::cout << "Perimeter           = " << aCriteria( 1 ) << " \n";
        std::cout << "HM Penalty          = " << aCriteria( 2 ) + aCriteria( 3 ) << " \n";

        std::cout << "Scaled power dissipation   = " << scaledPowDis << " \n";
        std::cout << "Scaled Perimeter           = " << scaledPerim << " \n";
        std::cout << "Scaled HM penalty          = " << scaledHM << " \n";

        std::cout << "Objective                  = " << tObjectives( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        std::cout << "min ADV             = " << aADVs.min() << " \n";
        std::cout << "max ADV             = " << aADVs.max() << " \n"
                  << std::flush;
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

        return tObjectives;
    }

    Matrix< DDRMat >
    compute_dobjective_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

        tDObjectiveDCriteria( 0 ) = sWeightPowDis / gRefPowDis;
        tDObjectiveDCriteria( 1 ) = sWeightPerim / gRefPerim;
        tDObjectiveDCriteria( 2 ) = sWeightHM / gRefHM;
        tDObjectiveDCriteria( 3 ) = sWeightHM / gRefHM;

        return tDObjectiveDCriteria;
    }

    Matrix< DDRMat >
    compute_constraints(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );

        if ( tCheckSensitivities )
        {
            tConstraints( 0 ) = aCriteria( 0 );
            tConstraints( 1 ) = aCriteria( 1 );
            tConstraints( 2 ) = aCriteria( 4 );
        }
        else
        {
            // fluid volume
            tConstraints( 0 ) = ( aCriteria( 4 ) / gRefFluidV ) - 1.0;

            std::cout << "Fluid volume                = " << aCriteria( 4 ) << " \n";
            std::cout << "Scaled Fluid volume         = " << aCriteria( 4 ) / gRefFluidV << " \n";
            std::cout << "Fluid volume constraint     = " << tConstraints( 0 ) << " \n";
            std::cout << "% --------------------------------- % \n"
                      << std::flush;
        }

        return tConstraints;
    }

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.size(), 0.0 );

        if ( tCheckSensitivities )
        {
            tDConstraintDCriteria( 0, 0 ) = 1.0;
            tDConstraintDCriteria( 1, 1 ) = 1.0;
            tDConstraintDCriteria( 2, 4 ) = 1.0;
        }
        else
        {
            tDConstraintDCriteria( 0, 4 ) = 1.0 / gRefFluidV;
        }
        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tSoFile );

        if ( tCheckSensitivities )
        {
            aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
            aParameterLists.set( "hdf5_path", tHdf5File );
            aParameterLists.set( "num_evaluations_per_adv", "1" );
            aParameterLists.set( "finite_difference_type", "all" );
            aParameterLists.set( "finite_difference_epsilons", tFDSweep );
            aParameterLists.set( "finite_difference_adv_indices", tAdvIndices );
        }
        else
        {
            aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
            aParameterLists.set( "step_size", tMMAStepSize );
            aParameterLists.set( "penalty", tMMAPenalty );
            aParameterLists.set( "max_its", tMMAMaxIter );
        }

        if ( tRestartId > 0 )
        {
            aParameterLists.set( "restart_file", "ADV_Alg_0_Iter_" + std::to_string( tRestartId ) + ".hdf5" );
            aParameterLists.set( "restart_index", tRestartId );
        }
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "domain_offset", tDomainOffset );
        aParameterLists.set( "domain_sidesets", tDomainSidesets );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", tLagrangePattern );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "initial_refinement", tInitialRef );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        aParameterLists.set( "use_advanced_T_matrix_scheme", 1 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        // aParameterLists.set( "cleanup_cut_mesh",            true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types",
                "IQIVolumePowDisp",
                "IQIPerimeterItf",
                "IQIHeatMethodPenaltyFluid",
                "IQIHeatMethodPenaltyInclusion",
                "IQIFluidVolume" );
        aParameterLists.set( "number_of_phases", 3 );
        aParameterLists.set( "phase_function_name", tGetPhaseIndex );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        // Plane 0 in y = 0
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Plane 1 in y = 10.0*sH
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 0.0 );
        aParameterLists.set( "center_y", 10.0 * sH );
        aParameterLists.set( "normal_x", 0.0 );
        aParameterLists.set( "normal_y", 1.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Plane 2 in x = 0.0*sH
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 2.0 * sH );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Plane 3 in x = 12.0*sH
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::LINE );
        aParameterLists.set( "center_x", 12.0 * sH );
        aParameterLists.set( "center_y", 0.0 );
        aParameterLists.set( "normal_x", 1.0 );
        aParameterLists.set( "normal_y", 0.0 );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Thin wall inner 4
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Thin_Wall" );
        aParameterLists( 1 ).insert( "variable_1", 6.5 * sH );
        aParameterLists( 1 ).insert( "variable_2", 5.0 * sH );
        aParameterLists( 1 ).insert( "variable_3", sH / 2.0 );
        aParameterLists( 1 ).insert( "variable_4", sH / 2.0 );
        aParameterLists( 1 ).insert( "variable_5", 2.0 );
        aParameterLists( 1 ).insert( "variable_6", 4.5 * sH );
        aParameterLists( 1 ).insert( "variable_7", 5.5 * sH );
        aParameterLists.set( "name", "Thin_Wall_Inner" );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Thin wall outer 5
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Thin_Wall" );
        aParameterLists( 1 ).insert( "variable_1", 6.5 * sH );
        aParameterLists( 1 ).insert( "variable_2", 5.0 * sH );
        aParameterLists( 1 ).insert( "variable_3", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_4", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_5", 2.0 );
        aParameterLists( 1 ).insert( "variable_6", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_7", 7.5 * sH );
        aParameterLists.set( "name", "Thin_Wall_Outer" );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );

        // Inclusions 6
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Func_Inclusion" );
        aParameterLists.set( "name", "Level_Set_Field" );
        aParameterLists.set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists.set( "refinement_mesh_index", 0 );
        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 1 );
            aParameterLists.set( "discretization_lower_bound", -tBSplineLimit );
            aParameterLists.set( "discretization_upper_bound", tBSplineLimit );
        }

        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );
        uint tParamCounter = 0;

        aParameterLists.set( "name", "LevelsetField" );
        aParameterLists.set( "dependencies", "Level_Set_Field" );
        aParameterLists.set( "scaling_factor", 1.0 );
        aParameterLists.set( "pdv_type", "LS1" );
        aParameterLists.set( "pdv_mesh_set_names", tFluid + "," + tInclusionSolid );
        tParamCounter++;
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        if ( par_rank() == 0 )
        {
            std::cout << "Reynolds       " << tDimReynolds << "------------------------------------" << '\n';

            std::cout << "velref         " << velref << "------------------------------------" << '\n';
            std::cout << "lenref         " << lenref << "------------------------------------" << '\n';
            std::cout << "rhoref         " << rhoref << "------------------------------------" << '\n';
            std::cout << "pressref       " << pressref << "------------------------------------" << '\n';
            std::cout << "timeref        " << timeref << "------------------------------------" << '\n';
            std::cout << "massref        " << massref << "------------------------------------" << '\n';

            std::cout << "tLengthScale   " << tLengthScale << "------------------------------------" << '\n';
            std::cout << "tTimeScale     " << tTimeScale << "------------------------------------" << '\n';
            std::cout << "tMassScale     " << tMassScale << "------------------------------------" << '\n';

            std::cout << "tPressureScale " << tPressureScale << "------------------------------------" << '\n';
            std::cout << "tDensityScale  " << tDensityScale << "------------------------------------" << '\n';

            std::cout << "tFluidDynViscosity     " << tFluidDynViscosity << "------------------------------------" << '\n';
            std::cout << "tFluidKinViscosity     " << tFluidKinViscosity << "------------------------------------" << '\n';
            std::cout << "tFluidDensity          " << tFluidDensity << "------------------------------------" << '\n';

            std::cout << "tInletKinViscosity       " << tInletKinViscosity << "------------------------------------" << '\n';
            std::cout << "tInletVelocity           " << tInletVelocity << "------------------------------------" << '\n';
            std::cout << "tElementEdgeLength       " << tElementEdgeLength << "------------------------------------" << '\n';

            std::cout << "tBSplineLimit            " << tBSplineLimit << " (rel wrt fem element : " << tBSplineLimit / tElementEdgeLength << " )\n";
            std::cout << "tPhiBandwidth            " << tPhiBandwidth << " (rel wrt fem element : " << tPhiBandwidth / tElementEdgeLength << " )\n";
            std::cout << "tPhiGradient             " << tPhiGradient << "\n";
            std::cout << "tPhiGamma                " << tPhiGamma << "\n";
        }

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "phase_indices", "1" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseInclusion" );
        aParameterLists.set( "phase_indices", "2" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseAll" );
        aParameterLists.set( "phase_indices", "1,2" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseHeatMethod" );
        aParameterLists.set( "phase_indices", "1,2" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid properties ------------------------------------------------------------
        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidDensity" );
        aParameterLists.set( "function_parameters", tFluidDensity );

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidViscosity" );
        aParameterLists.set( "function_parameters", tFluidDynViscosity );
        aParameterLists.set( "value_function", "Func_Ramp_Property" );

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidKinViscosity" );
        aParameterLists.set( "function_parameters", tFluidKinViscosity );
        aParameterLists.set( "value_function", "Func_Ramp_Property" );

        // create fluid pressure spring property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFluidPressureSpring" );
        aParameterLists.set( "function_parameters", tFluidPressureSpring );

        // BC properties ---------------------------------------------------------------
        // create inlet velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletU" );
        aParameterLists.set( "value_function", "Func_Inlet_U" );
        aParameterLists.set( "function_parameters", tInletVelocity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelectInletU" );
        aParameterLists.set( "value_function", "Func_Select_Inlet_U" );

        // create zero velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletZeroU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // create wall distance property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWallDistance" );
        aParameterLists.set( "dof_dependencies", "L2" );
        aParameterLists.set( "value_function", "Func_Wall_Distance" );
        aParameterLists.set( "dof_derivative_functions", "Func_Wall_Distance_Der" );

        // create wall distance property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInvWallDistanceSquare" );
        aParameterLists.set( "dof_dependencies", "PHID" );
        aParameterLists.set( "value_function", "Func_Wall_InvDistanceSquare" );
        aParameterLists.set( "dof_derivative_functions", "Func_Wall_InvDistanceSquare_Der" );

        // inlet viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletV" );
        aParameterLists.set( "value_function", "Func_Inlet_V" );
        aParameterLists.set( "function_parameters", tInletKinViscosity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelectInletV" );
        aParameterLists.set( "value_function", "Func_Select_Inlet_V" );

        // zero viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletZeroV" );
        aParameterLists.set( "function_parameters", "0.0" );

        // upwind
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInletUpwind" );
        aParameterLists.set( "value_function", "Func_Inlet_Upwind" );
        aParameterLists.set( "function_parameters", "1.0" );

        // time continuity weights for velocity and viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightUV" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Time_Weights" );

        // time continuity weights for velocity and viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightResUV" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Time_Weight_Res" );

        // initial condition velocity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialConditionU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        // initial condition viscosity
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialConditionV" );
        aParameterLists.set( "function_parameters", tInletKinViscosity );

        // select integration domain
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelectDomain" );
        aParameterLists.set( "value_function", "Func_Select_Int_Domain" );

        // Wall distance properties ----------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );

        // properties for Theta

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityTheta" );
        aParameterLists.set( "function_parameters", tDensityTheta );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityTheta" );
        aParameterLists.set( "function_parameters", tCapacityTheta );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescTheta" );
        aParameterLists.set( "function_parameters", tPrescTheta );

        // properties for phi problem

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityPhi" );
        aParameterLists.set( "function_parameters", tDensityPhi );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityPhi" );
        aParameterLists.set( "function_parameters", tCapacityPhi );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescPhi" );
        aParameterLists.set( "function_parameters", tPrescPhi );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropEigenStrainPhi" );
        aParameterLists.set( "function_parameters", "1.0" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "10.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "10.0" );

        // initial condition
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "function_parameters", "0.0" );

        // target level set
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetConst" );
        aParameterLists.set( "function_parameters", "1.0" );

        // target level set gradient
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists.set( "function_parameters", tLevelSetGradxConstant );

        // actual level set
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSet" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // actual level set gradient
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradx" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // incompressible NS fluid constitutive model
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::FLUID_TURBULENCE );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P;VISCOSITY", "Velocity,Pressure,Viscosity" ) );
        aParameterLists.set( "properties",
                "PropFluidViscosity   ,Viscosity;"
                "PropFluidKinViscosity,KinViscosity;"
                "PropFluidDensity     ,Density" );

        // Spalart Allmaras turbulence constitutive model
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMTurbulence" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        aParameterLists.set( "function_parameters", tCMTurbFt2 + "/" + tCMTurbAlpha );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        aParameterLists.set( "properties",
                "PropFluidKinViscosity,KinViscosity;"
                "PropWallDistance     ,WallDistance" );

        //------------------------------------------------------------------------------
        // create parameter list for constitutive model - Theta problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluidDiffusionTheta" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMInclusionSolidDiffusionTheta" );
        aParameterLists.set( "phase_name", "PhaseInclusion" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluidDiffusionPhi" );
        aParameterLists.set( "phase_name", "PhaseFluid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMInclusionSolidDiffusionPhi" );
        aParameterLists.set( "phase_name", "PhaseInclusion" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // SP for fluid phase ----------------------------------------------------------
        // SUPG_PSPG NS incompressible
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPIncFlow" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists.set( "function_parameters", tSupgFluidC1 + "/0.0" );
        aParameterLists.set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // Dirichlet Nitsche for fluid velocity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPFluidDirichletNitscheU" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sIfcNitscheFluid );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // SUPG Spalart-Allmaras
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSUPGSA" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "function_parameters", tSupgTurbPower + "/" + tSupgTurbSource );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );

        // Dirichlet Nitsche for fluid turbulent viscosity
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheV" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sIfcNitscheViscosity );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

        // SP for ghost ----------------------------------------------------------------
        if ( tUseGhost )
        {
            // ghost fluid viscous
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPViscous" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
            aParameterLists.set( "function_parameters", alpha_velo );
            aParameterLists.set( "leader_properties", "PropFluidViscosity,Viscosity" );

            // ghost fluid convective
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPVelocity" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
            aParameterLists.set( "function_parameters", alpha_conv );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties", "PropFluidDensity,Density" );

            // ghost fluid pressure
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPPressure" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
            aParameterLists.set( "function_parameters", alpha_pres );
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists.set( "leader_properties",
                    "PropFluidViscosity,Viscosity;"
                    "PropFluidDensity,Density" );

            // ghost fluid turbulence viscosity
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPViscosity" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", alpha_visc );
            aParameterLists.set( "leader_properties", "PropFluidKinViscosity,Material" );

            // create parameter list for projection of the distance
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPFluidL2" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", alpha_l2 );
            aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        }

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPFluidThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "follower_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", alpha_heat );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheFluidThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sIfcNitscheHeat );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPInclusionThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "follower_phase_name", "PhaseInclusion" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", alpha_heat );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheInclusionThetaPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", sIfcNitscheHeat );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // fluid bulk IWG --------------------------------------------------------------
        // NS incompressible for velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidVelocityBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // NS incompressible for pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidPressureBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // turbulence
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTurbulenceBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPSUPGSA,SUPG" );

        // L2 projection for distance
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDistance" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::L2 );
        aParameterLists.set( "dof_residual", "L2" );
        aParameterLists.set( "leader_properties", "PropInvWallDistanceSquare,Source" );

        // fluid BC IWG ----------------------------------------------------------------
        // incompressible NS velocity Dirichlet IWG for inlet
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select;PropInletUpwind,Upwind" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );

        // incompressible NS pressure Dirichlet IWG for inlet
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero velocity for velocity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidZeroVelocity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );

        // zero velocity for pressure
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidZeroPressure" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        if ( tUsePseudoTimeStepping )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionU,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual" );
        }

        // inlet viscosity
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletViscosity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_properties", "PropInletV,Dirichlet;PropSelectInletV,Select;PropInletUpwind,Upwind" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // zero viscosity for walls
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGZeroViscosity" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "VISCOSITY" );
        aParameterLists.set( "leader_properties", "PropDirichletZeroV,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        if ( tUsePseudoTimeStepping )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists.set( "dof_residual", "VISCOSITY" );
            aParameterLists.set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionV,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual" );
        }

        if ( tUseGhost )
        {
            // ghost fluid viscous
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPViscous" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscous,GhostSP" );

            // ghost fluid convective
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPConvective" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VX,VY" );
            aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );

            // ghost fluid pressure
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPPressure" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "P" );
            aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );

            // ghost fluid viscosity
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPFluidViscosity" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "VISCOSITY" );
            aParameterLists.set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost L2 projection distance
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerL2" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "L2" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidL2,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // theta problem

        // theta bulk in fluid
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidDiffusionThetaBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );

        // theta bulk in inclusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInclusionSolidDiffusionThetaBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "neighbor_phases", "PhaseFluid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerTheta" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "THETA" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPOuterTheta" );
            aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
            aParameterLists.set( "follower_phase_name", "PhaseInclusion" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "THETA" );
            aParameterLists.set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGFluidTimeContinuityTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInclusionTimeContinuityTheta" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        //------------------------------------------------------------------------------
        // phid problem

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionInnerBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "neighbor_phases", "PhaseFluid" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPInnerPhi" );
            aParameterLists.set( "leader_phase_name", "PhaseFluid" );
            aParameterLists.set( "follower_phase_name", "PhaseFluid" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "PHID" );
            aParameterLists.set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGPOuterPhi" );
            aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
            aParameterLists.set( "follower_phase_name", "PhaseInclusion" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "PHID" );
            aParameterLists.set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists.set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "THETA" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkPHID" );
        aParameterLists.set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // fluid volume
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIFluidVolume" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        // inclusion volume
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQISolidVolume" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        // inclusion perimeter
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeterItf" );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "neighbor_phases", "PhaseInclusion" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        // fluid power dissipation in volume
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumePowDisp" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION_BULK );
        aParameterLists.set( "leader_properties", "PropSelectDomain,Select" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid" );

        // level set
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSet" );
        aParameterLists.set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSet,Property" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatMethodPenaltyFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatMethodPenaltyInclusion" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIProjectedLevelSetFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 3 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIProjectedLevelSetInclusion" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 3 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatAlphaFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 5 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatAlphaInclusion" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 5 );
        aParameterLists.set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists.set( "leader_phase_name", "PhaseInclusion" );

        // viscosity dof
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVISCOSITY" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VISCOSITY" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // L2 projection distance dof
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "L2" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // wall distance
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIWallDistance" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropWallDistance,Property" );

        // turbulent dynamic viscosity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTurbDynVisc" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective dynamic viscosity
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkEffDynVisc" );
        aParameterLists.set( "leader_phase_name", "PhaseFluid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );
        /*
                // SUPG fluid
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkSUPGU" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists.set( "vectorial_field_index",      0 );
                aParameterLists.set( "stabilization_parameters",   "SPIncFlow,Stabilization" );

                // PSPG fluid
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkPSPGU" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists.set( "vectorial_field_index",      1 );
                aParameterLists.set( "stabilization_parameters",   "SPIncFlow,Stabilization" );

                // SUPG viscosity
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkSUPGSA" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists.set( "vectorial_field_index",      0 );
                aParameterLists.set( "stabilization_parameters",   "SPSUPGSA,Stabilization" );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UX" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists.set( "vectorial_field_index",      1 );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UY" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists.set( "vectorial_field_index",      2 );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_P" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists.set( "vectorial_field_index",      0 );

                // SA strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name",                   "IQIBulkSpalartAllmarasSF" );
                aParameterLists.set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists.set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_SA );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence") ;
                        // production coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name", "IQIBulkProductionCoeff" );
                aParameterLists.set( "leader_phase_name", "PhaseFluid" );
                aParameterLists.set( "vectorial_field_index", 0 );
                aParameterLists.set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // wall destruction coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name", "IQIBulkWallDestructionCoeff" );
                aParameterLists.set( "leader_phase_name", "PhaseFluid" );
                aParameterLists.set( "vectorial_field_index", 1 );
                aParameterLists.set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // diffusion coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name", "IQIBulkDiffusionCoeff" );
                aParameterLists.set( "leader_phase_name", "PhaseFluid" );
                aParameterLists.set( "vectorial_field_index", 2 );
                aParameterLists.set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                 // production term
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name", "IQIBulkProductionTerm" );
                aParameterLists.set( "leader_phase_name", "PhaseFluid" );
                aParameterLists.set( "vectorial_field_index", 3 );
                aParameterLists.set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // wall destruction term
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists.set( "IQI_name", "IQIBulkWallDestructionTerm" );
                aParameterLists.set( "leader_phase_name", "PhaseFluid" );
                aParameterLists.set( "vectorial_field_index", 4 );
                aParameterLists.set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists.set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
                */
        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );
        aParameterLists.set( "is_analytical_sensitivity", false );
        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", stod( tFDEpsilon ) );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        // Newton solve in forward problem used by Newton flow subproblem (includes load control)
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_ref_iter", tNLA_ref_its );

        if ( tSolverType == "Newton" && tUseNewtonLoadControl )
        {
            aParameterLists.set( "NLA_load_control_strategy", sol::SolverLoadControlType::Exponential );
            aParameterLists.set( "NLA_load_control_factor", tRampInitial );
            aParameterLists.set( "NLA_load_control_steps", tRampSteps );
            aParameterLists.set( "NLA_load_control_relres", tRampRelRes );
            aParameterLists.set( "NLA_load_control_exponent", 1.0 );
        }

        // for linear solve - forward: theta, phi, l2 - adjoint: all problems
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        // NLBGS for overall forward solve and adjoint in theta,phi,l2,(vx,vy,p,viscosity)
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        // algorithm not used
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        // NLBGS for flow subproblem using pseudo time stepping and (optional) load control)
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLBGS_rel_res );
        aParameterLists.set( "NLA_relaxation_parameter", tNLBGS_realx );
        aParameterLists.set( "NLA_max_iter", tNLBGS_max_itr );
        aParameterLists.set( "NLA_ref_iter", tNLBGS_ref_its );

        if ( tSolverType == "NLBGS" )
        {
            if ( tUsePseudoTimeStepping )
            {
                aParameterLists.set( "NLA_pseudo_time_control_strategy", sol::SolverPseudoTimeControlType::Expur );
                aParameterLists.set( "NLA_pseudo_time_initial_steps", tNLBGS_init_itr );
                aParameterLists.set( "NLA_pseudo_time_initial", tSteadyStateStepSize );
                aParameterLists.set( "NLA_pseudo_time_max_num_steps", tMaxNumTimeSteps );
                aParameterLists.set( "NLA_pseudo_time_max_step_size", tMaxTimeStepSize );
                aParameterLists.set( "NLA_pseudo_time_rel_res_norm_drop", tRelResNormDrop );
                aParameterLists.set( "NLA_pseudo_time_rel_res_norm_update", tRelResNormUpdate );
                aParameterLists.set( "NLA_pseudo_time_steady_rel_res_norm", tSteadyStateRelRes );
                aParameterLists.set( "NLA_pseudo_time_steady_step_size", tSteadyStateStepSize );
                aParameterLists.set( "NLA_pseudo_time_offset", tTimeOffSet );
                aParameterLists.set( "NLA_pseudo_time_constant", tConstantTimeStep );
                aParameterLists.set( "NLA_pseudo_time_step_index_factor", tIndexFactor );
                aParameterLists.set( "NLA_pseudo_time_step_index_exponent", tIndexExpo );
                aParameterLists.set( "NLA_pseudo_time_residual_factor", tResFactor );
                aParameterLists.set( "NLA_pseudo_time_residual_exponent", tResExpo );
                aParameterLists.set( "NLA_pseudo_time_comsol_1", tComsolParameter1 );
                aParameterLists.set( "NLA_pseudo_time_comsol_2", tComsolParameter2 );
            }
            if ( tUseNLBGSLoadControl )
            {
                aParameterLists.set( "NLA_load_control_strategy", sol::SolverLoadControlType::Exponential );
                aParameterLists.set( "NLA_load_control_factor", tRampInitial );
                aParameterLists.set( "NLA_load_control_steps", tRampSteps );
                aParameterLists.set( "NLA_load_control_relres", tRampRelRes );
                aParameterLists.set( "NLA_load_control_exponent", 1.0 );
            }
        }

        // Newton solve used within NLBGS for flow subproblem
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNewton_rel_res );
        aParameterLists.set( "NLA_max_iter", tNewton_max_iter );
        aParameterLists.set( "NLA_relaxation_strategy", sol::SolverRelaxationType::InvResNormAdaptive );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_relaxation_damping", 0.5 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "THETA" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "PHID" );

        if ( tSolverType == "NLBGS" )
        {
            if ( tFluidViscSolver == "Mono" )    // monolithic flow problem
            {
                aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
                aParameterLists.set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
                aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
            }
            else    // flow problem in vx,vy,p
            {
                aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
                aParameterLists.set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
                aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                aParameterLists.set( "NLA_DofTypes", "VX,VY,P" );
            }
        }
        else    // using Newton solver only monolithic
        {
            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );
            aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }

        // flow viscosity problem only
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "VISCOSITY" );

        // following settings only used when solving flow problem by NLBGS
        if ( tFluidViscSolver == "Mono" )
        {
            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "4" );
            aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "2" );
            aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }
        else
        {
            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "4" );
            aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "2,3" );
            aParameterLists.set( "NLA_DofTypes", "VX,VY,P;VISCOSITY" );
        }

        // L2 projection problem
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "L2" );

        // total forward problem
        if ( tSolverType == "NLBGS" )
        {
            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );
            aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,5,4" );
            aParameterLists.set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }
        else
        {
            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );
            aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,5,2" );
            aParameterLists.set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }

        // adjoint of flow problem
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        // adjoint of total problem
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,5,7" );    // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );

        // ----------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 6 );                // using NLBGS for forward problem
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 8 );    // using NLBGS for sensitivity problem

        //------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0;VX," + tVXInitial + ";VY,0.0;P,0.0;VISCOSITY," + tInletKinViscosity + ";L2,1.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "THETA", 0 );
        aParameterLists.set( "PHID", 0 );
        aParameterLists.set( "VX", 0 );
        aParameterLists.set( "VY", 0 );
        aParameterLists.set( "P", 0 );
        aParameterLists.set( "VISCOSITY", 0 );
        aParameterLists.set( "L2", 0 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tAllDomain );
        aParameterLists.set( "Field_Names",
                "THETA,PHID,L2,"
                "VX,VY,P,VISCOSITY,"
                "TURBDYNVISC,EFFVISC,"
                "LS,WALLDIST,PRJLEVSFLD,PRJLEVSINC,ALPHAFLD,ALPHAINC" );
        aParameterLists.set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "NODAL,NODAL,NODAL,NODAL,"
                "NODAL,NODAL,"
                "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTHETA,IQIBulkPHID,IQIBulkL2,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkVISCOSITY,"
                "IQIBulkTurbDynVisc,IQIBulkEffDynVisc,"
                "IQILevelSet,IQIWallDistance,IQIProjectedLevelSetFluid,IQIProjectedLevelSetInclusion,IQIHeatAlphaFluid,IQIHeatAlphaInclusion" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
