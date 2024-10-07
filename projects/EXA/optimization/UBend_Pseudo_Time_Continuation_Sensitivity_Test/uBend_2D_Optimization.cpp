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
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", tIsOpt );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tSoFile );

        if ( tCheckSensitivities )
        {
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
            aParameterLists( 2 ).set( "hdf5_path", tHdf5File );
            aParameterLists( 2 ).set( "num_evaluations_per_adv", "1" );
            aParameterLists( 2 ).set( "finite_difference_type", "all" );
            aParameterLists( 2 ).set( "finite_difference_epsilons", tFDSweep );
            aParameterLists( 2 ).set( "finite_difference_adv_indices", tAdvIndices );
        }
        else
        {
            aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
            aParameterLists( 2 ).set( "step_size", tMMAStepSize );
            aParameterLists( 2 ).set( "penalty", tMMAPenalty );
            aParameterLists( 2 ).set( "max_its", tMMAMaxIter );
        }

        if ( tRestartId > 0 )
        {
            aParameterLists( 0 ).set( "restart_file", "ADV_Alg_0_Iter_" + std::to_string( tRestartId ) + ".hdf5" );
            aParameterLists( 2 ).set( "restart_index", tRestartId );
        }
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", tLagrangeOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", tLagrangePattern );

        aParameterLists( 0 ).set( "bspline_orders", tBsplineOrder );
        aParameterLists( 0 ).set( "bspline_pattern", "0,1" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0,1" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", tRefineBuffer );
        aParameterLists( 0 ).set( "staircase_buffer", tRefineBuffer );

        aParameterLists( 0 ).set( "initial_refinement", tInitialRef );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0,1" );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "use_advanced_T_matrix_scheme", 1 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", false );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
        // aParameterLists( 0 ).set( "cleanup_cut_mesh",            true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types",
                "IQIVolumePowDisp",
                "IQIPerimeterItf",
                "IQIHeatMethodPenaltyFluid",
                "IQIHeatMethodPenaltyInclusion",
                "IQIFluidVolume" );
        aParameterLists( 0 ).set( "number_of_phases", 3 );
        aParameterLists( 0 ).set( "phase_function_name", tGetPhaseIndex );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );
        aParameterLists( 0 ).set( "time_offset", 10.0 );

        // Plane 0 in y = 0
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Plane 1 in y = 10.0*sH
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 10.0 * sH );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Plane 2 in x = 0.0*sH
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 2.0 * sH );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Plane 3 in x = 12.0*sH
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists( 1 ).set( "center_x", 12.0 * sH );
        aParameterLists( 1 ).set( "center_y", 0.0 );
        aParameterLists( 1 ).set( "normal_x", 1.0 );
        aParameterLists( 1 ).set( "normal_y", 0.0 );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Thin wall inner 4
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Thin_Wall" );
        aParameterLists( 1 ).insert( "variable_1", 6.5 * sH );
        aParameterLists( 1 ).insert( "variable_2", 5.0 * sH );
        aParameterLists( 1 ).insert( "variable_3", sH / 2.0 );
        aParameterLists( 1 ).insert( "variable_4", sH / 2.0 );
        aParameterLists( 1 ).insert( "variable_5", 2.0 );
        aParameterLists( 1 ).insert( "variable_6", 4.5 * sH );
        aParameterLists( 1 ).insert( "variable_7", 5.5 * sH );
        aParameterLists( 1 ).set( "name", "Thin_Wall_Inner" );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Thin wall outer 5
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Thin_Wall" );
        aParameterLists( 1 ).insert( "variable_1", 6.5 * sH );
        aParameterLists( 1 ).insert( "variable_2", 5.0 * sH );
        aParameterLists( 1 ).insert( "variable_3", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_4", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_5", 2.0 );
        aParameterLists( 1 ).insert( "variable_6", 2.5 * sH );
        aParameterLists( 1 ).insert( "variable_7", 7.5 * sH );
        aParameterLists( 1 ).set( "name", "Thin_Wall_Outer" );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );

        // Inclusions 6
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Func_Inclusion" );
        aParameterLists( 1 ).set( "name", "Level_Set_Field" );
        aParameterLists( 1 ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0 );
        if ( tIsOpt )
        {
            aParameterLists( 1 ).set( "discretization_mesh_index", 1 );
            aParameterLists( 1 ).set( "discretization_lower_bound", -tBSplineLimit );
            aParameterLists( 1 ).set( "discretization_upper_bound", tBSplineLimit );
        }

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gen_property_parameter_list( gen::Field_Type::SCALED_FIELD ) );
        uint tParamCounter = 0;

        aParameterLists( 2 ).set( "name", "LevelsetField" );
        aParameterLists( 2 ).set( "dependencies", "Level_Set_Field" );
        aParameterLists( 2 ).set( "scaling_factor", 1.0 );
        aParameterLists( 2 ).set( "pdv_type", "LS1" );
        aParameterLists( 2 ).set( "pdv_mesh_set_names", tFluid + "," + tInclusionSolid );
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

        // create a cell of cell of parameter list for fem
        uint tPropIndex  = 0;
        uint tCMIndex    = 1;
        uint tSPIndex    = 2;
        uint tIWGIndex   = 3;
        uint tIQIIndex   = 4;
        uint tFEMIndex   = 5;
        uint tPhaseIndex = 7;

        //------------------------------------------------------------------------------

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseVoid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseInclusion" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "2" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseAll" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1,2" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", "PhaseHeatMethod" );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1,2" );

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // fluid properties ------------------------------------------------------------
        // create parameter list for property 2
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidDensity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidDensity );

        // create parameter list for property 1
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidViscosity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidDynViscosity );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Ramp_Property" );

        // create parameter list for property 1
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidKinViscosity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidKinViscosity );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Ramp_Property" );

        // create fluid pressure spring property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropFluidPressureSpring" );
        aParameterLists( tPropIndex ).set( "function_parameters", tFluidPressureSpring );

        // BC properties ---------------------------------------------------------------
        // create inlet velocity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletU" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Inlet_U" );
        aParameterLists( tPropIndex ).set( "function_parameters", tInletVelocity );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSelectInletU" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Select_Inlet_U" );

        // create zero velocity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDirichletZeroU" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // create wall distance property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWallDistance" );
        aParameterLists( tPropIndex ).set( "dof_dependencies", "L2" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Wall_Distance" );
        aParameterLists( tPropIndex ).set( "dof_derivative_functions", "Func_Wall_Distance_Der" );

        // create wall distance property
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInvWallDistanceSquare" );
        aParameterLists( tPropIndex ).set( "dof_dependencies", "PHID" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Wall_InvDistanceSquare" );
        aParameterLists( tPropIndex ).set( "dof_derivative_functions", "Func_Wall_InvDistanceSquare_Der" );

        // inlet viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletV" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Inlet_V" );
        aParameterLists( tPropIndex ).set( "function_parameters", tInletKinViscosity );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSelectInletV" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Select_Inlet_V" );

        // zero viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDirichletZeroV" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0" );

        // upwind
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInletUpwind" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Inlet_Upwind" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // time continuity weights for velocity and viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightUV" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Time_Weights" );

        // time continuity weights for velocity and viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightResUV" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Time_Weight_Res" );

        // initial condition velocity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInitialConditionU" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0;0.0" );

        // initial condition viscosity
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInitialConditionV" );
        aParameterLists( tPropIndex ).set( "function_parameters", tInletKinViscosity );

        // select integration domain
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropSelectDomain" );
        aParameterLists( tPropIndex ).set( "value_function", "Func_Select_Int_Domain" );

        // Wall distance properties ----------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropConductivity" );
        aParameterLists( tPropIndex ).set( "function_parameters", tConductivity );

        // properties for Theta

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensityTheta" );
        aParameterLists( tPropIndex ).set( "function_parameters", tDensityTheta );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropCapacityTheta" );
        aParameterLists( tPropIndex ).set( "function_parameters", tCapacityTheta );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPrescTheta" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPrescTheta );

        // properties for phi problem

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropDensityPhi" );
        aParameterLists( tPropIndex ).set( "function_parameters", tDensityPhi );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropCapacityPhi" );
        aParameterLists( tPropIndex ).set( "function_parameters", tCapacityPhi );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropPrescPhi" );
        aParameterLists( tPropIndex ).set( "function_parameters", tPrescPhi );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropEigenStrainPhi" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // time continuity weights
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( tPropIndex ).set( "function_parameters", "10.0" );

        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( tPropIndex ).set( "function_parameters", "10.0" );

        // initial condition
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropInitialCondition" );
        aParameterLists( tPropIndex ).set( "function_parameters", "0.0" );

        // target level set
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLevelSetConst" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );

        // target level set gradient
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists( tPropIndex ).set( "function_parameters", tLevelSetGradxConstant );

        // actual level set
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLevelSet" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );
        aParameterLists( tPropIndex ).set( "value_function", "tLevelSetFunc" );
        aParameterLists( tPropIndex ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists( tPropIndex ).set( "dv_dependencies", "LS1" );

        // actual level set gradient
        aParameterLists( tPropIndex ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( tPropIndex ).set( "property_name", "PropLevelSetGradx" );
        aParameterLists( tPropIndex ).set( "function_parameters", "1.0" );
        aParameterLists( tPropIndex ).set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists( tPropIndex ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists( tPropIndex ).set( "dv_dependencies", "LS1" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // incompressible NS fluid constitutive model
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluid" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::FLUID_TURBULENCE );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P;VISCOSITY", "Velocity,Pressure,Viscosity" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropFluidViscosity   ,Viscosity;"
                "PropFluidKinViscosity,KinViscosity;"
                "PropFluidDensity     ,Density" );

        // Spalart Allmaras turbulence constitutive model
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMTurbulence" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        aParameterLists( tCMIndex ).set( "function_parameters", tCMTurbFt2 + "/" + tCMTurbAlpha );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropFluidKinViscosity,KinViscosity;"
                "PropWallDistance     ,WallDistance" );

        //------------------------------------------------------------------------------
        // create parameter list for constitutive model - Theta problem
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluidDiffusionTheta" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMInclusionSolidDiffusionTheta" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseInclusion" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMFluidDiffusionPhi" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseFluid" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( tCMIndex ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( tCMIndex ).set( "constitutive_name", "CMInclusionSolidDiffusionPhi" );
        aParameterLists( tCMIndex ).set( "phase_name", "PhaseInclusion" );
        aParameterLists( tCMIndex ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( tCMIndex ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists( tCMIndex ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // SP for fluid phase ----------------------------------------------------------
        // SUPG_PSPG NS incompressible
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPIncFlow" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        aParameterLists( tSPIndex ).set( "function_parameters", tSupgFluidC1 + "/0.0" );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // Dirichlet Nitsche for fluid velocity
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPFluidDirichletNitscheU" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", sIfcNitscheFluid );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists( tSPIndex ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );

        // SUPG Spalart-Allmaras
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPSUPGSA" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
        aParameterLists( tSPIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists( tSPIndex ).set( "function_parameters", tSupgTurbPower + "/" + tSupgTurbSource );
        aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );

        // Dirichlet Nitsche for fluid turbulent viscosity
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheV" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", sIfcNitscheViscosity );
        aParameterLists( tSPIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

        // SP for ghost ----------------------------------------------------------------
        if ( tUseGhost )
        {
            // ghost fluid viscous
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPViscous" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
            aParameterLists( tSPIndex ).set( "function_parameters", alpha_velo );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidViscosity,Viscosity" );

            // ghost fluid convective
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPVelocity" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
            aParameterLists( tSPIndex ).set( "function_parameters", alpha_conv );
            aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidDensity,Density" );

            // ghost fluid pressure
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPPressure" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
            aParameterLists( tSPIndex ).set( "function_parameters", alpha_pres );
            aParameterLists( tSPIndex ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            aParameterLists( tSPIndex ).set( "leader_properties",
                    "PropFluidViscosity,Viscosity;"
                    "PropFluidDensity,Density" );

            // ghost fluid turbulence viscosity
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPViscosity" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists( tSPIndex ).set( "function_parameters", alpha_visc );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropFluidKinViscosity,Material" );

            // create parameter list for projection of the distance
            aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPFluidL2" );
            aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists( tSPIndex ).set( "function_parameters", alpha_l2 );
            aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );
        }

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPFluidThetaPhi" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", alpha_heat );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheFluidThetaPhi" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", sIfcNitscheHeat );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPGPInclusionThetaPhi" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tSPIndex ).set( "follower_phase_name", "PhaseInclusion" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( tSPIndex ).set( "function_parameters", alpha_heat );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( tSPIndex ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( tSPIndex ).set( "stabilization_name", "SPNitscheInclusionThetaPhi" );
        aParameterLists( tSPIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tSPIndex ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( tSPIndex ).set( "function_parameters", sIfcNitscheHeat );
        aParameterLists( tSPIndex ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // fluid bulk IWG --------------------------------------------------------------
        // NS incompressible for velocity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidVelocityBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // NS incompressible for pressure
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidPressureBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );

        // turbulence
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGTurbulenceBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPSUPGSA,SUPG" );

        // L2 projection for distance
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDistance" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::L2 );
        aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInvWallDistanceSquare,Source" );

        // fluid BC IWG ----------------------------------------------------------------
        // incompressible NS velocity Dirichlet IWG for inlet
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletVelocity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select;PropInletUpwind,Upwind" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );

        // incompressible NS pressure Dirichlet IWG for inlet
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletPressure" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        // zero velocity for velocity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidZeroVelocity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );

        // zero velocity for pressure
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidZeroPressure" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );

        if ( tUsePseudoTimeStepping )
        {
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
            aParameterLists( tIWGIndex ).set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionU,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual" );
            }

        // inlet viscosity
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInletViscosity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "side_ordinals", "4" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropInletV,Dirichlet;PropSelectInletV,Select;PropInletUpwind,Upwind" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        // zero viscosity for walls
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGZeroViscosity" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropDirichletZeroV,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );

        if ( tUsePseudoTimeStepping )
        {
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
            aParameterLists( tIWGIndex ).set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionV,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual" );
            }

        if ( tUseGhost )
        {
            // ghost fluid viscous
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPViscous" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );

            // ghost fluid convective
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPConvective" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VX,VY" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );

            // ghost fluid pressure
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPPressure" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "P" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );

            // ghost fluid viscosity
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPFluidViscosity" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "VISCOSITY" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );

            // ghost L2 projection distance
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPInnerL2" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "L2" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidL2,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // theta problem

        // theta bulk in fluid
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidDiffusionThetaBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );

        // theta bulk in inclusion
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInclusionSolidDiffusionThetaBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPInnerTheta" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPOuterTheta" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseInclusion" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGFluidTimeContinuityTheta" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGInclusionTimeContinuityTheta" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( tIWGIndex ).set( "dof_residual", "THETA" );
        aParameterLists( tIWGIndex ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );

        //------------------------------------------------------------------------------
        // phid problem

        // create IWG - bulk diffusion
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionInnerBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );

        // create IWG - bulk diffusion
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );

        // create parameter list for single side interface condition
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );

        // create parameter list for single side interface condition
        aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( tIWGIndex ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIWGIndex ).set( "neighbor_phases", "PhaseFluid" );
        aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
        aParameterLists( tIWGIndex ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( tIWGIndex ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG - ghost
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPInnerPhi" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseFluid" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );

            // create IWG - ghost
            aParameterLists( tIWGIndex ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( tIWGIndex ).set( "IWG_name", "IWGGPOuterPhi" );
            aParameterLists( tIWGIndex ).set( "leader_phase_name", "PhaseInclusion" );
            aParameterLists( tIWGIndex ).set( "follower_phase_name", "PhaseInclusion" );
            aParameterLists( tIWGIndex ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists( tIWGIndex ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( tIWGIndex ).set( "dof_residual", "PHID" );
            aParameterLists( tIWGIndex ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVX" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVY" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VX,VY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkP" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "P" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "THETA" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkPHID" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // fluid volume
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIFluidVolume" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // inclusion volume
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQISolidVolume" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseInclusion" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // inclusion perimeter
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIPerimeterItf" );
        aParameterLists( tIQIIndex ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "neighbor_phases", "PhaseInclusion" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::VOLUME );

        // fluid power dissipation in volume
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIVolumePowDisp" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION_BULK );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropSelectDomain,Select" );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid" );

        // level set
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQILevelSet" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseHeatMethod" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropLevelSet,Property" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIHeatMethodPenaltyFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIHeatMethodPenaltyInclusion" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseInclusion" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIProjectedLevelSetFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 3 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIProjectedLevelSetInclusion" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 3 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseInclusion" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIHeatAlphaFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 5 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );

        // heat method penalty
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIHeatAlphaInclusion" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "PHID" );
        aParameterLists( tIQIIndex ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 5 );
        aParameterLists( tIQIIndex ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        aParameterLists( tIQIIndex ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseInclusion" );

        // viscosity dof
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkVISCOSITY" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "VISCOSITY" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // L2 projection distance dof
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkL2" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( tIQIIndex ).set( "dof_quantity", "L2" );
        aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );

        // wall distance
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIWallDistance" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( tIQIIndex ).set( "leader_properties", "PropWallDistance,Property" );

        // turbulent dynamic viscosity
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkTurbDynVisc" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );

        // effective dynamic viscosity
        aParameterLists( tIQIIndex ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkEffDynVisc" );
        aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
        aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY );
        aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );
        /*
                // SUPG fluid
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkSUPGU" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      0 );
                aParameterLists( tIQIIndex ).set( "stabilization_parameters",   "SPIncFlow,Stabilization" );

                // PSPG fluid
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkPSPGU" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      1 );
                aParameterLists( tIQIIndex ).set( "stabilization_parameters",   "SPIncFlow,Stabilization" );

                // SUPG viscosity
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkSUPGSA" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      0 );
                aParameterLists( tIQIIndex ).set( "stabilization_parameters",   "SPSUPGSA,Stabilization" );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UX" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      1 );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UY" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      2 );

                // incompressible NS strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_P" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
                aParameterLists( tIQIIndex ).set( "vectorial_field_index",      0 );

                // SA strong residual
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIBulkSpalartAllmarasSF" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_SA );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence") ;
                        // production coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkProductionCoeff" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index", 0 );
                aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // wall destruction coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkWallDestructionCoeff" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index", 1 );
                aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // diffusion coefficient
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkDiffusionCoeff" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index", 2 );
                aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                 // production term
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkProductionTerm" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index", 3 );
                aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );

                // wall destruction term
                aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
                aParameterLists( tIQIIndex ).set( "IQI_name", "IQIBulkWallDestructionTerm" );
                aParameterLists( tIQIIndex ).set( "leader_phase_name", "PhaseFluid" );
                aParameterLists( tIQIIndex ).set( "vectorial_field_index", 4 );
                aParameterLists( tIQIIndex ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
                aParameterLists( tIQIIndex ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
                */
        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( tFEMIndex ).add_parameter_list( prm::create_computation_parameter_list() );
        aParameterLists( tFEMIndex ).set( "print_physics_model", false );
        aParameterLists( tFEMIndex ).set( "is_analytical_sensitivity", false );
        aParameterLists( tFEMIndex ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists( tFEMIndex ).set( "finite_difference_perturbation_size", stod( tFDEpsilon ) );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
#ifdef MORIS_USE_MUMPS
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        // Newton solve in forward problem used by Newton flow subproblem (includes load control)
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_ref_iter", tNLA_ref_its );

        if ( tSolverType == "Newton" && tUseNewtonLoadControl )
        {
            aParameterLists( 2 ).set( "NLA_load_control_strategy", sol::SolverLoadControlType::Exponential );
            aParameterLists( 2 ).set( "NLA_load_control_factor", tRampInitial );
            aParameterLists( 2 ).set( "NLA_load_control_steps", tRampSteps );
            aParameterLists( 2 ).set( "NLA_load_control_relres", tRampRelRes );
            aParameterLists( 2 ).set( "NLA_load_control_exponent", 1.0 );
        }

        // for linear solve - forward: theta, phi, l2 - adjoint: all problems
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        // NLBGS for overall forward solve and adjoint in theta,phi,l2,(vx,vy,p,viscosity)
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        // algorithm not used
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        // NLBGS for flow subproblem using pseudo time stepping and (optional) load control)
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLBGS_rel_res );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLBGS_realx );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLBGS_max_itr );
        aParameterLists( 2 ).set( "NLA_ref_iter", tNLBGS_ref_its );

        if ( tSolverType == "NLBGS" )
        {
            if ( tUsePseudoTimeStepping )
            {
                aParameterLists( 2 ).set( "NLA_pseudo_time_control_strategy", sol::SolverPseudoTimeControlType::Expur );
                aParameterLists( 2 ).set( "NLA_pseudo_time_initial_steps", tNLBGS_init_itr );
                aParameterLists( 2 ).set( "NLA_pseudo_time_initial", tSteadyStateStepSize );
                aParameterLists( 2 ).set( "NLA_pseudo_time_max_num_steps", tMaxNumTimeSteps );
                aParameterLists( 2 ).set( "NLA_pseudo_time_max_step_size", tMaxTimeStepSize );
                aParameterLists( 2 ).set( "NLA_pseudo_time_rel_res_norm_drop", tRelResNormDrop );
                aParameterLists( 2 ).set( "NLA_pseudo_time_rel_res_norm_update", tRelResNormUpdate );
                aParameterLists( 2 ).set( "NLA_pseudo_time_steady_rel_res_norm", tSteadyStateRelRes );
                aParameterLists( 2 ).set( "NLA_pseudo_time_steady_step_size", tSteadyStateStepSize );
                aParameterLists( 2 ).set( "NLA_pseudo_time_offset", tTimeOffSet );
                aParameterLists( 2 ).set( "NLA_pseudo_time_constant", tConstantTimeStep );
                aParameterLists( 2 ).set( "NLA_pseudo_time_step_index_factor", tIndexFactor );
                aParameterLists( 2 ).set( "NLA_pseudo_time_step_index_exponent", tIndexExpo );
                aParameterLists( 2 ).set( "NLA_pseudo_time_residual_factor", tResFactor );
                aParameterLists( 2 ).set( "NLA_pseudo_time_residual_exponent", tResExpo );
                aParameterLists( 2 ).set( "NLA_pseudo_time_comsol_1", tComsolParameter1 );
                aParameterLists( 2 ).set( "NLA_pseudo_time_comsol_2", tComsolParameter2 );
            }
            if ( tUseNLBGSLoadControl )
            {
                aParameterLists( 2 ).set( "NLA_load_control_strategy", sol::SolverLoadControlType::Exponential );
                aParameterLists( 2 ).set( "NLA_load_control_factor", tRampInitial );
                aParameterLists( 2 ).set( "NLA_load_control_steps", tRampSteps );
                aParameterLists( 2 ).set( "NLA_load_control_relres", tRampRelRes );
                aParameterLists( 2 ).set( "NLA_load_control_exponent", 1.0 );
            }
        }

        // Newton solve used within NLBGS for flow subproblem
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNewton_rel_res );
        aParameterLists( 2 ).set( "NLA_max_iter", tNewton_max_iter );
        aParameterLists( 2 ).set( "NLA_relaxation_strategy", sol::SolverRelaxationType::InvResNormAdaptive );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_relaxation_damping", 0.5 );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "PHID" );

        if ( tSolverType == "NLBGS" )
        {
            if ( tFluidViscSolver == "Mono" )    // monolithic flow problem
            {
                aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
                aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
                aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
            }
            else    // flow problem in vx,vy,p
            {
                aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
                aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
                aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P" );
            }
        }
        else    // using Newton solver only monolithic
        {
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }

        // flow viscosity problem only
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VISCOSITY" );

        // following settings only used when solving flow problem by NLBGS
        if ( tFluidViscSolver == "Mono" )
        {
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "4" );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "2" );
            aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }
        else
        {
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "4" );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "2,3" );
            aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P;VISCOSITY" );
        }

        // L2 projection problem
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "L2" );

        // total forward problem
        if ( tSolverType == "NLBGS" )
        {
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,4" );
            aParameterLists( 3 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }
        else
        {
            aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
            aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,2" );
            aParameterLists( 3 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }

        // adjoint of flow problem
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        // adjoint of total problem
        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,7" );    // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );

        // ----------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 6 );                // using NLBGS for forward problem
        aParameterLists( 4 ).set( "TSA_Nonlinear_Sensitivity_Solver", 8 );    // using NLBGS for sensitivity problem

        //------------------------------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0;VX," + tVXInitial + ";VY,0.0;P,0.0;VISCOSITY," + tInletKinViscosity + ";L2,1.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        //------------------------------------------------------------------------------

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "THETA", 0 );
        aParameterLists( 0 ).set( "PHID", 0 );
        aParameterLists( 0 ).set( "VX", 0 );
        aParameterLists( 0 ).set( "VY", 0 );
        aParameterLists( 0 ).set( "P", 0 );
        aParameterLists( 0 ).set( "VISCOSITY", 0 );
        aParameterLists( 0 ).set( "L2", 0 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tAllDomain );
        aParameterLists( 0 ).set( "Field_Names",
                "THETA,PHID,L2,"
                "VX,VY,P,VISCOSITY,"
                "TURBDYNVISC,EFFVISC,"
                "LS,WALLDIST,PRJLEVSFLD,PRJLEVSINC,ALPHAFLD,ALPHAINC" );
        aParameterLists( 0 ).set( "Field_Type",
                "NODAL,NODAL,NODAL,"
                "NODAL,NODAL,NODAL,NODAL,"
                "NODAL,NODAL,"
                "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists( 0 ).set( "IQI_Names",
                "IQIBulkTHETA,IQIBulkPHID,IQIBulkL2,"
                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkVISCOSITY,"
                "IQIBulkTurbDynVisc,IQIBulkEffDynVisc,"
                "IQILevelSet,IQIWallDistance,IQIProjectedLevelSetFluid,IQIProjectedLevelSetInclusion,IQIHeatAlphaFluid,IQIHeatAlphaInclusion" );
        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
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
