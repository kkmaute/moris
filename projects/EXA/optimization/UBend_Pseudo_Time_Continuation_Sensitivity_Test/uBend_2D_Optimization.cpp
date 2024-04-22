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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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
    real tDensFld = 1.0;          // kg/m3
    real tKinVisc = 4.0e-5;       // m2/s
    
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
    
    bool tPowerlawProfile = false; // FIXME to be changed to power law for turbulence profile
    
    bool tRampViscBC      = false;
    
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
    std::string tSolverType            = "NLBGS";       // "Newton" or "NLBGS";
    std::string tFluidViscSolver       = "Staggered";   // "Staggered" or "Mono"; only monolithic scheme possible with "Newton" algorithm
    std::string tSubNewtonSolverOption = "5";           // "0": standard Newton solver  "5": Newton with adaptive relaxation scheme
     
    bool tUsePseudoTimeStepping = true;                 // NLBGS with pseudo time stepping
    
    // Ramping of flow properties
    sint tPropRampSteps     =  3;       // number of NLBFGS steps for ramping properties (>=1); here dynamic viscosity; nominal value reached in iteration tPropRampSteps+1
    real tPropRampScaling   = 100.0;   // factor by which dynamic viscosity is increased
    
    // NLBGS settings for fluid subproblem
    sint tNLBGS_max_itr  = 35;                             // max number of step in NLBGS for flow subproblem
    real tNLBGS_rel_res  = 1.0e-9;                         // required convergence on pseudo-dynamic residual
    real tNLBGS_realx    = 1.0;                            // relaxation parameter
    sint tNLBGS_init_itr = tPropRampSteps;                 // initialization phase during which property has not reached nominal value
    sint tNLBGS_ref_its  = std::max(tPropRampSteps+1,2);   // iteration id when reference residual to be calculated
                                                           // determines minimum number of steps in NLBFGS; at least two iteration need to be performed
    
    // NLBGS - Load control
    bool tUseNLBGSLoadControl = false;

    // Newton - Load control
    bool tUseNewtonLoadControl = false;
    
    // Load control settings used for Newton and NLBGS
    real tRampInitial = 0.01;        // initial load factor
    real tRampRelRes  = 0.01;        // max rel. residual below which load factor is increased 
    sint tRampSteps   = 3.;          // number of ramping steps
    
    // NLBGS - Pseudo time stepping
     bool tUseLoadTimeRamp       = false;
    bool tUseGlobalTimeStep     = false;
    
    real tPseudoTimeIndexCflOld   = -1;
    real tPseudoTimeIndexPropOld  = -1;

    sint tMaxNumTimeSteps     =  35;        // max. number of steps in PTC
    real tMaxTimeStepSize     =  1.0e6;     // max. CFL
    real tRampTotalTime       =  1.0;       // total CFL corresponding to end of load ramping (only used if tUseLoadTimeRamp = true)
    real tRelResNormDrop      =  1.0e-9;    // required rel. static res drop for convergence 
    real tRelResNormUpdate    =  1.0e1;     // required rel. static res drop for stepping forward in time
    real tRelResInexNewton    = -1.0;       // required rel. static res drop for switching to inexact Newton
    real tSteadyStateRelRes   = -1.0;       // required rel. static res drop for switching to exact Newton
    real tSteadyStateStepSize =  1.0e6;     // CFL in exact Newton
    real tTimeOffSet          =  0.0;       // time offset for writing PTC steps (= 0: no output)
    
    real tConstantTimeStep = 1.0;
    real tIndexFactor      = 2.0;
    real tIndexExpo        = 1.0;
    
    real tResFactor        = 1.8;    // for Expur: increase factor
    real tResExpo          = 0.33;   // for Expur: decrease factor
    
    real tComsolParameter1  = 20.0;
    real tComsolParameter2  = 30.0;
    
    // Newton paramters when using NLBGS
    real tNewton_rel_res   = 2.5e-1;
    real tNewton_relax     = 1.0;
    sint tNewton_max_iter  = 15;  
    
    // Newton paramters without NLBGS
    moris::real tNLA_rel_res_norm_drop    = 1.0e-9;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 100;  
    sint        tNLA_ref_its              = tSolverType == "Newton" ? std::max(tPropRampSteps+1,2) : 1;

    // Time solver parameters
    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0;
    
    std::string tVXInitial = tUseNLBGSLoadControl ? std::to_string(tRampInitial) : tUseLoadTimeRamp ? std::to_string(1.0/tRampTotalTime) : "0.001";
  
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
    moris::real tOffsetInOutlet  = 2.0; // original
    moris::real tOffsetSide = 1.0;
    moris::real tOffsetCut = 0.0214;

    moris::real tDimX = ( tOffsetInOutlet + 12.0 + tOffsetSide ) * sH + tOffsetCut;    // x-dimension of computational domain
    moris::real tDimY = ( tOffsetSide + 10.0 + tOffsetSide ) * sH + tOffsetCut;    // y-dimension of computational domain

    moris::real tOffsetX = -tOffsetInOutlet * sH - tOffsetCut;
    moris::real tOffsetY = -tOffsetSide     * sH - tOffsetCut;

    moris::real tApproxEleSize = 1.0 / ( tNumElements  ) * tLengthScale;

    /* ------------------------------------------------------------------------ */
    // background mesh
    std::string tNumElemX = moris_to_string( std::ceil( tDimX / tApproxEleSize ) + 1);
    std::string tNumElemY = moris_to_string( std::ceil( tDimY / tApproxEleSize ) + 1);

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

    std::string tInterfaceRefinementInclusion = "0";

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

    moris::real sWeightPowDis = 1.0;     // weight on temperature obj
    moris::real sWeightPerim  = 0.0; //0.1;    // weight on perimeter
    moris::real sWeightHM     = 0.0; //0.5;    // weight on heat method
    
    moris::real sFldVolPenalty = 1.0;   // constraint contribution to objective

    // formulation of constraint
    moris::real sWeightFluidV   = 1.0;      // weight on max fluid volume
    moris::real sFractionFluidV = 0.3;      // weight on max fluid volume

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
    std::string tFDEpsilon            = "1e-5";
    std::string tFDSweep              = "1e-5";
    
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
    real tDynVis  = tDensFld*tKinVisc;

    // reynolds number based on reference values and fluid properties see below
    real tDimReynolds = lenref * velref / tKinVisc;

    std::string tFluidDensity      = moris_to_string( tDensFld * tDensityScale );                              // kg/m3   air at 300 K and 1 atm
    std::string tFluidDynViscosity = moris_to_string( tDynVis  * tPressureScale * tTimeScale );                // N s/m2  air at 300 K
    std::string tFluidKinViscosity = moris_to_string( tKinVisc * tLengthScale * tLengthScale / tTimeScale );   // m2/s air at 300 K

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if( tY >= 5.5 * sH && tY <=7.5*sH )
        {
            aPropMatrix( 0, 0 ) = 1.0;
            aPropMatrix( 1, 1 ) = 1.0;
        }
    }
    
    /* ------------------------------------------------------------------------ */
    void
    Func_Inlet_Upwind(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );
        
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if( tY >= 5.5 * sH && tY <=7.5*sH )
        {
            aPropMatrix( 0, 0 ) = aParameters( 0 )( 0 );
        }
    }
    
    /* ------------------------------------------------------------------------ */
    void
    Func_Select_Inlet_V(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1, 0.0 );
        
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        if( tY >= 5.5 * sH && tY <=7.5*sH )
        {
            aPropMatrix( 0, 0 ) = 1.0;
        }
    }

    /* ------------------------------------------------------------------------ */
    // pseudo time stepping function
    void
    Func_Time_Weights(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1 );
        
        if ( gLogger.exists("MDL","Model","Perform Forward Analysis") )
        {
            real tWeight = aParameters( 0 )( 0 );
           
            real tCFL = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "PseudoTimeStep" );
    
            const Matrix<DDRMat> & tXhat = aFIManager->get_IP_geometry_interpolator()->get_space_coeff();
          
            real tXmax = tXhat.get_column(0).max();
            real tXmin = tXhat.get_column(0).min();
          
            real tEleLength = tXmax - tXmin;

            real tVelocNorm = 1.0;           

            if ( ! tUseGlobalTimeStep )
            {
                moris::fem::Field_Interpolator_Manager* tPrevFIManager = aFIManager->get_field_interpolator_manager_previous();

                const Matrix<DDRMat> & tVelocity = tPrevFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX )->val();
           
                tVelocNorm =  norm( tVelocity );
            }
    
            aPropMatrix( 0 ) = tWeight * tVelocNorm /std::max(MORIS_REAL_EPS,tCFL * tEleLength);
	    
            real tPseudoTimeStep = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "PseudoTimeStep" );
	    
            if ( tPseudoTimeStep < 1.0) 
            {
                aPropMatrix( 0 ) = 0.0;
            }
            
          // get pseudo time indes of most recent NLBGS iteration (last parameter set to true)
           uint tPseudoTimeIndex = gLogger.get_iteration( "NonLinearAlgorithm", "NLBGS", "Solve", true );

            
           if ( std::abs(tPseudoTimeIndex-tPseudoTimeIndexCflOld)>1e-8)
           {
               MORIS_LOG_INFO("PseudoTimeIndex = %d  in Func_Time_Weights: tPseudoTimeStep = %f  tCFL = %e",tPseudoTimeIndex,tPseudoTimeStep,tCFL);
               
               tPseudoTimeIndexCflOld = tPseudoTimeIndex; 
           }
        }
        else if( gLogger.exists("MDL","Model","Perform Sensitivity Analysis") )
        {
            aPropMatrix( 0 ) = 0.0;
        }
        else
        {
            MORIS_ERROR( false, "Func_Time_Weights - Neither FA, nor SA.");
            aPropMatrix( 0 ) = 0.0;
        }
    } 
    
    /* ------------------------------------------------------------------------ */
    // To turn on/off time continuity residual for inexact Newton iterations
    // in this case time continuitiy residual is omitted but jacobian is still
    // computed
    void
    Func_Time_Weight_Res(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
       aPropMatrix.set_size( 1, 1 );
        
        if ( gLogger.exists("MDL","Model","Perform Forward Analysis") )
        {
           real tStaticRelRes = gLogger.get_action_data( "NonLinearAlgorithm", "NLBGS", "Solve", "RelativeStaticResidual" );
           aPropMatrix( 0 ) = tStaticRelRes < tRelResInexNewton ? 0.0 : 1.0;
        }
        else if( gLogger.exists("MDL","Model","Perform Sensitivity Analysis") )
        {
            aPropMatrix( 0 ) = 0.0;
        }
        else
        {
            MORIS_ERROR( false, "Func_Time_Weights - Neither FA, nor SA.");
            aPropMatrix( 0 ) = 0.0;
        }
    }
    
    /* ------------------------------------------------------------------------ */
    // To decrease viscosity in initialization phase of NLBGS
    void
    Func_Ramp_Property(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
       real tProperty = aParameters( 0 )( 0 ); 
       
       if ( gLogger.exists("MDL","Model","Perform Forward Analysis") )
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
               real tFactor = std::max(0.0,((real)tPropRampSteps-(real)tPseudoTimeIndex+1.0)/((real)tPropRampSteps));

               if ( tLogInterpolation )
               {
                   aPropMatrix = std::pow(10.0,tFactor*std::log10(tMaxProperty)+(1.0-tFactor)*std::log10(tProperty));
               }
               else
               {
                   aPropMatrix = tFactor*tMaxProperty+(1.0-tFactor)*tProperty;
               }
               if ( std::abs(tPseudoTimeIndex-tPseudoTimeIndexPropOld)>1e-8)
               {
                   MORIS_LOG_INFO("PseudoTimeIndex = %d  factor = %e  nominal Property = %e  actual Property = %e",tPseudoTimeIndex,tFactor,tProperty,aPropMatrix(0));
               
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
    	// init property container
        aPropMatrix.set_size( 1, 1 , 0.0 );

    	// grab x-coordinate
    	real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
	
	// if within selected domain in x
	if( tX >= 2.0 && tX <= 12.0*sH )
	{
    	    // grab y-coordinate
	    real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            // if within selected domain in y
	    if( tY >= 0.0 && tY<= 10.0*sH )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
                real tPower = 2.0 * std::log10( tDimReynolds / 10.0 );
                aPropMatrix( 0 ) = tScalingFactor * tInVelocity * ( std::pow( 1.0 - ( ( tY - 6.5 * sH ) / std::abs( tY - 6.5 * sH ) ) * ( ( tY - 6.5 * sH ) / sH ) , 1.0 / tPower ) );
            }
            else
            {
                aPropMatrix( 0 ) = tScalingFactor * tInVelocity *  ( 1.0 - std::pow( ( tY - 6.5 * sH ) / sH, 2.0 ) );
            }
        }
    }

    /* ------------------------------------------------------------------------ */
    // inlet viscosity function
    void
    Func_Inlet_V(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 1, 1 , 0.0 );
  
        moris::real tInVisc = aParameters( 0 )( 0 );

        real tScalingFactor = tRampViscBC ? std::min( 1.0, gLogger.get_action_data( "NonLinearAlgorithm", tSolverType, "Solve", "LoadFactor" ) ) : 1.0;
  
        aPropMatrix( 0 ) = tScalingFactor * tInVisc;
    }

    /* ------------------------------------------------------------------------ */
    // wall distance function
    void
    Func_Wall_Distance(
            moris::Matrix< moris::DDRMat >&    aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tL2Val = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->val()( 0 );
        aPropMatrix.set_size( 1, 1, std::max( std::sqrt( std::pow( 10.0, - tL2Val ) ), 1e-9 ) );
    }

    /* ------------------------------------------------------------------------ */
    // wall distance derivative function
    void
    Func_Wall_Distance_Der( 
            moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&  aParameters,
            moris::fem::Field_Interpolator_Manager*         aFIManager )
    {
        real tL2Val = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->val()( 0 );

        real tFuncVal = std::sqrt( std::pow( 10.0, - tL2Val ) );

        if( tFuncVal > 1e-9 )
        {
            aPropMatrix = -log(10) * std::pow(10.0,-tL2Val) / 2.0 / tFuncVal * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::L2 )->N();
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
            moris::Matrix< moris::DDRMat >&    aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tPhiDVal = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );
        
        aPropMatrix.set_size( 1, 1, std::log10( 1.0 / std::pow( std::max( 1e-9, tPhiDVal ), 2.0 ) ) );
    }
    
    // inverse wall distance squared function derivative
    void
    Func_Wall_InvDistanceSquare_Der( 
            moris::Matrix< moris::DDRMat >&    aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tPhiDVal = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );
        if( tPhiDVal > 1e-9 )
        {
            aPropMatrix = - 2.0  / tPhiDVal / log( 10.0 ) * aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->N();
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
    moris::real tInclusionRadius = sH / 1.0;       // inclusion radius
    moris::real tLeftBound       = 2.0  * sH + sH; //2.0 * tInclusionRadius;     // min x coord for bounding box
    moris::real tRightBound      = 12.0 * sH - sH; //2.0 * tInclusionRadius;     // max x coord for bounding box
    moris::real tBottomBound     = 0.0  * sH + sH; //2.0 * tInclusionRadius;     // min y coord for bounding box
    moris::real tTopBound        = 10.0 * sH - sH; //2.0 * tInclusionRadius;     // max y coord for bounding box

    uint tNumXRods = 3;    // number of inclusions along x
    uint tNumYRods = 3;    // number of inclusions along y

    moris::real tInclusionExponent = 12.0;    // inclusion exponent

    moris::real
    Func_Inclusion(
            const moris::Matrix< DDRMat >&    aCoordinates,
            const Vector< moris::real >& aGeometryParameters )
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

            //if ( fmod( iYRod, 2 ) != 0 )
            //{
            //    tXOffset = ( tRightBound - tLeftBound ) / ( 2 * ( tNumXRods - 1 ) );
            //    tTempNumXRods--;
            //}

            for ( uint iXRod = 0; iXRod < tTempNumXRods; iXRod++ )
            // for( uint iXRod = 0; iXRod < tNumXRods; iXRod++ )
            {
                    real tXCenter = tLeftBound + tXOffset + iXRod * ( tRightBound - tLeftBound ) / ( tNumXRods - 1 );
                    // real tXCenter = tLeftBound + iXRod * ( tRightBound - tLeftBound ) / ( tNumXRods - 1 );

                // skip area near the thin wall
                if ( tY > 6.0 * sH || tY < 4.0 *sH || true)
                {
                    real tCylValue = - tInclusionRadius 
                                     + std::pow( std::pow( tX - tXCenter, tInclusionExponent ) 
                                     + std::pow( tY - tYCenter, tInclusionExponent ), 1.0 / tInclusionExponent );

                    tReturnValue = std::min( tReturnValue, tCylValue );
                }
                else
                {
                    // skip area near the thin wall
                    if ( tX > 2.5 * sH ) // 8.5 * sH )
                    {
                        real tCylValue = - tInclusionRadius 
                                         + std::pow( std::pow( tX - tXCenter, tInclusionExponent ) 
                                         + std::pow( tY - tYCenter, tInclusionExponent ), 1.0 / tInclusionExponent );

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
            const moris::Matrix< DDRMat >&    aCoordinates,
            const Vector< moris::real >& aGeometryParameters )
    {
        // get coordinates
        real tX = aCoordinates( 0 );
        real tY = aCoordinates( 1 );
        
        // parameters
        real tXctr = aGeometryParameters( 0 ); //6.5 * sH;
        real tYctr = aGeometryParameters( 1 ); //5.0 * sH;
        real tXRad = aGeometryParameters( 2 ); //sH;
        real tYRad = aGeometryParameters( 3 ); //sH;
        real tExpn = aGeometryParameters( 4 ); //2.0;
        real tBot  = aGeometryParameters( 5 ); //4.5 * sH;
        real tTop  = aGeometryParameters( 6 ); //5.5 * sH;

        // init return value
        real tReturnValue = 1.0;

        // compute value of 3 LS (2 planes, 1 circle)
        real tLS0Value = - tY + tBot;
        real tLS1Value = tY - tTop;

        // combine the LS
        tReturnValue = std::max( tLS0Value, tLS1Value );
        
        if( tX >= tXctr )
        {
            real tLS2Value = - 1.0 + std::pow( std::pow( ( tX - tXctr ) / tXRad, tExpn ) + std::pow( ( tY - tYctr ) / tYRad, tExpn ), 1.0 / tExpn );
            tReturnValue = std::max( tReturnValue, tLS2Value );
        }
        
        return tReturnValue;
    }


    /* ------------------------------------------------------------------------ */
    // Phase assignement
    
    std::string tGetPhaseIndex = "get_phase_index_sharp";
    
    uint get_phase_index_sharp( const Bitset<9>& aGeometrySigns )
    {
        // Phase Void
        if( !aGeometrySigns.test( 0 ) )
        {
            return 0;
        }
        
        if( aGeometrySigns.test( 3 ) )
        {
            return 0;
        }
        
        if( aGeometrySigns.test( 1 ) )
        {
            return 0;
        }
        
        if( !aGeometrySigns.test( 2 ) &&
             aGeometrySigns.test( 5 ) )
        {
            return 0;
        }
        
        if( !aGeometrySigns.test( 4 ) )
        {
            return 0;
        }
        
        
        if( !aGeometrySigns.test( 2 ) &&
             aGeometrySigns.test( 5 ) )
        {
            return 0;
        }
        
                        
        // Phase Inclusion
        if(  aGeometrySigns.test( 0 ) &&
            !aGeometrySigns.test( 1 ) &&
             aGeometrySigns.test( 2 ) &&
            !aGeometrySigns.test( 3 ) &&
            !aGeometrySigns.test( 6 ) )
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
        tObjectives( 0 ) = scaledPowDis + scaledPerim + scaledHM ;

        std::cout << "Contributions to objective function:\n";

        std::cout << "Power dissipation   = " << aCriteria( 0 ) << " \n";
        std::cout << "Perimeter           = " << aCriteria( 1 ) << " \n";
        std::cout << "HM Penalty          = " << aCriteria( 2 ) + aCriteria( 3 ) << " \n";

        std::cout << "Scaled power dissipation   = " << scaledPowDis << " \n";
        std::cout << "Scaled Perimeter           = " << scaledPerim  << " \n";
        std::cout << "Scaled HM penalty          = " << scaledHM     << " \n";	

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
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).push_back( prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOpt );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tSoFile );

        tParameterlist( 1 ).resize( 0 );

        if ( tCheckSensitivities )
        {
            tParameterlist( 2 ).push_back( moris::prm::create_sweep_parameter_list() );
            tParameterlist( 2 )( 0 ).set( "hdf5_path", tHdf5File );
            tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
            tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
            tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tFDSweep );
            tParameterlist( 2 )( 0 ).set( "finite_difference_adv_indices", tAdvIndices ); 
        }
        else
        {
            tParameterlist( 2 ).push_back( moris::prm::create_gcmma_parameter_list() );
            tParameterlist( 2 )( 0 ).set( "step_size", tMMAStepSize );
            tParameterlist( 2 )( 0 ).set( "penalty", tMMAPenalty );
            tParameterlist( 2 )( 0 ).set( "max_its", tMMAMaxIter );  
        }

        if (tRestartId > 0)
        {
            tParameterlist( 0 )( 0 ).set("restart_file", "ADV_Alg_0_Iter_"+std::to_string(tRestartId)+".hdf5" );
            tParameterlist( 2 )( 0 ).set("restart_index", tRestartId );
        }
    }

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", tLagrangePattern );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRef );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "use_advanced_T_matrix_scheme", 1 );
    }

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", false );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        // tParameterlist( 0 )( 0 ).set( "cleanup_cut_mesh",            true );
    }

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "IQI_types",
                "IQIVolumePowDisp,"
                "IQIPerimeterItf,"
                "IQIHeatMethodPenaltyFluid,"
                "IQIHeatMethodPenaltyInclusion,"
                "IQIFluidVolume" );
        tParameterlist( 0 )( 0 ).set( "number_of_phases", 3 );
        tParameterlist( 0 )( 0 ).set( "phase_function_name", tGetPhaseIndex );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );

        uint tGeoCounter = 0;

        // Plane 0 in y = 0
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "line");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 0.0, 0.0, 0.0, 1.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;
        
        // Plane 1 in y = 10.0*sH
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "line");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 0.0, 10.0*sH, 0.0, 1.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;
        
        // Plane 2 in x = 0.0*sH
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "line");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 2.0*sH,0.0, 1.0, 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;
        
        // Plane 3 in x = 12.0*sH
        tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_type", "line");
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", (12.0*sH), 0.0, 1.0, 0.0 );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;
        
        // Thin wall inner 4
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name",   "Func_Thin_Wall" );
        tParameterlist( 1 )( tGeoCounter ).set( "name",                  "Thin_Wall_Inner" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 6.5*sH, 5.0*sH, sH / 2.0, sH / 2.0, 2.0, 4.5*sH, 5.5*sH ) ;
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;

        // Thin wall outer 5
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name",   "Func_Thin_Wall" );
        tParameterlist( 1 )( tGeoCounter ).set( "name",                  "Thin_Wall_Outer" );
        tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", 6.5*sH, 5.0*sH, 2.5*sH, 2.5*sH, 2.0, 2.5*sH, 7.5*sH ) ;
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        tGeoCounter++;
        
        // Inclusions 6
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Inclusion" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Level_Set_Field" );
        tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", tInterfaceRefinementInclusion );
        tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0" );
        if ( tIsOpt )
        {
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 1 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -tBSplineLimit );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", tBSplineLimit );
        }
        tGeoCounter++;

        tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list() );
        uint tParamCounter = 0;

        tParameterlist( 2 )( tParamCounter ).set( "name", "LevelsetField" );
        tParameterlist( 2 )( tParamCounter ).set( "dependencies", "Level_Set_Field" );
        tParameterlist( 2 )( tParamCounter ).set( "field_type", "scaled_field" );
        tParameterlist( 2 )( tParamCounter ).set( "constant_parameters", 1.0 );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tFluid + "," + tInclusionSolid );
        tParamCounter++;
    }

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        if ( par_rank() == 0 )
        { 
             std::cout << "Reynolds       " << tDimReynolds << "------------------------------------" << std::endl;

            std::cout << "velref         " << velref << "------------------------------------" << std::endl;
            std::cout << "lenref         " << lenref << "------------------------------------" << std::endl;
            std::cout << "rhoref         " << rhoref << "------------------------------------" << std::endl;
            std::cout << "pressref       " << pressref << "------------------------------------" << std::endl;
            std::cout << "timeref        " << timeref << "------------------------------------" << std::endl;
            std::cout << "massref        " << massref << "------------------------------------" << std::endl;

            std::cout << "tLengthScale   " << tLengthScale << "------------------------------------" << std::endl;
            std::cout << "tTimeScale     " << tTimeScale << "------------------------------------" << std::endl;
            std::cout << "tMassScale     " << tMassScale << "------------------------------------" << std::endl;

            std::cout << "tPressureScale " << tPressureScale << "------------------------------------" << std::endl;
            std::cout << "tDensityScale  " << tDensityScale << "------------------------------------" << std::endl;

            std::cout << "tFluidDynViscosity     " << tFluidDynViscosity << "------------------------------------" << std::endl;
            std::cout << "tFluidKinViscosity     " << tFluidKinViscosity << "------------------------------------" << std::endl;
            std::cout << "tFluidDensity          " << tFluidDensity << "------------------------------------" << std::endl;

            std::cout << "tInletKinViscosity       " << tInletKinViscosity << "------------------------------------" << std::endl;
            std::cout << "tInletVelocity           " << tInletVelocity << "------------------------------------" << std::endl;
            std::cout << "tElementEdgeLength       " << tElementEdgeLength << "------------------------------------" << std::endl;
	    
            std::cout << "tBSplineLimit            " << tBSplineLimit <<  " (rel wrt fem element : " << tBSplineLimit/tElementEdgeLength << " )\n";
	        std::cout << "tPhiBandwidth            " << tPhiBandwidth <<  " (rel wrt fem element : " << tPhiBandwidth/tElementEdgeLength << " )\n";
            std::cout << "tPhiGradient             " << tPhiGradient  <<"\n";
            std::cout << "tPhiGamma                " << tPhiGamma     <<"\n";

        }

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
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseVoid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "0" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseInclusion" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "2" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseAll" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1,2" );
        tPhaseCounter++;

        tParameterList( tPhaseIndex ).push_back( prm::create_phase_parameter_list() );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_name", "PhaseHeatMethod" );
        tParameterList( tPhaseIndex )( tPhaseCounter ).set( "phase_indices", "1,2" );
        tPhaseCounter++;

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // init property counter
        uint tPropCounter = 0;

        // fluid properties ------------------------------------------------------------
        // create parameter list for property 2
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidDensity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidDensity );
        tPropCounter++;

        // create parameter list for property 1
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidViscosity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidDynViscosity );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Ramp_Property" );
        tPropCounter++;

        // create parameter list for property 1
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidKinViscosity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidKinViscosity );
	tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Ramp_Property" );
        tPropCounter++;

        // create fluid pressure spring property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropFluidPressureSpring" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tFluidPressureSpring );
        tPropCounter++;

        // BC properties ---------------------------------------------------------------
        // create inlet velocity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletU" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Inlet_U" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tInletVelocity );
        tPropCounter++;
        
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSelectInletU" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Select_Inlet_U" );
        tPropCounter++;

        // create zero velocity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDirichletZeroU" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        // create wall distance property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWallDistance" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dof_dependencies", "L2" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Wall_Distance" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dof_derivative_functions", "Func_Wall_Distance_Der" );
        tPropCounter++;

        // create wall distance property
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInvWallDistanceSquare" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dof_dependencies", "PHID" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Wall_InvDistanceSquare" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dof_derivative_functions", "Func_Wall_InvDistanceSquare_Der" );
        tPropCounter++;

        // inlet viscosity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Inlet_V" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tInletKinViscosity );
        tPropCounter++;
        
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSelectInletV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Select_Inlet_V" );
        tPropCounter++;

        // zero viscosity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDirichletZeroV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0" );
        tPropCounter++;

        // upwind
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInletUpwind" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Inlet_Upwind" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // time continuity weights for velocity and viscosity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightUV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Time_Weights");
        tPropCounter++;
        
	    // time continuity weights for velocity and viscosity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightResUV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Time_Weight_Res");
        tPropCounter++;

        // initial condition velocity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialConditionU" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        // initial condition viscosity
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialConditionV" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tInletKinViscosity );
        tPropCounter++;
	
	// select integration domain
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropSelectDomain" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "Func_Select_Int_Domain");
        tPropCounter++;

        // Wall distance properties ----------------------------------------------------
        // common properties for theta and phi problems

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropConductivity" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tConductivity );
        tPropCounter++;

        // properties for Theta

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensityTheta" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tDensityTheta );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCapacityTheta" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tCapacityTheta );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPrescTheta" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPrescTheta );
        tPropCounter++;

        // properties for phi problem

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropDensityPhi" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tDensityPhi );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropCapacityPhi" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tCapacityPhi );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropPrescPhi" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tPrescPhi );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropEigenStrainPhi" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // time continuity weights
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightCurrent" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "10.0" );
        tPropCounter++;

        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropWeightPrevious" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "10.0" );
        tPropCounter++;

        // initial condition
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropInitialCondition" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "0.0" );
        tPropCounter++;

        // target level set
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLevelSetConst" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tPropCounter++;

        // target level set gradient
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLevelSetGradxConst" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", tLevelSetGradxConstant );
        tPropCounter++;

        // actual level set
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLevelSet" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "tLevelSetFunc" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // actual level set gradient
        tParameterList( tPropIndex ).push_back( prm::create_property_parameter_list() );
        tParameterList( tPropIndex )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( tPropIndex )( tPropCounter ).set( "property_name", "PropLevelSetGradx" );
        tParameterList( tPropIndex )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( tPropIndex )( tPropCounter ).set( "value_function", "tLevelSetGradxFunc" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        tParameterList( tPropIndex )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // init CM counter
        uint tCMCounter = 0;

        // incompressible NS fluid constitutive model
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::FLUID_TURBULENCE );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P;VISCOSITY", "Velocity,Pressure,Viscosity" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
            "PropFluidViscosity   ,Viscosity;"
            "PropFluidKinViscosity,KinViscosity;"
            "PropFluidDensity     ,Density" );
        tCMCounter++;
        
        // Spalart Allmaras turbulence constitutive model
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMTurbulence" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        tParameterList( tCMIndex )( tCMCounter ).set( "function_parameters", tCMTurbFt2 + "/" + tCMTurbAlpha );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
            "PropFluidKinViscosity,KinViscosity;"
            "PropWallDistance     ,WallDistance" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // create parameter list for constitutive model - Theta problem
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluidDiffusionTheta" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );
        tCMCounter++;

        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMInclusionSolidDiffusionTheta" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseInclusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );
        tCMCounter++;

        // create parameter list for constitutive model - Phi problem
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMFluidDiffusionPhi" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseFluid" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );
        tCMCounter++;

        // create parameter list for constitutive model - Phi problem
        tParameterList( tCMIndex ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_name", "CMInclusionSolidDiffusionPhi" );
        tParameterList( tCMIndex )( tCMCounter ).set( "phase_name", "PhaseInclusion" );
        tParameterList( tCMIndex )( tCMCounter ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        tParameterList( tCMIndex )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        tParameterList( tCMIndex )( tCMCounter ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // init SP counter
        uint tSPCounter = 0;

        // SP for fluid phase ----------------------------------------------------------
        // SUPG_PSPG NS incompressible
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPIncFlow" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tSupgFluidC1 + "/0.0" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        tSPCounter++;

        // Dirichlet Nitsche for fluid velocity
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPFluidDirichletNitscheU" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", sIfcNitscheFluid );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                "PropFluidViscosity,Viscosity;"
                "PropFluidDensity,Density" );
        tSPCounter++;

        // SUPG Spalart-Allmaras
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPSUPGSA" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", tSupgTurbPower + "/" + tSupgTurbSource );    
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;VISCOSITY", "Velocity,Viscosity" ) );
        tSPCounter++;
        
        // Dirichlet Nitsche for fluid turbulent viscosity
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheV" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", sIfcNitscheViscosity );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tSPCounter++;

        // SP for ghost ----------------------------------------------------------------
        if ( tUseGhost )
        {
            // ghost fluid viscous
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPViscous" );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::VISCOUS_GHOST );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_velo );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidViscosity,Viscosity" );
            tSPCounter++;

            // ghost fluid convective
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPVelocity" );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::CONVECTIVE_GHOST );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_conv );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidDensity,Density" );
            tSPCounter++;

            // ghost fluid pressure
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPPressure" );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::PRESSURE_GHOST );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_pres );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties",
                    "PropFluidViscosity,Viscosity;"
                    "PropFluidDensity,Density" );
            tSPCounter++;

            // ghost fluid turbulence viscosity
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPViscosity" );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_visc );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropFluidKinViscosity,Material" );
            tSPCounter++;
            
            // create parameter list for projection of the distance
            tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPFluidL2" );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_l2 );
            tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
            tSPCounter++;
         }

        // create parameter list for ghost stabilization parameter for theta and phi problems
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPFluidThetaPhi" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_heat );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for DBC on interface for theta problem
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheFluidThetaPhi" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", sIfcNitscheHeat );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for theta and phi problems
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPGPInclusionThetaPhi" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tSPIndex )( tSPCounter ).set( "follower_phase_name", "PhaseInclusion" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", alpha_heat );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for DBC on interface for theta problem
        tParameterList( tSPIndex ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_name", "SPNitscheInclusionThetaPhi" );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tSPIndex )( tSPCounter ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( tSPIndex )( tSPCounter ).set( "function_parameters", sIfcNitscheHeat );
        tParameterList( tSPIndex )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // init IWG counter
        uint tIWGCounter = 0;

        // fluid bulk IWG --------------------------------------------------------------
        // NS incompressible for velocity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidVelocityBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tIWGCounter++;

        // NS incompressible for pressure
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidPressureBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropFluidPressureSpring,PressureSpring" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        tIWGCounter++;

        // turbulence
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGTurbulenceBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPSUPGSA,SUPG" );
        tIWGCounter++;
        
        // L2 projection for distance
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDistance" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::L2 );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "L2" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInvWallDistanceSquare,Source" );
        tIWGCounter++;

        // fluid BC IWG ----------------------------------------------------------------
        // incompressible NS velocity Dirichlet IWG for inlet
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletVelocity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select;PropInletUpwind,Upwind" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );
        tIWGCounter++;

        // incompressible NS pressure Dirichlet IWG for inlet
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletU,Dirichlet;PropSelectInletU,Select" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tIWGCounter++;

        // zero velocity for velocity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidZeroVelocity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPFluidDirichletNitscheU,DirichletNitsche" );
        tIWGCounter++;

        // zero velocity for pressure
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidZeroPressure" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropDirichletZeroU,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        tIWGCounter++;
        
        if( tUsePseudoTimeStepping )
        {
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionU,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual");
            tIWGCounter++;
        }

        // inlet viscosity
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInletViscosity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "side_ordinals", "4" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropInletV,Dirichlet;PropSelectInletV,Select;PropInletUpwind,Upwind" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );
        tIWGCounter++;
        
        // zero viscosity for walls
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGZeroViscosity" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseVoid,PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VISCOSITY" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropDirichletZeroV,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheV,Nitsche" );
        tIWGCounter++;

        
        if( tUsePseudoTimeStepping )
        {
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VISCOSITY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                    "PropWeightUV,             WeightCurrent;"
                    "PropWeightUV,             WeightPrevious;"
                    "PropInitialConditionV,    InitialCondition;"
                    "PropWeightResUV,          WeightResidual" );
            tIWGCounter++;
        }

        if ( tUseGhost )
        {
            // ghost fluid viscous
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPViscous" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPViscous,GhostSP" );
            tIWGCounter++;

            // ghost fluid convective
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPConvective" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VX,VY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
            tIWGCounter++;

            // ghost fluid pressure
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPPressure" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "P" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPPressure,GhostSP" );
            tIWGCounter++;
            
            // ghost fluid viscosity
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPFluidViscosity" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "VISCOSITY" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPViscosity,GhostSP" );
            tIWGCounter++;
            
            // ghost L2 projection distance
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPInnerL2" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "L2" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPFluidL2,GhostSP" );
            tIWGCounter++;
        }
        
        //------------------------------------------------------------------------------
        // theta problem

        // theta bulk in fluid
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidDiffusionThetaBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        tIWGCounter++;

        // theta bulk in inclusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInclusionSolidDiffusionThetaBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionTheta,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusionTheta,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPInnerTheta" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );
            tIWGCounter++;

            // create IWG - ghost
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPOuterTheta" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseInclusion" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
            tIWGCounter++;
        }

        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGFluidTimeContinuityTheta" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tIWGCounter++;

        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGInclusionTimeContinuityTheta" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // phid problem

        // create IWG - bulk diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionInnerBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        tIWGCounter++;

        // create IWG - bulk diffusion
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseInclusion,PhaseVoid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMFluidDiffusionPhi,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheFluidThetaPhi,DirichletNitsche" );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "neighbor_phases", "PhaseFluid" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_constitutive_models", "CMInclusionSolidDiffusionPhi,Diffusion" );
        tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheInclusionThetaPhi,DirichletNitsche" );
        tIWGCounter++;

        if ( tUseGhost )
        {
            // create IWG - ghost
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPInnerPhi" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseFluid" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPFluidThetaPhi,GhostSP" );
            tIWGCounter++;

            // create IWG - ghost
            tParameterList( tIWGIndex ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_name", "IWGGPOuterPhi" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "leader_phase_name", "PhaseInclusion" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "follower_phase_name", "PhaseInclusion" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "dof_residual", "PHID" );
            tParameterList( tIWGIndex )( tIWGCounter ).set( "stabilization_parameters", "SPGPInclusionThetaPhi,GhostSP" );
            tIWGCounter++;
        }


        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVX" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VX,VY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tIQICounter++;

        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkP" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTHETA" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseHeatMethod" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "THETA" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkPHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseHeatMethod" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // fluid volume
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIFluidVolume" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tIQICounter++;

        // inclusion volume
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQISolidVolume" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseInclusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tIQICounter++;

        // inclusion perimeter
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIPerimeterItf" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases", "PhaseInclusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tIQICounter++;

        // fluid power dissipation in volume
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIVolumePowDisp" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::POWER_DISSIPATION_BULK );
 	tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropSelectDomain,Select" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid" );
        tIQICounter++;
        
         // level set 
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQILevelSet" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseHeatMethod" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropLevelSet,Property" );
        tIQICounter++;

        // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIHeatMethodPenaltyFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tIQICounter++;

        // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIHeatMethodPenaltyInclusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseInclusion" );
        tIQICounter++;

        // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIProjectedLevelSetFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 3 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tIQICounter++;

       // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIProjectedLevelSetInclusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 3 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseInclusion" );
        tIQICounter++;
	
       // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIHeatAlphaFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 5 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/1.0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tIQICounter++;

       // heat method penalty
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIHeatAlphaInclusion" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 5 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties",
                "PropLevelSet,L2_Reference;"
                "PropLevelSetGradx,H1S_Reference" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) 
                                                                               + "/" + tWPhi1 + "/" + tWPhi2 + "/" + tWGradPhi1 + "/" + tWGradPhi2 + "/-1.0");
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseInclusion" );
        tIQICounter++;

        // viscosity dof
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkVISCOSITY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "VISCOSITY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // L2 projection distance dof
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkL2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity", "L2" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tIQICounter++;

        // wall distance
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIWallDistance" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_properties", "PropWallDistance,Property" );
        tIQICounter++;

        // turbulent dynamic viscosity
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkTurbDynVisc" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::TURBULENT_DYNAMIC_VISCOSITY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );
        tIQICounter++;
    
        // effective dynamic viscosity
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkEffDynVisc" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,Fluid_Turbulence" );
        tIQICounter++;
/*
        // SUPG fluid
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkSUPGU" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPIncFlow,Stabilization" );
        tIQICounter++; 
        
        // PSPG fluid
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkPSPGU" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      1 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPIncFlow,Stabilization" );
        tIQICounter++; 
        
        // SUPG viscosity
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkSUPGSA" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STABILIZATION );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "stabilization_parameters",   "SPSUPGSA,Stabilization" );
        tIQICounter++;
	
        // incompressible NS strong residual
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UX" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      1 );
        tIQICounter++; 
	
        // incompressible NS strong residual
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_UY" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      2 );
        tIQICounter++; 
	
        // incompressible NS strong residual
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkIncompressibleNaveierStokesSF_P" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMFluid,IncompressibleFluid") ;
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      0 );
        tIQICounter++; 

        // SA strong residual
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIBulkSpalartAllmarasSF" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   ( uint ) fem::IQI_Type::STRONG_RESIDUAL_SA );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence") ;
        tIQICounter++; 
        // production coefficient
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkProductionCoeff" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tIQICounter++;
        
        // wall destruction coefficient
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkWallDestructionCoeff" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tIQICounter++;
        
        // diffusion coefficient
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkDiffusionCoeff" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 2 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tIQICounter++;
        
         // production term
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkProductionTerm" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 3 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tIQICounter++;
        
        // wall destruction term
        tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name", "IQIBulkWallDestructionTerm" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name", "PhaseFluid" );
        tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index", 4 );
        tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type", fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT );
        tParameterList( tIQIIndex )( tIQICounter ).set( "leader_constitutive_models", "CMTurbulence,SpalartAllmarasTurbulence" );
        tIQICounter++;
*/
        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( tFEMIndex ).push_back( prm::create_computation_parameter_list() );
        tParameterList( tFEMIndex )( 0 ).set( "print_physics_model", false );
        tParameterList( tFEMIndex )( 0 ).set( "is_analytical_sensitivity", false );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        tParameterList( tFEMIndex )( 0 ).set( "finite_difference_perturbation_size", stod( tFDEpsilon ) );
    }

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        //------------------------------------------------------------------------------

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );

        //------------------------------------------------------------------------------

        tParameterlist( 2 ).resize( 6 );

        // Newton solve in forward problem used by Newton flow subproblem (includes load control)
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );
        tParameterlist( 2 )( 0 ).set( "NLA_ref_iter", tNLA_ref_its );        

        if ( tSolverType == "Newton" && tUseNewtonLoadControl )
        {
            tParameterlist( 2 )( 0 ).set( "NLA_load_control_strategy",  sol::SolverLoadControlType::Exponential ) ;
            tParameterlist( 2 )( 0 ).set( "NLA_load_control_factor", tRampInitial );
            tParameterlist( 2 )( 0 ).set( "NLA_load_control_steps",  tRampSteps   );
            tParameterlist( 2 )( 0 ).set( "NLA_load_control_relres", tRampRelRes  );
            tParameterlist( 2 )( 0 ).set( "NLA_load_control_exponent", 1.0 );
        }

        // for linear solve - forward: theta, phi, l2 - adjoint: all problems
        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 1 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 1 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 1 );

        // NLBGS for overall forward solve and adjoint in theta,phi,l2,(vx,vy,p,viscosity)
        tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    
        tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 2 )( 2 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0 );
        tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 1 );

        // algorithm not used 
        tParameterlist( 2 )( 3 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 3 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 3 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 3 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 3 ).set( "NLA_max_iter", 1 );
        
        // NLBGS for flow subproblem using pseudo time stepping and (optional) load control)
        tParameterlist( 2 )( 4 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    
        tParameterlist( 2 )( 4 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 2 )( 4 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 4 ).set( "NLA_rel_res_norm_drop",    tNLBGS_rel_res );
        tParameterlist( 2 )( 4 ).set( "NLA_relaxation_parameter", tNLBGS_realx );
        tParameterlist( 2 )( 4 ).set( "NLA_max_iter",             tNLBGS_max_itr );
        tParameterlist( 2 )( 4 ).set( "NLA_ref_iter",             tNLBGS_ref_its );        

        if ( tSolverType == "NLBGS" )
        {
            if ( tUsePseudoTimeStepping )
            {
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_control_strategy",     sol::SolverPseudoTimeControlType::Expur ) ;
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_initial_steps",       tNLBGS_init_itr );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_initial",             tSteadyStateStepSize );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_max_num_steps",       tMaxNumTimeSteps	);
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_max_step_size",       tMaxTimeStepSize );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_rel_res_norm_drop",   tRelResNormDrop );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_rel_res_norm_update", tRelResNormUpdate  );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_steady_rel_res_norm", tSteadyStateRelRes  );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_steady_step_size",    tSteadyStateStepSize  );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_offset",              tTimeOffSet );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_constant",            tConstantTimeStep );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_step_index_factor",   tIndexFactor );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_step_index_exponent", tIndexExpo );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_residual_factor",     tResFactor );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_residual_exponent",   tResExpo );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_comsol_1",            tComsolParameter1 );
                tParameterlist( 2 )( 4 ).set( "NLA_pseudo_time_comsol_2",            tComsolParameter2 );
            }
            if ( tUseNLBGSLoadControl )
            {
                tParameterlist( 2 )( 4 ).set( "NLA_load_control_strategy",  sol::SolverLoadControlType::Exponential ) ;
                tParameterlist( 2 )( 4 ).set( "NLA_load_control_factor",   tRampInitial );
                tParameterlist( 2 )( 4 ).set( "NLA_load_control_steps",    tRampSteps   );
                tParameterlist( 2 )( 4 ).set( "NLA_load_control_relres",   tRampRelRes  );
                tParameterlist( 2 )( 4 ).set( "NLA_load_control_exponent", 1.0 );
            }
        }

        // Newton solve used within NLBGS for flow subproblem
        tParameterlist( 2 )( 5 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    
        tParameterlist( 2 )( 5 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 2 )( 5 ).set( "NLA_Linear_solver", 0 );
        tParameterlist( 2 )( 5 ).set( "NLA_rel_res_norm_drop",    tNewton_rel_res );
        tParameterlist( 2 )( 5 ).set( "NLA_max_iter",             tNewton_max_iter );               
        tParameterlist( 2 )( 5 ).set( "NLA_relaxation_strategy",   sol::SolverRelaxationType::InvResNormAdaptive ) ;
        tParameterlist( 2 )( 5 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 5 ).set( "NLA_relaxation_damping",   0.5 );
 
        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 9 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();   
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "1" );            
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "THETA" );

        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();    
        tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             
        tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "PHID" );

        if ( tSolverType == "NLBGS" )
        {
            if ( tFluidViscSolver == "Mono" )  // monolithic flow problem
            {
                tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();    
                tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );           
                tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
            }
            else // flow problem in vx,vy,p
            {
                tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();    
                tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );             
                tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
                tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "VX,VY,P" );
            }
        }
        else  // using Newton solver only monolithic 
        {
            tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();   
            tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             
            tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
            tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }

        // flow viscosity problem only
        tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();   
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", tSubNewtonSolverOption );           
        tParameterlist( 3 )( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "VISCOSITY" );

        // following settings only used when solving flow problem by NLBGS
        if ( tFluidViscSolver == "Mono" )
        {
            tParameterlist( 3 )( 4 ) = moris::prm::create_nonlinear_solver_parameter_list();    
            tParameterlist( 3 )( 4 ).set( "NLA_Nonlinear_solver_algorithms", "4" );             
            tParameterlist( 3 )( 4 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            tParameterlist( 3 )( 4 ).set( "NLA_Sub_Nonlinear_Solver", "2" );    
            tParameterlist( 3 )( 4 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );
        }
        else
        {
            tParameterlist( 3 )( 4 ) = moris::prm::create_nonlinear_solver_parameter_list();    
            tParameterlist( 3 )( 4 ).set( "NLA_Nonlinear_solver_algorithms", "4" );            
            tParameterlist( 3 )( 4 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            tParameterlist( 3 )( 4 ).set( "NLA_Sub_Nonlinear_Solver", "2,3" );   
            tParameterlist( 3 )( 4 ).set( "NLA_DofTypes", "VX,VY,P;VISCOSITY" );
        }

        // L2 projection problem
        tParameterlist( 3 )( 5 ) = moris::prm::create_nonlinear_solver_parameter_list();    
        tParameterlist( 3 )( 5 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             
        tParameterlist( 3 )( 5 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 5 ).set( "NLA_DofTypes", "L2" );

        // total forward problem
        if ( tSolverType == "NLBGS" )
        {
            tParameterlist( 3 )( 6 ) = moris::prm::create_nonlinear_solver_parameter_list();    
            tParameterlist( 3 )( 6 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             
            tParameterlist( 3 )( 6 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            tParameterlist( 3 )( 6 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,4" );   
            tParameterlist( 3 )( 6 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }
        else
        {
            tParameterlist( 3 )( 6 ) = moris::prm::create_nonlinear_solver_parameter_list();    
            tParameterlist( 3 )( 6 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             
            tParameterlist( 3 )( 6 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
            tParameterlist( 3 )( 6 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,2" );    
            tParameterlist( 3 )( 6 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        }

        // adjoint of flow problem
        tParameterlist( 3 )( 7 ) = moris::prm::create_nonlinear_solver_parameter_list();    
        tParameterlist( 3 )( 7 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             
        tParameterlist( 3 )( 7 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        tParameterlist( 3 )( 7 ).set( "NLA_DofTypes", "VX,VY,P,VISCOSITY" );

        // adjoint of total problem
        tParameterlist( 3 )( 8 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 8 ).set( "NLA_Nonlinear_solver_algorithms", "2" );
        tParameterlist( 3 )( 8 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        tParameterlist( 3 )( 8 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,5,7" );    // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 8 ).set( "NLA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );

        // ----------------------------------------------------------

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 6 );                      // using NLBGS for forward problem
        tParameterlist( 4 )( 0 ).set( "TSA_nonlinear_solver_for_adjoint_solve", 8 );    // using NLBGS for sensitivity problem

        //------------------------------------------------------------------------------

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "THETA;PHID;L2;VX,VY,P,VISCOSITY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "THETA,0.0;PHID,0.0;VX,"+tVXInitial+";VY,0.0;P,0.0;VISCOSITY," + tInletKinViscosity + ";L2,1.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //------------------------------------------------------------------------------

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        
        //------------------------------------------------------------------------------

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::NONE);
    }

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).push_back( prm::create_msi_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "THETA", 0 );
        tParameterlist( 0 )( 0 ).set( "PHID", 0 );
        tParameterlist( 0 )( 0 ).set( "VX", 0 );
        tParameterlist( 0 )( 0 ).set( "VY", 0 );
        tParameterlist( 0 )( 0 ).set( "P", 0 );
        tParameterlist( 0 )( 0 ).set( "VISCOSITY", 0 );
        tParameterlist( 0 )( 0 ).set( "L2", 0 );
    }

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );

        tParameterlist( 0 ).push_back( prm::create_vis_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tExoFile ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 )( 0 ).set( "Set_Names", tAllDomain );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "THETA,PHID,L2,"
                                                    "VX,VY,P,VISCOSITY,"
                                                    "TURBDYNVISC,EFFVISC,"
                                                    "LS,WALLDIST,PRJLEVSFLD,PRJLEVSINC,ALPHAFLD,ALPHAINC");
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,"
                                                    "NODAL,NODAL,NODAL,NODAL,"
                                                    "NODAL,NODAL,"
                                                    "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkTHETA,IQIBulkPHID,IQIBulkL2,"
                                                "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkVISCOSITY,"
                                                "IQIBulkTurbDynVisc,IQIBulkEffDynVisc,"
                                                "IQILevelSet,IQIWallDistance,IQIProjectedLevelSetFluid,IQIProjectedLevelSetInclusion,IQIHeatAlphaFluid,IQIHeatAlphaInclusion" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
