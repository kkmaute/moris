/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Zienkiewicz_Zhu_Adaptive_Refinement.cpp
 *
 */

#include <string>
#include <iostream>
#include <sstream>
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
#include "cl_MTK_Field.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
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

    std::string tName = "Zienkiewicz_Zhu_Adaptive_Refinement";

    std::string tLibraryName   = "Zienkiewicz_Zhu_Adaptive_Refinement.so";
    std::string tGENOutputFile = "GEN_" + tName + ".exo";

    std::string tOutputFileName = tName + ".exo";

    bool tIs3D             = false;
    bool tIsOpt            = true;
    bool tUseGhost         = true;
    bool tUseAbsoluteValue = true;

    bool tUseDensityShift = true;

    static real tDensityShift          = 0.0;
    uint        tDensityShiftStart     = 20;
    uint        tDensityShiftIntervall = 150;

    static uint tItarationCounter = 0;

    //-------------------------------
    // Opt constant_parameters

    real tMMAPenalty  = 100.0;
    real tMMAStepSize = 0.02;
    int  tMMAMaxIter  = 19;

    real tInitialStrainEnergy = 90.00;
    real tStrainEnergyPen     = 0.1;
    real tMaxMass             = 4;

    real tPerimeterPenalty = 0.1;
    real tInitialPerimeter = 8.0;

    real tRegularizationPenalty = 0.2;    // 0.02
    real tInitialRegularization = 250;

    //-------------------------------
    // geometry details

    moris::real tFirstHoleX = -2.0;
    moris::real tFirstHoleY = 0.00;

    moris::real tSecondHoleX = 2.0;
    moris::real tSecondHoleY = 0.0;

    moris::real tRadius = 1.237;

    moris::real tRadiusHolesInner = 0.517;
    moris::real tRadiusHolesOuter = 0.7173;

    //-------------------------------

    real tElementEdgeLength = 0.80;

    real phi_sh = 1.20;

    real phi_sh_scale = 1.2 / ( 3.0 * tElementEdgeLength );
    real phi_rt       = 3.0;

    int tDispOrder = 1;

    real tBsplineLimitTop    = 2.4;
    real tBsplineLimitBottom = 0.0;

    real tInitialDensity = 0.0001;

    moris::real tBSplineLimit = 1.2;
    moris::real tPhiBandwidth = 3.0 * tElementEdgeLength;
    moris::real tPhiGradient  = tBSplineLimit * std::log( 199.0 ) / ( 2.0 * tPhiBandwidth );
    moris::real tPhiGamma     = 2.0 * std::log( 10.0 ) / std::pow( 2.0 / ( std::exp( -2.0 * tPhiBandwidth * tPhiGradient / tBSplineLimit ) + 1.0 ) - 1.0, 2.0 );

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

    real tVMMaxDofNValue = 2.0;    // exponent in max_dof formulation
    real tVMShift        = 1.0;

    real tStressGhost = 0.001;

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
    //-------------------------------
    std::string tBulkSets  = "HMR_dummy_n_p9,HMR_dummy_c_p9";
    std::string tVoidSets  = "HMR_dummy_n_p8,HMR_dummy_c_p8";
    std::string tFrameSets = "HMR_dummy_n_p11,HMR_dummy_c_p11";

    std::string tDirichletSets = "iside_b0_11_b1_15";
    std::string tLoadSets      = "iside_b0_11_b1_15";

    std::string tInterfaceVoidSets = "iside_b0_9_b1_8";

    std::string tVoidInterfaceSets = "iside_b0_8_b1_9";

    std::string tFrameInteriorDSets = "dbl_iside_p0_9_p1_11";

    std::string tInteriorGhost = "ghost_p9";
    std::string tVoidGhost     = "ghost_p8";
    std::string tFrameGhost    = "ghost_p11";

    std::string tTotalDomain = tBulkSets + "," + tVoidSets;

    std::string tTotalDomainAGhost = tTotalDomain + "," + tInteriorGhost;

    //------------------------------------------------------------------------------
    enum hmr::ElementalRefienmentIndicator
    refinement_function(
            mtk::Cell*                    aElement,
            std::shared_ptr< mtk::Field > aField,
            uint                          tActivationPattern,
            uint&                         aMaxLevel )
    {
        // current refinement level of element
        uint tLevel = aElement->get_level();

        // moris_index tElementIndex = aElement->get_index();

        uint tMaxLevel = 2;

        //         if( tActivationPattern == 0)
        //         {
        //             tMaxLevel = 2;
        //         }
        //         else if( tActivationPattern == 1)
        //         {
        //             tMaxLevel = 1;
        //         }

        real lsth    = 1.2;
        real lsbwabs = 0.2;

        Matrix< IndexMat > tVertexInds = aElement->get_vertex_inds();

        Matrix< DDRMat >   tValues;
        Matrix< IndexMat > tFieldIndex( 1, 1, 0 );

        aField->get_value(
                tVertexInds,
                tValues,
                tFieldIndex );

        enum hmr::ElementalRefienmentIndicator tRefine = hmr::ElementalRefienmentIndicator::DROP;

        // refinement strategy
        if ( tValues.max() >= lsth - lsbwabs )
        {
            // for volume refinement
            if ( tValues.min() >= lsth + lsbwabs )
            {
            }
            // for interface refinement
            else
            {
                if ( tLevel < tMaxLevel )
                {
                    tRefine = hmr::ElementalRefienmentIndicator::REFINE;
                }
                else
                {
                    tRefine = hmr::ElementalRefienmentIndicator::HOLD;
                }
            }
        }
        else
        {
            //             if( curlevel <  minlevel )
            //             {
            //                 tRefine = 1; // refine
            //             }
            //             else if ( curlevel == minlevel )
            //             {
            //                 tRefine = 0; // keep
            //             }
        }

        return tRefine;
    }

    //------------------------------------------------------------------------------
    enum hmr::ElementalRefienmentIndicator
    refinement_function_stress(
            mtk::Cell*                    aElement,
            std::shared_ptr< mtk::Field > aField,
            uint                          tActivationPattern,
            uint&                         aMaxLevel )
    {
        // current refinement level of element
        uint tLevel = aElement->get_level();

        moris_index tElementIndex = aElement->get_index();

        uint tMaxLevel = 2;

        Matrix< IndexMat > tVertexInds( 1, 1, tElementIndex );

        Matrix< DDRMat >   tValues;
        Matrix< IndexMat > tFieldIndex( 1, 1, 0 );

        aField->get_value(
                tVertexInds,
                tValues,
                tFieldIndex );

        // real tMaxStressValue = aField->get_max_value();
        real tMaxStressValue = 0.1;

        enum hmr::ElementalRefienmentIndicator tRefine = hmr::ElementalRefienmentIndicator::DROP;

        if ( tValues( 0 ) != MORIS_REAL_MIN )
        {
            // refinement strategy
            if ( tValues( 0 ) / tMaxStressValue >= 0.7 )
            {
                if ( tLevel < tMaxLevel )
                {
                    tRefine = hmr::ElementalRefienmentIndicator::REFINE;
                }
            }
            else if ( tValues( 0 ) / tMaxStressValue >= 0.2 )
            {
                tRefine = hmr::ElementalRefienmentIndicator::HOLD;
            }
            else
            {
                tRefine = hmr::ElementalRefienmentIndicator::COARSEN;
            }
        }

        return tRefine;
    }

    //------------------------------------------------------------------------------
    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.0 && aFIManager->get_IG_geometry_interpolator()->valx()( 0 ) > 0.0 )
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
    void
    Func_Select(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );

        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 0 ) < 0.0 )
        {
            if ( tIs3D )
            {
                // aPropMatrix = {{ 0.0 },{ aParameters( 0 )( 0 ) },{ 0.0 }};
            }
            else
            {
                aPropMatrix( 0, 0 ) = 1.0;
                aPropMatrix( 1, 1 ) = 1.0;
            }
        }
    }
    //------------------------------------------------------------------------------
    void
    PhiD_Prop(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tPhid = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );

        real tVal = std::exp( -2.0 * tPhiGradient * tPhid / tBSplineLimit );

        aPropMatrix.set_size( 1, 1, 0.0 );
        aPropMatrix( 0 ) = ( 2.0 / ( 1.0 + tVal ) - 1.0 ) * tBSplineLimit;
    }
    //------------------------------------------------------------------------------

    real
    Const_Geometry(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )
    {
        return 2.0;
    }

    //------------------------------------------------------------------------------
    // Hole pattern
    real
    general_shape_outer(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )

    {
        real tLSval;

        if ( aCoordinates( 0 ) <= tFirstHoleX )
        {
            tLSval = std::pow( std::pow( ( aCoordinates( 0 ) - tFirstHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tFirstHoleY ), 2 ), 0.5 ) - tRadius;
        }
        else if ( aCoordinates( 0 ) >= tSecondHoleX )
        {
            tLSval = std::pow( std::pow( ( aCoordinates( 0 ) - tSecondHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tSecondHoleY ), 2 ), 0.5 ) - tRadius;
        }
        else if ( aCoordinates( 1 ) <= 0.0 )
        {
            tLSval = -1.0 * ( aCoordinates( 1 ) + tRadius );
        }
        else if ( aCoordinates( 1 ) >= 0.0 )
        {
            tLSval = aCoordinates( 1 ) - tRadius;
        }
        else
        {
            MORIS_ERROR( false, "should never be in here" );
        }

        // clean return value to return non-zero value
        return std::abs( tLSval ) < 1e-8 ? ( 1e-8 + tBSplineLimit ) : ( -1.0 * tLSval + tBSplineLimit );
    }
    //------------------------------------------------------------------------------
    // Hole pattern
    real
    Hole_inner(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )
    {
        real tLSval1 = std::pow( std::pow( ( aCoordinates( 0 ) - tFirstHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tFirstHoleY ), 2 ), 0.5 ) - tRadiusHolesInner;
        real tLSval2 = std::pow( std::pow( ( aCoordinates( 0 ) - tSecondHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tSecondHoleY ), 2 ), 0.5 ) - tRadiusHolesInner;

        real tLSval = std::min( tLSval1, tLSval2 );

        // clean return value to return non-zero value
        return std::abs( tLSval ) < 1e-8 ? ( 1e-8 + tBSplineLimit ) : ( -1.0 * tLSval + tBSplineLimit );
    }
    //------------------------------------------------------------------------------
    real
    Hole_outer(
            const Matrix< DDRMat >& aCoordinates,
            const Cell< real* >&    aGeometryParameters )
    {
        real tLSval1 = std::pow( std::pow( ( aCoordinates( 0 ) - tFirstHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tFirstHoleY ), 2 ), 0.5 ) - tRadiusHolesOuter;
        real tLSval2 = std::pow( std::pow( ( aCoordinates( 0 ) - tSecondHoleX ), 2 ) + std::pow( ( aCoordinates( 1 ) - tSecondHoleY ), 2 ), 0.5 ) - tRadiusHolesOuter;

        real tLSval = std::min( tLSval1, tLSval2 );

        // clean return value to return non-zero value
        return std::abs( tLSval ) < 1e-8 ? ( 1e-8 + tBSplineLimit ) : ( -1.0 * tLSval + tBSplineLimit );
    }

    //------------------------------------------------------------------------------
    void
    compute_density_shift()
    {
        if ( tUseDensityShift )
        {
            if ( tItarationCounter >= tDensityShiftStart && tItarationCounter < ( tDensityShiftStart + tDensityShiftIntervall ) )
            {
                real tPow     = std::pow( ( ( (real)tItarationCounter - (real)tDensityShiftStart ) / ( (real)tDensityShiftIntervall ) ), 2.0 );
                tDensityShift = tInitialDensity + ( 1 - tInitialDensity ) * tPow;
            }
            else if ( tItarationCounter < tDensityShiftStart )
            {
                tDensityShift = 0.0;
            }
            else
            {
                tDensityShift = 1.0;
            }
        }
    }

    //------------------------------------------------------------------------------
    void
    tYoungsFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        moris::Matrix< moris::DDRMat > tHCT  = aParameters( 0 );
        moris::real                    tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 ) / 2.40;

        real tDensity = ( tLevelSet - phi_sh_scale ) / ( 1 - phi_sh_scale );

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix = tHCT * std::pow( tDensity, tBeta );
    }

    //------------------------------------------------------------------------------
    void
    tDerYoungsFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        moris::Matrix< moris::DDRMat > tHCT  = aParameters( 0 );
        moris::real                    tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        // FIXME density shift missing

        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->N() *    //
                      tBeta * tHCT( 0 ) * std::pow( tDensity, tBeta - 1 ) / ( 1 - phi_sh );
    }

    //------------------------------------------------------------------------------

    void
    tDensityFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tLevelSet = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 ) / 2.40;

        real tDensity = ( tLevelSet - phi_sh_scale ) / ( 1 - phi_sh_scale );

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix.set_size( 1, 1, tDensity );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->N() / ( 1 - phi_sh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tLevelSetFuncReal(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val();
    }

    /* ------------------------------------------------------------------------ */

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // Matrix< DDRMat > tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

        // bool tBool = is_in_vicinity_of_load( tCoords );

        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
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

        value = value - phi_sh;

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
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

        value = value - phi_sh;

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
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

        value = value - phi_sh;

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( PDV_Type::LS1 )->dnNdxn( 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------
    /*
Matrix<DDRMat> compute_objectives(Matrix<DDRMat> aADVs, Matrix<DDRMat> aCriteria)
{
    Matrix<DDRMat> tObjectives( 1, 1 );

    real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
    real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
    real obj3 = tPerimeterPenalty * aCriteria( 3 ) / tInitialPerimeter;
    real obj4 = tRegularizationPenalty* aCriteria( 4 ) / tInitialRegularization;
    //real obj5 = tRegularizationPenalty* aCriteria( 5 ) ;

    tObjectives( 0, 0 ) = obj1 + obj2 + obj3 + obj4;// + obj2 + obj3 + obj4 + obj5;

    std::cout << "% --------------------------------- % \n";
    std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
    std::cout << "Strain Energy            = " << aCriteria( 0 )    << " ( " <<  obj1 / tObjectives( 0, 0 ) << " )\n";
    std::cout << "Strain Energy (Interior) = " << aCriteria( 1 )      << " ( " <<  obj2 / tObjectives( 0, 0 ) << " )\n";
    std::cout << "Perimeter                = " << aCriteria( 3 )      << " ( " <<  obj3 / tObjectives( 0, 0 ) << " )\n";
    std::cout << "Regularization            = " << aCriteria( 4 )      << " ( " <<  obj4 / tObjectives( 0, 0 ) << " )\n";
    //std::cout << "H1Error PDV              = " << aCriteria( 5 )      << " ( " <<  aCriteria( 1 ) / tObjectives( 0, 0 ) << " )\n";
    std::cout << " \n";

    std::cout << "min ADV                  = " << aADVs.min()         << " \n";
    std::cout << "max ADV                  = " << aADVs.max()         << " \n" << std::flush;

    return tObjectives;
}
     */

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( Matrix< DDRMat > aADVs, Matrix< DDRMat > aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real tStrainEnergyInner = aCriteria( 0 );
        real tStrainEnergyOuter = aCriteria( 1 );
        real tMass              = aCriteria( 2 );
        real tPerimeter         = aCriteria( 3 );
        real tHeatMethodPen     = aCriteria( 4 );
        // real tStress            = aCriteria( 5 );

        real obj1  = tStrainEnergyPen * tStrainEnergyInner / tInitialStrainEnergy;
        real obj11 = tStrainEnergyPen * tStrainEnergyOuter / tInitialStrainEnergy;
        real obj2  = tMass;
        real obj3  = tPerimeterPenalty * tPerimeter / tInitialPerimeter;
        real obj4  = tRegularizationPenalty * tHeatMethodPen / tInitialRegularization;
        // real obj5 = tRegularizationPenalty* tStress ;

        tObjectives( 0, 0 ) = obj1 + obj11 + obj2 + obj3 + obj4;    // + obj2 + obj3 + obj4 + obj5;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy            = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj11 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Mass                     = " << aCriteria( 2 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Regularization            = " << aCriteria( 4 ) << " ( " << obj4 / tObjectives( 0, 0 ) << " )\n";
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
        tConstraints( 0 ) = aCriteria( 5 );

        std::cout << "Volume     = " << aCriteria( 1 ) << " \n";
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

        tDObjectiveDCriteria( 0 ) = tStrainEnergyPen * 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 1 ) = tStrainEnergyPen * 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 2 ) = 1.0;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenalty / tInitialPerimeter;
        tDObjectiveDCriteria( 4 ) = tRegularizationPenalty / tInitialRegularization;
        // tDObjectiveDCriteria( 4 ) = tRegularizationPenalty / tInitialRegularization;
        // tDObjectiveDCriteria( 5 ) = tRegularizationPenalty* 1.0;

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

        tDConstraintDCriteria( 5 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        compute_density_shift();

        tItarationCounter++;

        return true;
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
        tParameterlist( 0 )( 0 ).set( "restart_file", "" );
        tParameterlist( 0 )( 0 ).set( "reinitialize_interface_iter", 18 );

        tParameterlist( 2 )( 0 ) = moris::prm::create_gcmma_parameter_list();
        tParameterlist( 2 )( 0 ).set( "step_size", tMMAStepSize );
        tParameterlist( 2 )( 0 ).set( "penalty", tMMAPenalty );
        tParameterlist( 2 )( 0 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        tParameterlist( 2 )( 0 ).set( "restart_index", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "60,30" );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", "8,4" );
        tParameterlist( 0 )( 0 ).set( "domain_offset", "-4.0,-2.0" );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", "2" );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", "1,2" );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "1,2" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0,0" );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        // tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        // tParameterlist( 0 )( 0 ).set( "use_refine_low_level_elements", true );
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
        tParameterlist( 0 )( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy,IQIBulkVolume,IQIPerimeter_InterfaceVoid,IQIHeatMethodPenalty,IQIMaxStress" );
        // tParameterlist(0)(0).set("output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );
        // FIXME     this has to change
        // tParameterlist( 0 )( 0 ).set("PDV_types"         , "LS1");

        Matrix< DDUMat > tPhaseMap( 16, 1, 0 );
        tPhaseMap( 1 )  = 1;
        tPhaseMap( 2 )  = 2;
        tPhaseMap( 3 )  = 3;
        tPhaseMap( 4 )  = 4;
        tPhaseMap( 5 )  = 5;
        tPhaseMap( 6 )  = 6;
        tPhaseMap( 7 )  = 7;
        tPhaseMap( 8 )  = 8;
        tPhaseMap( 9 )  = 9;
        tPhaseMap( 10 ) = 11;
        tPhaseMap( 11 ) = 11;
        tPhaseMap( 12 ) = 12;
        tPhaseMap( 13 ) = 13;
        tPhaseMap( 14 ) = 14;
        tPhaseMap( 15 ) = 15;
        tParameterlist( 0 )( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );
        tParameterlist( 0 )( 0 ).set( "print_phase_table", true );

        // init geometry counter
        uint tGeoCounter = 0;

        // outer frame
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "general_shape_outer" );
        // tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "0,0" );
        // tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 10e-14 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 1.2 );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Outer_shape" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Hole_inner" );
        // tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "0,0" );
        // tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 10e-14 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 1.2 );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Inner_Holes" );
        tGeoCounter++;

        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Hole_outer" );
        // tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "0,0" );
        // tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 10e-14 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 1.2 );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Outer_Holes" );
        tGeoCounter++;

        // initialize fins as swiss cheese geometry
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Const_Geometry" );
        tParameterlist( 1 )( tGeoCounter ).set( "name", "Level_Set_Field" );
        // tParameterlist( 1 )( tGeoCounter ).set( "number_of_refinements", "0,0" );
        // tParameterlist( 1 )( tGeoCounter ).set( "refinement_mesh_index", "0,1" );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_tolerance", 10e-14 );
        tParameterlist( 1 )( tGeoCounter ).set( "isocontour_threshold", 1.2 );
        tParameterlist( 1 )( tGeoCounter ).set( "multilinear_intersections", false );

        if ( tIsOpt )
        {
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 1 );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", tBsplineLimitBottom );
            tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", tBsplineLimitTop );
        }
        tGeoCounter++;

        uint tParamCounter = 0;
        tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list() );
        tParameterlist( 2 )( tParamCounter ).set( "name", "LvL_Set_Field" );
        tParameterlist( 2 )( tParamCounter ).set( "dependencies", "Level_Set_Field" );
        tParameterlist( 2 )( tParamCounter ).set( "type", "scaled_field" );
        tParameterlist( 2 )( tParamCounter ).set( "constant_parameters", "1.0" );
        tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
        // tParameterlist( 2 )( tParamCounter ).set("discretization_mesh_index",   -1);
        // tParameterlist( 2 )( tParamCounter ).set("discretization_lower_bound", 0.001);
        // tParameterlist( 2 )( tParamCounter ).set("discretization_upper_bound", 1.0);

        if ( tUseGhost )
        {
            tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomainAGhost );
        }
        else
        {
            tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomain );
        }

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
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tDensityFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerDensityFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensityFrame" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0,3.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tYoungsFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerYoungsFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungsFrame" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0e-8" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelect" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Select" );
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

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropUnitValue" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
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
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0;1.0" );
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

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropTruncPHID" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "PhiD_Prop" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIsoFrame" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungsFrame,YoungsModulus;PropPoisson,PoissonRatio" );
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

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungsFrame,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungs,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhostFrame" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.001" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungsFrame,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGhost_VM" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::to_string( tStressGhost ) + "/0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungsFrame,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropUnitValue,Material" );
        tSPCounter++;

        // create parameter list for ghost stabilization parameter for theta and phi problems
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "0.001" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", "PropConductivity,Material" );
        tSPCounter++;

        // create parameter list for DBC on interface for theta problem
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropConductivity,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // create IWG  - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionThetaBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVoidInterfaceSets );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityThetaInterior" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID;UX,UY;STRESS_DOF" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulkSets );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTimeContinuityThetaVoid" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVoidSets );
        tParameterList( 3 )( tIWGCounter ).set( "time_continuity", true );
        tIWGCounter++;

        // create IWG - bulk diffusion
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tTotalDomain );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tVoidInterfaceSets );
        tIWGCounter++;

        // create parameter list for single side interface condition
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInterfaceVoidSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulkSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkUFrame" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet;PropSelect,Select" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tDirichletSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLoadSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dv_dependencies", "LS1" );
        tParameterList( 3 )( tIWGCounter ).set( "follower_dv_dependencies", "LS1" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMStrucLinIso,ElastLinIso" ) );
        tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMStrucLinIsoFrame,ElastLinIso" ) );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameInteriorDSets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGStress" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_VON_MISES_STRESS ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "STRESS_DOF" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "STRESS_DOF" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulkSets );
        tIWGCounter++;

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhost" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dv_dependencies", "LS1" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dv_dependencies", "LS1" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostFrame" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "UX,UY" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhostFrame,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tFrameGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGSStress" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "STRESS_DOF" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "STRESS_DOF" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "STRESS_DOF" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_VM,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhost" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "PHID" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhost" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", "THETA" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets + "," + tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets + "," + tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStress" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tFrameSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tInterfaceVoidSets );
        tIQICounter++;

        // smooth fuction
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQILevelSet" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSetReal,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQITruncPHID" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropTruncPHID,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQILevelSetHeatMethod" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkTHETA" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // Nodal PHID IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkPHID" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // H1 Error if reference is constant
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1ErrorConst" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // H1 Error if reference is design dependent
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIH1Error" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::H1_ERROR ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        // H1 Error if reference is design dependent
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIHeatMethodPenalty" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::HEAT_METHOD_PENALTY ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "PHID" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "THETA;PHID" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + " / 0.0/ 0.2 / 1.0 / 1.0 " );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomain );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIMaxStress" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "function_parameters", std::to_string( 5.0 ) + "/" + std::to_string( tVMMaxDofNValue ) + "/" + std::to_string( tVMShift ) + "/1" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkZienkiewiczZhu" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::ZIENKIEWICZ_ZHU_VON_MISES_STRESS ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulkSets );
        tIQICounter++;

        /*
    tParameterList( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
    tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_name",                   "IQIMaxStress2");
    tParameterList( tIQIIndex )( tIQICounter ).set( "leader_phase_name",          "PhaseLBracket" );
    tParameterList( tIQIIndex )( tIQICounter ).set( "neighbor_phases",            "PhaseVoid2" );
    tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_bulk_type",              static_cast< uint >( fem::Element_Type::SIDESET ) );
    tParameterList( tIQIIndex )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::MAX_DOF ) );
    tParameterList( tIQIIndex )( tIQICounter ).set( "dof_quantity",               "STRESS_DOF");
    tParameterList( tIQIIndex )( tIQICounter ).set( "leader_dof_dependencies",    "STRESS_DOF");
    tParameterList( tIQIIndex )( tIQICounter ).set( "vectorial_field_index",      0 );
    tParameterList( tIQIIndex )( tIQICounter ).set( "function_parameters",
        std::to_stri
         */

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        // tParameterList( 5 )( 0 ).set( "print_physics_model",      true );

        uint tFieldCounter = 0;

        tParameterList( 6 ).push_back( prm::create_fem_field_parameter_list() );
        tParameterList( 6 )( tFieldCounter ).set( "field_name", "FieldZienkiewiczZhu" );
        tParameterList( 6 )( tFieldCounter ).set( "field_entity_type", "ELEMENTAL" );
        tParameterList( 6 )( tFieldCounter ).set( "field_type", "FIELD_1" );
        tParameterList( 6 )( tFieldCounter ).set( "IQI_Name", "IQIBulkZienkiewiczZhu" );
        // tParameterList( 6 )( tFieldCounter ).set( "field_output_to_file",      "Field_Zienkiewicz_Zhu.hdf5" );
        tFieldCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

        /*
    // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
    //              Pseudoblock TFQMR, Seed GMRES, Seed CG
    tParameterlist( 0 )( 0 ).set( "Solver Type" ,  "GMRES" );

    // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
    sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
    tParameterlist( 0 )( 0 ).set( "Verbosity" , tVerbosity );

    // Maximum number of blocks in Krylov factorization
    tParameterlist( 0 )( 0 ).set( "Num Blocks", 500   );

    // Block size to be used by iterative solver
    tParameterlist( 0 )( 0 ).set( "Block Size", 1   );

    // Allowable Belos solver iterations
    tParameterlist( 0 )( 0 ).set( "Maximum Iterations" , 500 );

    // Allowable Belos solver iterations
    //tParameterlist( 0 )( 0 ).set( "Maximum Restarts" ,  );

    // Convergence criteria
    tParameterlist( 0 )( 0 ).set( "Convergence Tolerance" ,  1e-12 );

    // Preconditioner
    tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILU");
    tParameterlist( 0 )( 0 ).set( "fact: level-of-fill", 10  );

    //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILUT");
    //tParameterlist( 0 )( 0 ).set( "fact: ilut level-of-fill", 15.0 );
    //tParameterlist( 0 )( 0 ).set( "fact: drop tolerance", 1e-12 );
         */

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 ).resize( 6 );
        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 0
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 1 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();    // nonlinear algorithm index 1
        tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 2 )( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 2 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 3 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 3 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 3 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 3 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 4 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 4 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 4 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 4 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 4 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 )( 5 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 5 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 2 )( 5 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        tParameterlist( 2 )( 5 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 5 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        tParameterlist( 3 ).resize( 6 );
        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 0
        tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "THETA" );
        tParameterlist( 3 )( 0 ).set( "NLA_Secondary_DofTypes", "" );

        tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 1
        tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "PHID" );
        tParameterlist( 3 )( 1 ).set( "NLA_Secondary_DofTypes", "" );

        tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "UX,UY" );
        tParameterlist( 3 )( 2 ).set( "NLA_Secondary_DofTypes", "STRESS_DOF" );

        tParameterlist( 3 )( 3 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "4" );    // set nonlinear algorithm with index 0
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "STRESS_DOF" );
        tParameterlist( 3 )( 3 ).set( "NLA_Secondary_DofTypes", "" );

        tParameterlist( 3 )( 4 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 2
        tParameterlist( 3 )( 4 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 3 )( 4 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             // set nonlinear algorithm with index 1.
        tParameterlist( 3 )( 4 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2,3" );              // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 4 ).set( "NLA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );

        tParameterlist( 3 )( 5 ) = moris::prm::create_nonlinear_solver_parameter_list();    // nonlinear solver index 2
        tParameterlist( 3 )( 5 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 3 )( 5 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             // set nonlinear algorithm with index 1.
        tParameterlist( 3 )( 5 ).set( "NLA_Sub_Nonlinear_Solver", "2,3" );                  // set sub nonlinear solvers with index 0 and 1
        tParameterlist( 3 )( 5 ).set( "NLA_DofTypes", "UX,UY;STRESS_DOF" );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 4 );
        // tParameterlist( 4 )( 0 ).set("TSA_nonlinear_solver_for_adjoint_solve", 5 );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
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
        tParameterlist( 0 )( 0 ).set( "STRESS_DOF", 0 );
        tParameterlist( 0 )( 0 ).set( "THETA", 0 );
        tParameterlist( 0 )( 0 ).set( "PHID", 0 );
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
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomain + "," + tFrameSets );

        tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,VOL,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT,TRUNCPHID,STRESS,ZienKiewiczZhu" ) );
        tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,ELEMENTAL_AVG" ) );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod,IQITruncPHID,IQIBulkStress,IQIBulkZienkiewiczZhu" ) );

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
        tParameterlist( 0 )( 0 ).set( "mode", "advanced" );
        tParameterlist( 0 )( 0 ).set( "remeshing_field_names", "Level_Set_Field,Outer_shape,Inner_Holes,Outer_Holes,FieldZienkiewiczZhu" );
        tParameterlist( 0 )( 0 ).set( "refinement_function_name", "refinement_function,refinement_function,refinement_function,refinement_function,refinement_function_stress" );

        tParameterlist( 0 )( 0 ).set( "remeshing_copy_old_pattern_to_pattern", "0, 3; 1, 4; 2, 5" );
        tParameterlist( 0 )( 0 ).set( "remeshing_refinement_pattern", "0;0;0;0;0,1" );

        tParameterlist( 0 )( 0 ).set( "minimum_refinement_level", "0, 0, 30, 10, 100, 0, 500; 1, 0, 30, 0, 500;2, 0, 500" );

        tParameterlist( 0 )( 0 ).set( "output_meshes", true );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
