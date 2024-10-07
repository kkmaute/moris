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
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            mtk::Cell*                           aElement,
            const std::shared_ptr< mtk::Field >& aField,
            uint                                 tActivationPattern,
            uint&                                aMaxLevel )
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
            mtk::Cell*                           aElement,
            const std::shared_ptr< mtk::Field >& aField,
            uint                                 tActivationPattern,
            uint&                                aMaxLevel )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            const Vector< real >&   aGeometryParameters )
    {
        return 2.0;
    }

    //------------------------------------------------------------------------------
    // Hole pattern
    real
    general_shape_outer(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )

    {
        real tLSval = 0.0;

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
            const Vector< real >&   aGeometryParameters )
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
            const Vector< real >&   aGeometryParameters )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 ) / 2.40;

        real tDensity = ( tLevelSet - phi_sh_scale ) / ( 1 - phi_sh_scale );

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix = tHCT * std::pow( tDensity, tBeta );
    }

    //------------------------------------------------------------------------------
    void
    tDerYoungsFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        // FIXME density shift missing

        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() *    //
                      tBeta * tHCT( 0 ) * std::pow( tDensity, tBeta - 1 ) / ( 1 - phi_sh );
    }

    //------------------------------------------------------------------------------

    void
    tDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 ) / 2.40;

        real tDensity = ( tLevelSet - phi_sh_scale ) / ( 1 - phi_sh_scale );

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix.set_size( 1, 1, tDensity );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() / ( 1 - phi_sh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tLevelSetFuncReal(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val();
    }

    /* ------------------------------------------------------------------------ */

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // Matrix< DDRMat > tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

        // bool tBool = is_in_vicinity_of_load( tCoords );

        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
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

        value = value - phi_sh;

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
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

        value = value - phi_sh;

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->dnNdxn( 1 );
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
Matrix<DDRMat> compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
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
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
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
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
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
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

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
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );

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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", tIsOpt );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );
        aParameterLists( 0 ).set( "restart_file", "" );
        aParameterLists( 0 ).set( "reinitialize_interface_iter", 18 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
        aParameterLists( 2 ).set( "step_size", tMMAStepSize );
        aParameterLists( 2 ).set( "penalty", tMMAPenalty );
        aParameterLists( 2 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        aParameterLists( 2 ).set( "restart_index", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", "60,30" );
        aParameterLists( 0 ).set( "domain_dimensions", "8,4" );
        aParameterLists( 0 ).set( "domain_offset", "-4.0,-2.0" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "2" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );

        aParameterLists( 0 ).set( "bspline_orders", "1,2" );
        aParameterLists( 0 ).set( "bspline_pattern", "1,2" );

        aParameterLists( 0 ).set( "initial_refinement", "0,0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0,1" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0,1" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 1 );
        aParameterLists( 0 ).set( "staircase_buffer", 1 );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        // aParameterLists( 0 ).set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        // aParameterLists( 0 ).set( "use_refine_low_level_elements", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0,1" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy", "IQIBulkVolume", "IQIPerimeter_InterfaceVoid", "IQIHeatMethodPenalty", "IQIMaxStress" );
        // aParameterLists(0).set("output_mesh_file", tGENOutputFile );
        aParameterLists( 0 ).set( "time_offset", 10.0 );
        // FIXME     this has to change
        // aParameterLists( 0 ).set("PDV_types"         , "LS1");

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
        aParameterLists( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );
        aParameterLists( 0 ).set( "print_phase_table", true );

        // outer frame
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "general_shape_outer" );
        // aParameterLists( 1 ).set( "number_of_refinements", "0,0" );
        // aParameterLists( 1 ).set( "refinement_mesh_index", "0,1" );
        aParameterLists( 1 ).set( "isocontour_tolerance", 10e-14 );
        aParameterLists( 1 ).set( "isocontour_threshold", 1.2 );
        aParameterLists( 1 ).set( "name", "Outer_shape" );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Hole_inner" );
        // aParameterLists( 1 ).set( "number_of_refinements", "0,0" );
        // aParameterLists( 1 ).set( "refinement_mesh_index", "0,1" );
        aParameterLists( 1 ).set( "isocontour_tolerance", 10e-14 );
        aParameterLists( 1 ).set( "isocontour_threshold", 1.2 );
        aParameterLists( 1 ).set( "name", "Inner_Holes" );

        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Hole_outer" );
        // aParameterLists( 1 ).set( "number_of_refinements", "0,0" );
        // aParameterLists( 1 ).set( "refinement_mesh_index", "0,1" );
        aParameterLists( 1 ).set( "isocontour_tolerance", 10e-14 );
        aParameterLists( 1 ).set( "isocontour_threshold", 1.2 );
        aParameterLists( 1 ).set( "name", "Outer_Holes" );

        // initialize fins as swiss cheese geometry
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Const_Geometry" );
        aParameterLists( 1 ).set( "name", "Level_Set_Field" );
        // aParameterLists( 1 ).set( "number_of_refinements", "0,0" );
        // aParameterLists( 1 ).set( "refinement_mesh_index", "0,1" );
        aParameterLists( 1 ).set( "isocontour_tolerance", 10e-14 );
        aParameterLists( 1 ).set( "isocontour_threshold", 1.2 );
        aParameterLists( 1 ).set( "use_multilinear_interpolation", false );

        if ( tIsOpt )
        {
            aParameterLists( 1 ).set( "discretization_mesh_index", 1 );
            aParameterLists( 1 ).set( "discretization_lower_bound", tBsplineLimitBottom );
            aParameterLists( 1 ).set( "discretization_upper_bound", tBsplineLimitTop );
        }

        uint tParamCounter = 0;
        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gen_property_parameter_list( gen::Field_Type::SCALED_FIELD ) );
        aParameterLists( 2 ).set( "name", "LvL_Set_Field" );
        aParameterLists( 2 ).set( "dependencies", "Level_Set_Field" );
        aParameterLists( 2 ).set( "scaling_factor", 1.0 );
        aParameterLists( 2 ).set( "pdv_type", "LS1" );
        // aParameterLists( 2 ).set("discretization_mesh_index",   -1);
        // aParameterLists( 2 ).set("discretization_lower_bound", 0.001);
        // aParameterLists( 2 ).set("discretization_upper_bound", 1.0);

        if ( tUseGhost )
        {
            aParameterLists( 2 ).set( "pdv_mesh_set_names", tTotalDomainAGhost );
        }
        else
        {
            aParameterLists( 2 ).set( "pdv_mesh_set_names", tTotalDomain );
        }

        tParamCounter++;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tDensityFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerDensityFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        // create parameter list for property 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityFrame" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", "1.0,3.0" );
        aParameterLists( 0 ).set( "value_function", "tYoungsFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerYoungsFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungsFrame" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", "1.0e-8" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichletU" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropSelect" );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Select" );

        // create parameter list for property 10
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropTraction" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropUnitValue" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity" );
        aParameterLists( 0 ).set( "function_parameters", tConductivity );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties for Theta
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityTheta" );
        aParameterLists( 0 ).set( "function_parameters", tDensityTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityTheta" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPrescTheta" );
        aParameterLists( 0 ).set( "function_parameters", tPrescTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties for phi problem

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityPhi" );
        aParameterLists( 0 ).set( "function_parameters", tDensityPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityPhi" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPrescPhi" );
        aParameterLists( 0 ).set( "function_parameters", tPrescPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropEigenStrainPhi" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetConst" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );
        aParameterLists( 0 ).set( "function_parameters", "1.0;1.0" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSet" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetReal" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetFuncReal" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetGradx" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropTruncPHID" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "PhiD_Prop" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIsoFrame" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungsFrame,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model - Theta problem
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionTheta" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionPhi" );
        aParameterLists( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsFrame,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPGhost" ) );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", std::string( "0.005" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPGhostFrame" ) );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", std::string( "0.001" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungsFrame,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGhost_VM" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", std::to_string( tStressGhost ) + "/0" );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );
        aParameterLists( 2 ).set( "leader_properties", "PropUnitValue,Material" );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists( 2 ).set( "function_parameters", "0.001" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG  - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionThetaBulk" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityThetaInterior" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID;UX,UY;STRESS_DOF" );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", tBulkSets );
        aParameterLists( 3 ).set( "time_continuity", true );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityThetaVoid" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidSets );
        aParameterLists( 3 ).set( "time_continuity", true );

        // create IWG - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tBulkSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkUFrame" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletU" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichletU,Dirichlet;PropSelect,Select" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tDirichletSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTraction" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_properties", "PropTraction,Traction" );
        aParameterLists( 3 ).set( "mesh_set_names", tLoadSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists( 3 ).set( "leader_dv_dependencies", "LS1" );
        aParameterLists( 3 ).set( "follower_dv_dependencies", "LS1" );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMStrucLinIso,ElastLinIso" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMStrucLinIsoFrame,ElastLinIso" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameInteriorDSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGStress" );
        aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_VON_MISES_STRESS );
        aParameterLists( 3 ).set( "dof_residual", "STRESS_DOF" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists( 3 ).set( "mesh_set_names", tBulkSets );

        if ( tUseGhost )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostDisp" ) );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dv_dependencies", "LS1" );
            aParameterLists( 3 ).set( "follower_dv_dependencies", "LS1" );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostFrame" ) );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhostFrame,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tFrameGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGSStress" ) );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "STRESS_DOF" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "STRESS_DOF" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "STRESS_DOF" );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhost_VM,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostPhi" ) );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "PHID" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "PHID" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "PHID" );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostTemp" ) );
            aParameterLists( 3 ).set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists( 3 ).set( "dof_residual", "THETA" );
            aParameterLists( 3 ).set( "leader_dof_dependencies", "THETA" );
            aParameterLists( 3 ).set( "follower_dof_dependencies", "THETA" );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUX" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets + "," + tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUY" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets + "," + tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStress" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists( 4 ).set( "dof_quantity", "STRESS_DOF" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIsoFrame,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists( 4 ).set( "mesh_set_names", tInterfaceVoidSets );

        // smooth fuction
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQILevelSet" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSetReal,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQITruncPHID" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( 4 ).set( "leader_properties", "PropTruncPHID,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQILevelSetHeatMethod" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSet,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "THETA" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Nodal PHID IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkPHID" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is constant
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIH1ErrorConst" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIH1Error" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIHeatMethodPenalty" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + " / 0.0/ 0.2 / 1.0 / 1.0 " );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxStress" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists( 4 ).set( "dof_quantity", "STRESS_DOF" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "function_parameters", std::to_string( 5.0 ) + "/" + std::to_string( tVMMaxDofNValue ) + "/" + std::to_string( tVMShift ) + "/1" );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkZienkiewiczZhu" );
        aParameterLists( 4 ).set( "IQI_type", fem::IQI_Type::ZIENKIEWICZ_ZHU_VON_MISES_STRESS );
        aParameterLists( 4 ).set( "dof_quantity", "STRESS_DOF" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists( 4 ).set( "mesh_set_names", tBulkSets );

        /*
    aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
    aParameterLists( tIQIIndex ).set( "IQI_name",                   "IQIMaxStress2");
    aParameterLists( tIQIIndex ).set( "leader_phase_name",          "PhaseLBracket" );
    aParameterLists( tIQIIndex ).set( "neighbor_phases",            "PhaseVoid2" );
    aParameterLists( tIQIIndex ).set( "IQI_bulk_type",               fem::Element_Type::SIDESET ) ;
    aParameterLists( tIQIIndex ).set( "IQI_type",                    fem::IQI_Type::MAX_DOF ) ;
    aParameterLists( tIQIIndex ).set( "dof_quantity",               "STRESS_DOF");
    aParameterLists( tIQIIndex ).set( "leader_dof_dependencies",    "STRESS_DOF");
    aParameterLists( tIQIIndex ).set( "vectorial_field_index",      0 );
    aParameterLists( tIQIIndex ).set( "function_parameters",
        std::to_stri
         */

        // create computation  parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
        // aParameterLists( 5 ).set( "print_physics_model",      true );

        aParameterLists( 6 ).add_parameter_list( prm::create_fem_field_parameter_list() );
        aParameterLists( 6 ).set( "field_name", "FieldZienkiewiczZhu" );
        aParameterLists( 6 ).set( "field_entity_type", "ELEMENTAL" );
        aParameterLists( 6 ).set( "field_type", "FIELD_1" );
        aParameterLists( 6 ).set( "IQI_Name", "IQIBulkZienkiewiczZhu" );
        // aParameterLists( 6 ).set( "field_output_to_file",      "Field_Zienkiewicz_Zhu.hdf5" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        /*
    // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
    //              Pseudoblock TFQMR, Seed GMRES, Seed CG
    aParameterLists( 0 ).set( "Solver Type" ,  "GMRES" );

    // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
    sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
    aParameterLists( 0 ).set( "Verbosity" , tVerbosity );

    // Maximum number of blocks in Krylov factorization
    aParameterLists( 0 ).set( "Num Blocks", 500   );

    // Block size to be used by iterative solver
    aParameterLists( 0 ).set( "Block Size", 1   );

    // Allowable Belos solver iterations
    aParameterLists( 0 ).set( "Maximum Iterations" , 500 );

    // Allowable Belos solver iterations
    //aParameterLists( 0 ).set( "Maximum Restarts" ,  );

    // Convergence criteria
    aParameterLists( 0 ).set( "Convergence Tolerance" ,  1e-12 );

    // Preconditioner
    aParameterLists( 0 ).set( "ifpack_prec_type",  "ILU");
    aParameterLists( 0 ).set( "fact: level-of-fill", 10  );

    //aParameterLists( 0 ).set( "ifpack_prec_type",  "ILUT");
    //aParameterLists( 0 ).set( "fact: ilut level-of-fill", 15.0 );
    //aParameterLists( 0 ).set( "fact: drop tolerance", 1e-12 );
         */

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 1
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 1
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );    // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "PHID" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "STRESS_DOF" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "4" );    // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "STRESS_DOF" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 2
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2,3" );     // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 2
        aParameterLists( 3 ).set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "2,3" );         // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY;STRESS_DOF" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 4 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "UX", 0 );
        aParameterLists( 0 ).set( "UY", 0 );
        aParameterLists( 0 ).set( "STRESS_DOF", 0 );
        aParameterLists( 0 ).set( "THETA", 0 );
        aParameterLists( 0 ).set( "PHID", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists( 0 ).set( "Set_Names", tTotalDomain + "," + tFrameSets );

        aParameterLists( 0 ).set( "Field_Names", std::string( "UX,UY,VOL,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT,TRUNCPHID,STRESS,ZienKiewiczZhu" ) );
        aParameterLists( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,ELEMENTAL_AVG" ) );
        aParameterLists( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod,IQITruncPHID,IQIBulkStress,IQIBulkZienkiewiczZhu" ) );

        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_remeshing_parameterlist( aParameterLists( 0 )( 0 ) );
        aParameterLists( 0 ).set( "mode", "advanced" );
        aParameterLists( 0 ).set( "remeshing_field_names", "Level_Set_Field,Outer_shape,Inner_Holes,Outer_Holes,FieldZienkiewiczZhu" );
        aParameterLists( 0 ).set( "refinement_function_name", "refinement_function,refinement_function,refinement_function,refinement_function,refinement_function_stress" );

        aParameterLists( 0 ).set( "remeshing_copy_old_pattern_to_pattern", "0, 3; 1, 4; 2, 5" );
        aParameterLists( 0 ).set( "remeshing_refinement_pattern", "0;0;0;0;0,1" );

        aParameterLists( 0 ).set( "minimum_refinement_level", "0, 0, 30, 10, 100, 0, 500; 1, 0, 30, 0, 500;2, 0, 500" );

        aParameterLists( 0 ).set( "output_meshes", true );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
