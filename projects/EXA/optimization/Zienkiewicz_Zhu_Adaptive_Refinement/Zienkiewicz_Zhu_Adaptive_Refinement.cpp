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
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );
        aParameterLists.set( "reinitialize_interface_iter", 18 );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "step_size", tMMAStepSize );
        aParameterLists.set( "penalty", tMMAPenalty );
        aParameterLists.set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        aParameterLists.set( "restart_index", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", "60,30" );
        aParameterLists.set( "domain_dimensions", "8,4" );
        aParameterLists.set( "domain_offset", "-4.0,-2.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "2" );
        aParameterLists.set( "lagrange_pattern", "0" );

        aParameterLists.set( "bspline_orders", "1,2" );
        aParameterLists.set( "bspline_pattern", "1,2" );

        aParameterLists.set( "initial_refinement", "0,0" );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        // aParameterLists.set( "write_lagrange_output_mesh", "HMRLagrangeMesh.vtk" );

        // aParameterLists.set( "use_refine_low_level_elements", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy", "IQIBulkVolume", "IQIPerimeter_InterfaceVoid", "IQIHeatMethodPenalty", "IQIMaxStress" );
        // aParameterLists.set("output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );
        // FIXME     this has to change
        // aParameterLists.set("PDV_types"         , "LS1");

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
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );
        aParameterLists.set( "print_phase_table", true );

        // outer frame
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "general_shape_outer" );
        // aParameterLists.set( "number_of_refinements", "0,0" );
        // aParameterLists.set( "refinement_mesh_index", "0,1" );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 1.2 );
        aParameterLists.set( "name", "Outer_shape" );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Hole_inner" );
        // aParameterLists.set( "number_of_refinements", "0,0" );
        // aParameterLists.set( "refinement_mesh_index", "0,1" );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 1.2 );
        aParameterLists.set( "name", "Inner_Holes" );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Hole_outer" );
        // aParameterLists.set( "number_of_refinements", "0,0" );
        // aParameterLists.set( "refinement_mesh_index", "0,1" );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 1.2 );
        aParameterLists.set( "name", "Outer_Holes" );

        // initialize fins as swiss cheese geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( gen::Field_Type::USER_DEFINED );
        aParameterLists.set( "field_function_name", "Const_Geometry" );
        aParameterLists.set( "name", "Level_Set_Field" );
        // aParameterLists.set( "number_of_refinements", "0,0" );
        // aParameterLists.set( "refinement_mesh_index", "0,1" );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 1.2 );
        aParameterLists.set( "use_multilinear_interpolation", false );

        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 1 );
            aParameterLists.set( "discretization_lower_bound", tBsplineLimitBottom );
            aParameterLists.set( "discretization_upper_bound", tBsplineLimitTop );
        }

        uint tParamCounter = 0;
        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );
        aParameterLists.set( "name", "LvL_Set_Field" );
        aParameterLists.set( "dependencies", "Level_Set_Field" );
        aParameterLists.set( "scaling_factor", 1.0 );
        aParameterLists.set( "pdv_type", "LS1" );
        // aParameterLists.set("discretization_mesh_index",   -1);
        // aParameterLists.set("discretization_lower_bound", 0.001);
        // aParameterLists.set("discretization_upper_bound", 1.0);

        if ( tUseGhost )
        {
            aParameterLists.set( "pdv_mesh_set_names", tTotalDomainAGhost );
        }
        else
        {
            aParameterLists.set( "pdv_mesh_set_names", tTotalDomain );
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
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tDensityFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerDensityFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityFrame" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0,3.0" );
        aParameterLists.set( "value_function", "tYoungsFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerYoungsFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsFrame" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-8" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSelect" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Select" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropUnitValue" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for Theta
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityTheta" );
        aParameterLists.set( "function_parameters", tDensityTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityTheta" );
        aParameterLists.set( "function_parameters", tCapacityTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescTheta" );
        aParameterLists.set( "function_parameters", tPrescTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for phi problem

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityPhi" );
        aParameterLists.set( "function_parameters", tDensityPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityPhi" );
        aParameterLists.set( "function_parameters", tCapacityPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescPhi" );
        aParameterLists.set( "function_parameters", tPrescPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropEigenStrainPhi" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetConst" );
        aParameterLists.set( "value_function", "Func_Const" );
        aParameterLists.set( "function_parameters", "1.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists.set( "value_function", "Func_Const" );
        aParameterLists.set( "function_parameters", "1.0;1.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSet" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetReal" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFuncReal" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradx" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTruncPHID" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "PhiD_Prop" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIsoFrame" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungsFrame,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model - Theta problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionTheta" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionPhi" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungsFrame,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhostFrame" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.001" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungsFrame,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost_VM" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::to_string( tStressGhost ) + "/0" );
        aParameterLists.set( "follower_properties", std::string( "PropYoungsFrame,Material" ) );
        aParameterLists.set( "leader_properties", "PropUnitValue,Material" );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.001" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG  - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionThetaBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityThetaInterior" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID;UX,UY;STRESS_DOF" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists.set( "mesh_set_names", tBulkSets );
        aParameterLists.set( "time_continuity", true );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityThetaVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists.set( "mesh_set_names", tVoidSets );
        aParameterLists.set( "time_continuity", true );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkUFrame" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet;PropSelect,Select" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoFrame,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tDirichletSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tLoadSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_dv_dependencies", "LS1" );
        aParameterLists.set( "follower_dv_dependencies", "LS1" );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMStrucLinIso,ElastLinIso" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMStrucLinIsoFrame,ElastLinIso" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tFrameInteriorDSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGStress" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_VON_MISES_STRESS );
        aParameterLists.set( "dof_residual", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostDisp" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "leader_dv_dependencies", "LS1" );
            aParameterLists.set( "follower_dv_dependencies", "LS1" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostFrame" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhostFrame,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tFrameGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGSStress" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "STRESS_DOF" );
            aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
            aParameterLists.set( "follower_dof_dependencies", "STRESS_DOF" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_VM,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostPhi" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "PHID" );
            aParameterLists.set( "leader_dof_dependencies", "PHID" );
            aParameterLists.set( "follower_dof_dependencies", "PHID" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostTemp" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "THETA" );
            aParameterLists.set( "leader_dof_dependencies", "THETA" );
            aParameterLists.set( "follower_dof_dependencies", "THETA" );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGPTemp,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost + "," + tVoidGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets + "," + tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSets + "," + tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStress" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists.set( "dof_quantity", "STRESS_DOF" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIsoFrame,Elast" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        // smooth fuction
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSet" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSetReal,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITruncPHID" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropTruncPHID,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSetHeatMethod" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSet,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Nodal PHID IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkPHID" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is constant
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIH1ErrorConst" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        aParameterLists.set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIH1Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatMethodPenalty" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + " / 0.0/ 0.2 / 1.0 / 1.0 " );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIMaxStress" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::MAX_DOF );
        aParameterLists.set( "dof_quantity", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "function_parameters", std::to_string( 5.0 ) + "/" + std::to_string( tVMMaxDofNValue ) + "/" + std::to_string( tVMShift ) + "/1" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkZienkiewiczZhu" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::ZIENKIEWICZ_ZHU_VON_MISES_STRESS );
        aParameterLists.set( "dof_quantity", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        /*
    aParameterLists( tIQIIndex ).push_back( prm::create_IQI_parameter_list() );
    aParameterLists.set( "IQI_name",                   "IQIMaxStress2");
    aParameterLists.set( "leader_phase_name",          "PhaseLBracket" );
    aParameterLists.set( "neighbor_phases",            "PhaseVoid2" );
    aParameterLists.set( "IQI_bulk_type",               fem::Element_Type::SIDESET ) ;
    aParameterLists.set( "IQI_type",                    fem::IQI_Type::MAX_DOF ) ;
    aParameterLists.set( "dof_quantity",               "STRESS_DOF");
    aParameterLists.set( "leader_dof_dependencies",    "STRESS_DOF");
    aParameterLists.set( "vectorial_field_index",      0 );
    aParameterLists.set( "function_parameters",
        std::to_stri
         */

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
        // aParameterLists.set( "print_physics_model",      true );

        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", "FieldZienkiewiczZhu" );
        aParameterLists.set( "field_entity_type", "ELEMENTAL" );
        aParameterLists.set( "field_type", "FIELD_1" );
        aParameterLists.set( "IQI_Name", "IQIBulkZienkiewiczZhu" );
        // aParameterLists.set( "field_output_to_file",      "Field_Zienkiewicz_Zhu.hdf5" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        /*
    // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
    //              Pseudoblock TFQMR, Seed GMRES, Seed CG
    aParameterLists.set( "Solver Type" ,  "GMRES" );

    // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
    sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
    aParameterLists.set( "Verbosity" , tVerbosity );

    // Maximum number of blocks in Krylov factorization
    aParameterLists.set( "Num Blocks", 500   );

    // Block size to be used by iterative solver
    aParameterLists.set( "Block Size", 1   );

    // Allowable Belos solver iterations
    aParameterLists.set( "Maximum Iterations" , 500 );

    // Allowable Belos solver iterations
    //aParameterLists.set( "Maximum Restarts" ,  );

    // Convergence criteria
    aParameterLists.set( "Convergence Tolerance" ,  1e-12 );

    // Preconditioner
    aParameterLists.set( "ifpack_prec_type",  "ILU");
    aParameterLists.set( "fact: level-of-fill", 10  );

    //aParameterLists.set( "ifpack_prec_type",  "ILUT");
    //aParameterLists.set( "fact: ilut level-of-fill", 15.0 );
    //aParameterLists.set( "fact: drop tolerance", 1e-12 );
         */

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 1
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "THETA" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 1
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "PHID" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "STRESS_DOF" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "4" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "STRESS_DOF" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 2
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,2,3" );     // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 2
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "2,3" );         // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "UX,UY;STRESS_DOF" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 4 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "THETA;PHID;UX,UY;STRESS_DOF" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
        aParameterLists.set( "STRESS_DOF", 0 );
        aParameterLists.set( "THETA", 0 );
        aParameterLists.set( "PHID", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain + "," + tFrameSets );

        aParameterLists.set( "Field_Names", std::string( "UX,UY,VOL,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT,TRUNCPHID,STRESS,ZienKiewiczZhu" ) );
        aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,ELEMENTAL_AVG" ) );
        aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod,IQITruncPHID,IQIBulkStress,IQIBulkZienkiewiczZhu" ) );

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_remeshing_parameterlist( aParameterLists( 0 )( 0 ) );
        aParameterLists.set( "mode", "advanced" );
        aParameterLists.set( "remeshing_field_names", "Level_Set_Field,Outer_shape,Inner_Holes,Outer_Holes,FieldZienkiewiczZhu" );
        aParameterLists.set( "refinement_function_name", "refinement_function,refinement_function,refinement_function,refinement_function,refinement_function_stress" );

        aParameterLists.set( "remeshing_copy_old_pattern_to_pattern", "0, 3; 1, 4; 2, 5" );
        aParameterLists.set( "remeshing_refinement_pattern", "0;0;0;0;0,1" );

        aParameterLists.set( "minimum_refinement_level", "0, 0, 30, 10, 100, 0, 500; 1, 0, 30, 0, 500;2, 0, 500" );

        aParameterLists.set( "output_meshes", true );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
