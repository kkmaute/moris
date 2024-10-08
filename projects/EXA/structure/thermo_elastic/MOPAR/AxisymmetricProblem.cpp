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
#include "fn_PRM_STK_Parameters.hpp"
#include "fn_equal_to.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "fn_interp1.hpp"
#include "paths.hpp"

#include "cl_Ascii.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_assert.hpp"
#include "catch.hpp"

#include "AztecOO.h"

extern moris::Logger gLogger;

#ifdef __cplusplus
extern "C" {

#endif
//------------------------------------------------------------------------------
namespace moris
{
    // Input file
    std::string tOutputFileName = "AxisymmetricProblem.exo";
    std::string tRoot           = get_base_moris_dir();
    std::string tFieldRefPath   = tRoot + "projects/EXA/structure/thermo_elastic/MOPAR/traction.e";

    bool byGeometry = false;

    /* ------------------------------------------------------------------------ */
    // material parameters (Quasi-isotropic carbon fiber)

    // thermally dependent modulus {{T1,E1},{T2,E2},etc)
    Matrix< DDRMat > tInterpTemp = { { 0, 293, 1000, 10000 } };

    Matrix< DDRMat > tInterpE_c  = { { 50.00 * 1e9, 50.00 * 1e9, 50.00 * 1e9, 50.00 * 1e9 } };
    Matrix< DDRMat > tInterpE_al = { { 6.89 * 1e9, 6.89 * 1e9, 6.89 * 1e9, 6.89 * 1e9 } };

    real tPois_c  = 0.10;
    real tPois_al = 0.33;

    real tDens_c  = 1.879 * 1e3;
    real tDens_al = 2.700 * 1e3;

    real tCTE_c  = 2.0 * 1e-6;
    real tCTE_al = 23.6 * 1e-6;

    std::string tBedding = std::to_string( 1.0 * 1e-6 );

    std::string tRefTemp  = "293.0";
    std::string tPropTemp = "0.0";

    //------------------------------------------------------------------------------

    // FUNCTIONS

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Density(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // aluminum
        if ( 0.0160091 / 0.09 * ( tX( 0 ) - 0.01 ) - tX( 1 ) > 0 )
        {
            aPropMatrix = { { tDens_al } };
        }
        // ablator
        else
        {
            aPropMatrix = { { tDens_c } };
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Poisson(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // aluminum
        if ( 0.0160091 / 0.09 * ( tX( 0 ) - 0.01 ) - tX( 1 ) > 0 )
        {
            aPropMatrix = { { tPois_al } };
        }
        // ablator
        else
        {
            aPropMatrix = { { tPois_c } };
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_CTE(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // aluminum
        if ( 0.0160091 / 0.09 * ( tX( 0 ) - 0.01 ) - tX( 1 ) > 0 )
        {
            aPropMatrix = { { tCTE_al } };
        }
        // ablator
        else
        {
            aPropMatrix = { { tCTE_c } };
        }
    }

    /* ------------------------------------------------------------------------ */

    // Young's Modulus Interpolator
    void
    Func_Thermal_Modulus_Block(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get temperature
        Matrix< DDRMat > tTemperature =
                aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();

        real tNewE = -1.0;

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // linearly interpolating modulus using a given temperature.
        for ( uint iT = 0; iT < tInterpTemp.numel(); iT++ )
        {
            if ( tTemperature( 0, 0 ) < tInterpTemp( iT ) )
            {
                // aluminum
                if ( aParameters( 0 )( 0 ) < 1.5 )
                {
                    tNewE = tInterpE_c( iT - 1 ) + ( tInterpE_c( iT ) - tInterpE_c( iT - 1 ) ) * ( tTemperature( 0, 0 ) - tInterpTemp( iT - 1 ) ) / ( tInterpTemp( iT ) - tInterpTemp( iT - 1 ) );
                }
                // ablator
                else
                {
                    tNewE = tInterpE_al( iT - 1 ) + ( tInterpE_al( iT ) - tInterpE_al( iT - 1 ) ) * ( tTemperature( 0, 0 ) - tInterpTemp( iT - 1 ) ) / ( tInterpTemp( iT ) - tInterpTemp( iT - 1 ) );
                }
                break;
            }
        }

        MORIS_ERROR( tNewE > 0, "Negative Young's modulus." );

        aPropMatrix = { { tNewE } };
    }

    /* ------------------------------------------------------------------------ */

    // Young's Modulus Interpolator
    void
    Func_Thermal_Modulus_Geo(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get temperature
        Matrix< DDRMat > tTemperature =
                aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();

        real tNewE = -1.0;

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // linearly interpolating modulus using a given temperature.
        for ( uint iT = 0; iT < tInterpTemp.numel(); iT++ )
        {
            if ( tTemperature( 0, 0 ) < tInterpTemp( iT ) )
            {
                // aluminum
                if ( 0.0160091 / 0.09 * ( tX( 0 ) - 0.01 ) - tX( 1 ) > 0 )
                {
                    tNewE = tInterpE_c( iT - 1 ) + ( tInterpE_c( iT ) - tInterpE_c( iT - 1 ) ) * ( tTemperature( 0, 0 ) - tInterpTemp( iT - 1 ) ) / ( tInterpTemp( iT ) - tInterpTemp( iT - 1 ) );
                }
                // ablator
                else
                {
                    tNewE = tInterpE_al( iT - 1 ) + ( tInterpE_al( iT ) - tInterpE_al( iT - 1 ) ) * ( tTemperature( 0, 0 ) - tInterpTemp( iT - 1 ) ) / ( tInterpTemp( iT ) - tInterpTemp( iT - 1 ) );
                }
                break;
            }
        }

        MORIS_ERROR( tNewE > 0, "Negative Young's modulus." );

        aPropMatrix = { { tNewE } };
    }

    /* ------------------------------------------------------------------------ */

    // calculates the outward normal vector for distance from rotation axis to point of interest
    // and the radius so aPropMatrix = {{2*pi*r},{r},{n1},{n2}}
    void
    Func_AxiRotation(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ASSERT( aParameters( 0 ).n_cols() == 2 and aParameters( 0 ).n_rows() == 2,
                "Axisymmetric rotation axis incorrectly defined. use {{x1,y1},{x2,y2}}." );

        // vector from point 1 to point 2
        moris::Matrix< moris::DDRMat > tRotVec = { //
            { aParameters( 0 )( 1, 0 ) - aParameters( 0 )( 0, 0 ) },
            { aParameters( 0 )( 1, 1 ) - aParameters( 0 )( 0, 1 ) }
        };

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // vector from point 1 to location of interest
        moris::Matrix< moris::DDRMat > tPntVec = { //
            { tX( 0 ) - aParameters( 0 )( 0, 0 ) },
            { tX( 1 ) - aParameters( 0 )( 0, 1 ) }
        };

        // minimum distance vector from line to point of interest
        moris::Matrix< moris::DDRMat > tRadVec =                 //
                tPntVec - dot( tPntVec, tRotVec ) * tRotVec /    //
                                  ( tRotVec( 0 ) * tRotVec( 0 ) + tRotVec( 1 ) * tRotVec( 1 ) );

        // init aPropMatrixSize
        aPropMatrix.set_size( 4, 1, 0.0 );

        // radius
        aPropMatrix( 1 ) = std::max( norm( tRadVec ), MORIS_REAL_EPS );

        // 2*pi*r for IWG
        aPropMatrix( 0 ) = aPropMatrix( 1 ) * 4.0 * std::acos( 0 );

        // outward radial normal vector
        aPropMatrix( { 2, 3 } ) = tRadVec / aPropMatrix( 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Field(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_Import_Traction(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // Input file
        std::string tBaseDir      = get_base_moris_dir();
        std::string tDataFilePath = tBaseDir + "projects/EXA/structure/thermo_elastic/MOPAR/wall_data.txt";

        // create ascii class
        Ascii tAscii = Ascii( tDataFilePath, FileMode::OPEN_RDONLY );

        // get num rows and columns of data file
        uint             tNumRows = tAscii.length();
        Matrix< DDRMat > tRow;
        string_to_matrix( tAscii.line( 0 ), tRow );
        uint tNumCols = tRow.n_cols();

        // set size of returned matrix
        Matrix< DDRMat > tDataMatrix( tNumRows, tNumCols );

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // setting prop matrix
        aPropMatrix.set_size( 2, 1 );

        real tTractionX0 = 0.0;
        real tTractionY0 = 0.0;
        real tTractionX1 = 0.0;
        real tTractionY1 = 0.0;

        for ( uint i = 0; i < tNumRows; i++ )
        {
            // storing ascii row to matrix
            string_to_matrix( tAscii.line( i ), tRow );

            // are the number of columns the same?
            MORIS_ASSERT( tRow.n_cols() == tNumCols, "Inconsistent number of columns in input file" );

            // store data at row of matrix
            tDataMatrix( { i, i }, { 0, tNumCols - 1 } ) = tRow( { 0, 0 }, { 0, tNumCols - 1 } );

            if ( i == 0 or i == tNumRows - 1 )
            {
                aPropMatrix( 0 ) = ( tDataMatrix( i, 14 ) - tDataMatrix( i, 13 ) ) * tDataMatrix( i, 3 ) + tDataMatrix( i, 15 ) * tDataMatrix( i, 4 );
                aPropMatrix( 1 ) = ( tDataMatrix( i, 17 ) - tDataMatrix( i, 13 ) ) * tDataMatrix( i, 4 ) + tDataMatrix( i, 15 ) * tDataMatrix( i, 3 );
                break;
            }

            // find 2 points to use linear interpolation.
            // this structure is specific to the data used
            if ( tDataMatrix( i, 2 ) > tX( 1 ) )
            {
                tTractionX0 = ( tDataMatrix( i - 1, 14 ) - tDataMatrix( i - 1, 13 ) ) * tDataMatrix( i - 1, 3 ) + tDataMatrix( i - 1, 15 ) * tDataMatrix( i - 1, 4 );
                tTractionY0 = ( tDataMatrix( i - 1, 17 ) - tDataMatrix( i - 1, 13 ) ) * tDataMatrix( i - 1, 4 ) + tDataMatrix( i - 1, 15 ) * tDataMatrix( i - 1, 3 );
                tTractionX1 = ( tDataMatrix( i, 14 ) - tDataMatrix( i, 13 ) ) * tDataMatrix( i, 3 ) + tDataMatrix( i, 15 ) * tDataMatrix( i, 4 );
                tTractionY1 = ( tDataMatrix( i, 17 ) - tDataMatrix( i, 13 ) ) * tDataMatrix( i, 4 ) + tDataMatrix( i, 15 ) * tDataMatrix( i, 3 );

                // getting x and y traction
                aPropMatrix( 0 ) = tTractionX0 + ( tTractionX1 - tTractionX0 ) * ( tX( 1 ) - tDataMatrix( i - 1, 2 ) ) / ( tDataMatrix( i, 2 ) - tDataMatrix( i - 1, 2 ) );

                aPropMatrix( 1 ) = tTractionY0 + ( tTractionY1 - tTractionY0 ) * ( tX( 1 ) - tDataMatrix( i - 1, 2 ) ) / ( tDataMatrix( i, 2 ) - tDataMatrix( i - 1, 2 ) );
                break;
            }
        }
        // aPropMatrix = 0.0 * aPropMatrix;
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Select_X(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 1.0;
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( Module_Parameter_Lists& aParameterLists )
    {

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "/projects/EXA/structure/thermo_elastic/MOPAR/exomesh.e";

        aParameterLists( 0 ).add_parameter_list( prm::create_stk_parameter_list() );;
        aParameterLists.set( "input_file", tMeshFileName );

        gLogger.set_severity_level( 0 );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Density" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity1" );
        aParameterLists.set( "function_parameters", std::to_string( tDens_c ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity2" );
        aParameterLists.set( "function_parameters", std::to_string( tDens_al ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Thermal_Modulus_Geo" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs1" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Thermal_Modulus_Block" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs2" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Thermal_Modulus_Block" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Poisson" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson1" );
        aParameterLists.set( "function_parameters", std::to_string( tPois_c ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson2" );
        aParameterLists.set( "function_parameters", std::to_string( tPois_al ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropCTE" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_CTE" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropCTE1" );
        aParameterLists.set( "function_parameters", std::to_string( tCTE_c ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropCTE2" );
        aParameterLists.set( "function_parameters", std::to_string( tCTE_al ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropRefTemp" );
        aParameterLists.set( "function_parameters", tRefTemp );

        // properties of dummy L2 projection of stress fields (supression for RBMs)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "DummySource" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "DummyL2Coefficient" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "DummyDiffusionCoefficient" );
        aParameterLists.set( "function_parameters", "0.00001" );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties of boundary conditions
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropSelectX" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Select_X" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropNeumann" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Import_Traction" );

        // create rotation axis property
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropRotAxis" );
        aParameterLists.set( "function_parameters", "0.0,0.0;1.0,0.0" );
        aParameterLists.set( "value_function", "Func_AxiRotation" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropField" );
        aParameterLists.set( "value_function", "Func_Field" );
        aParameterLists.set( "field_dependencies", "FIELD_1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", fem::Model_Type::AXISYMMETRIC );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );

        if ( byGeometry )
        {
            aParameterLists.set( "properties",
                    "PropYoungs,YoungsModulus;"
                    "PropPoisson,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        else
        {
            aParameterLists.set( "properties",
                    "PropYoungs1,YoungsModulus;"
                    "PropPoisson1,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE1,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso2" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", fem::Model_Type::AXISYMMETRIC );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        if ( byGeometry )
        {
            aParameterLists.set( "properties",
                    "PropYoungs,YoungsModulus;"
                    "PropPoisson,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        else
        {
            aParameterLists.set( "properties",
                    "PropYoungs2,YoungsModulus;"
                    "PropPoisson2,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE2,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs1,Material" );

        //------------------------------------------------------------------------------
        // Bulk Terms
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropRotAxis,Thickness" );
        aParameterLists.set( "mesh_set_names", "block_1" );

        // Bulk Terms
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU_2" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropRotAxis,Thickness" );
        aParameterLists.set( "mesh_set_names", "block_2" );

        // Dirichlet Boundary
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletU_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropRotAxis,Thickness" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "side_2" );

        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletU_2" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropRotAxis,Thickness" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "side_3" );

        // Neumann Boundary
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_properties", "PropNeumann,Traction;PropRotAxis,Thickness" );
        aParameterLists.set( "mesh_set_names", "side_1" );

        // Stress projection for block 1
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGStress_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_VON_MISES_STRESS );
        aParameterLists.set( "dof_residual", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", "block_1" );

        // Dummy equation for STRESS_DOF in block_2
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGL2_2" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::L2 );
        aParameterLists.set( "dof_residual", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_properties",    //
                "DummySource,Source;DummyL2Coefficient,L2coefficient;DummyDiffusionCoefficient,Diffusion" );
        aParameterLists.set( "mesh_set_names", "block_2" );

        // Stress projection for block 2
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGStress_2" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_VON_MISES_STRESS );
        aParameterLists.set( "dof_residual", "VX" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", "block_2" );

        // Dummy equation for VX in block_1
        aParameterLists( FEM::IWG ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGL2_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::L2 );
        aParameterLists.set( "dof_residual", "VX" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "leader_properties",    //
                "DummySource,Source;DummyL2Coefficient,L2coefficient;DummyDiffusionCoefficient,Diffusion" );
        aParameterLists.set( "mesh_set_names", "block_1" );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "block_1,block_2" );

        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", "block_1,block_2" );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIVMStress1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "STRESS_DOF" );
        aParameterLists.set( "leader_dof_dependencies", "STRESS_DOF" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "block_1" );

        // nodal von-mises stresses for shell
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIVMStress2" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "VX" );
        aParameterLists.set( "leader_dof_dependencies", "VX" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", "block_2" );

        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQITempField" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropField,Property" );
        aParameterLists.set( "mesh_set_names", "block_1,block_2" );

        // Emod
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIEmod1" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropYoungs1,Property" );
        aParameterLists.set( "mesh_set_names", "block_1" );

        // Emod
        aParameterLists( FEM::IQI ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIEmod2" );
        aParameterLists.set( "IQI_type", ( fem::IQI_Type::PROPERTY ) );
        aParameterLists.set( "leader_properties", "PropYoungs2,Property" );
        aParameterLists.set( "mesh_set_names", "block_2" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION ).add_parameter_list( prm::create_computation_parameter_list() );

        //------------------------------------------------------------------------------

        aParameterLists( FEM::FIELDS ).add_parameter_list( prm::create_fem_field_parameter_list() );
        aParameterLists.set( "field_name", "FieldNodalTEMP" );
        aParameterLists.set( "field_entity_type", "NODAL" );
        aParameterLists.set( "field_type", "FIELD_1" );
        aParameterLists.set( "field_create_from_file", tFieldRefPath );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-09 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 10.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "STRESS_DOF" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", "VX" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,2" );
        aParameterLists.set( "NLA_DofTypes", "UX,UY;STRESS_DOF;VX" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists.set( "TSA_Nonlinear_Solver", 3 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", "UX,UY;STRESS_DOF;VX" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", "block_1,block_2" );
        aParameterLists.set( "Field_Names", "UX,UY,STRESS1,STRESS2,Temp,EMOD1,EMOD2" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,ELEMENTAL_AVG,ELEMENTAL_AVG" );
        aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIVMStress1,IQIVMStress2,IQITempField,IQIEmod1,IQIEmod2" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
