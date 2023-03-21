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
    std::string tFieldRefPath   = tRoot + "projects/EXA/structure/thermo_elastic/MOPAR/exomesh.e";

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Density(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        MORIS_ASSERT( aParameters( 0 ).n_cols() == 2 and aParameters( 0 ).n_rows() == 2,
                "Axisymmetric rotation axis incorrectly defined. use {{x1,y1},{x2,y2}}." );

        // vector from point 1 to point 2
        moris::Matrix< moris::DDRMat > tRotVec = { { aParameters( 0 )( 1, 0 ) - aParameters( 0 )( 0, 0 ) },
            { aParameters( 0 )( 1, 1 ) - aParameters( 0 )( 0, 1 ) } };

        // spatial location of interest
        moris::Matrix< moris::DDRMat > tX = aFIManager->get_IP_geometry_interpolator()->valx();

        // vector from point 1 to location of interest
        moris::Matrix< moris::DDRMat > tPntVec = { { tX( 0 ) - aParameters( 0 )( 0, 0 ) }, { tX( 1 ) - aParameters( 0 )( 0, 1 ) } };

        // minimum distance vector from line to point of interest
        moris::Matrix< moris::DDRMat > tRadVec = tPntVec - dot( tPntVec, tRotVec ) * tRotVec / ( tRotVec( 0 ) * tRotVec( 0 ) + tRotVec( 1 ) * tRotVec( 1 ) );

        // init aPropMatrixSize
        aPropMatrix.set_size( 4, 1, 0.0 );

        // radius
        aPropMatrix( 1 ) = norm( tRadVec );

        // 2*pi*r for IWG
        aPropMatrix( 0 ) = aPropMatrix( 1 ) * 4.0 * std::acos( 0 );

        // outward radial normal vector
        aPropMatrix( { 2, 3 } ) = tRadVec / aPropMatrix( 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Field(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->val();
    }

    /* ------------------------------------------------------------------------ */
    void
    Func_Import_Traction(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // Input file
        std::string tBaseDir      = get_base_moris_dir();
        std::string tDataFilePath = tBaseDir + "projects/EXA/structure/thermo_elastic/MOPAR/wall_data.txt";

        // create ascii class
        Ascii tAscii = Ascii( tDataFilePath, FileMode::OPEN_RDONLY );

        // get num rows and columns of data file
        uint             tNumRows = tAscii.length();
        Matrix< DDRMat > tRow;
        string_to_mat( tAscii.line( 0 ), tRow );
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
            string_to_mat( tAscii.line( i ), tRow );

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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
    OPTParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
        tParameterlist( 0 )( 0 ).set( "workflow", "STK_FEM" );
    }

    void
    STKParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "/projects/EXA/structure/thermo_elastic/MOPAR/exomesh.e";

        tParameterlist( 0 )( 0 ) = prm::create_stk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "input_file", tMeshFileName );

        gLogger.set_severity_level( 0 );
    }

    //------------------------------------------------------------------------------

    void
    FEMParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Density" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tDens_c ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tDens_al ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Thermal_Modulus_Geo" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Thermal_Modulus_Block" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Thermal_Modulus_Block" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Poisson" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tPois_c ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tPois_al ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCTE" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_CTE" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCTE1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tCTE_c ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropCTE2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::to_string( tCTE_al ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropRefTemp" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tRefTemp );
        tPropCounter++;

        // properties of dummy L2 projection of stress fields (supression for RBMs)
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "DummySource" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "DummyL2Coefficient" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "DummyDiffusionCoefficient" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.00001" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // properties of boundary conditions
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichlet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropSelectX" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Select_X" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropNeumann" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Import_Traction" );
        tPropCounter++;

        // create rotation axis property
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropRotAxis" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0,0.0;1.0,0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_AxiRotation" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropField" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Field" );
        tParameterList( 0 )( tPropCounter ).set( "field_dependencies", "FIELD_1" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::AXISYMMETRIC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );

        if ( byGeometry )
        {
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropYoungs,YoungsModulus;"
                    "PropPoisson,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        else
        {
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropYoungs1,YoungsModulus;"
                    "PropPoisson1,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE1,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso2" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::AXISYMMETRIC ) );
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        if ( byGeometry )
        {
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropYoungs,YoungsModulus;"
                    "PropPoisson,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        else
        {
            tParameterList( 1 )( tCMCounter ).set( "properties",
                    "PropYoungs2,YoungsModulus;"
                    "PropPoisson2,PoissonRatio;"
                    "PropRotAxis,AxisymRotationAxis;"
                    "PropCTE2,CTE;"
                    "PropField,PropertyTemperature;"
                    "PropRefTemp,ReferenceTemperature" );
        }
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitsche" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
        tParameterList( 2 )( tSPCounter ).set( "master_properties", "PropYoungs1,Material" );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        // Bulk Terms
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropRotAxis,Thickness" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_1" );
        tIWGCounter++;

        // Bulk Terms
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropRotAxis,Thickness" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_2" );
        tIWGCounter++;

        // Dirichlet Boundary
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropRotAxis,Thickness" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "side_2" );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU_2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropDirichlet,Dirichlet;PropSelectX,Select;PropRotAxis,Thickness" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "side_3" );
        tIWGCounter++;

        // Neumann Boundary
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGNeumannFlux" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties", "PropNeumann,Traction;PropRotAxis,Thickness" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "side_1" );
        tIWGCounter++;

        // Stress projection for block 1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGStress_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_VON_MISES_STRESS ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "STRESS_DOF" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_1" );
        tIWGCounter++;

        // Dummy equation for STRESS_DOF in block_2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGL2_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::L2 );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "STRESS_DOF" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",    //
                "DummySource,Source;DummyL2Coefficient,L2coefficient;DummyDiffusionCoefficient,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_2" );
        tIWGCounter++;

        // Stress projection for block 2
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGStress_2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_VON_MISES_STRESS ) );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMStrucLinIso2,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_2" );
        tIWGCounter++;

        // Dummy equation for VX in block_1
        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGL2_1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type", (uint)fem::IWG_Type::L2 );
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies", "UX,UY;STRESS_DOF;VX" );
        tParameterList( 3 )( tIWGCounter ).set( "master_properties",    //
                "DummySource,Source;DummyL2Coefficient,L2coefficient;DummyDiffusionCoefficient,Diffusion" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", "block_1" );
        tIWGCounter++;

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_1,block_2" );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "UX,UY" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_1,block_2" );
        tIQICounter++;

        // nodal von-mises stresses for shell
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIVMStress1" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "STRESS_DOF" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_1" );
        tIQICounter++;

        // nodal von-mises stresses for shell
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIVMStress2" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::DOF ) );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "VX" );
        tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies", "VX" );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_2" );
        tIQICounter++;

        // Nodal Temperature IQI
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQITempField" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropField,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_1,block_2" );
        tIQICounter++;

        // Emod
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIEmod1" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropYoungs1,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_1" );
        tIQICounter++;

        // Emod
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIEmod2" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::PROPERTY ) );
        tParameterList( 4 )( tIQICounter ).set( "master_properties", "PropYoungs2,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", "block_2" );
        tIQICounter++;


        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tFieldCounter = 0;

        tParameterList( 6 ).push_back( prm::create_fem_field_parameter_list() );
        tParameterList( 6 )( tFieldCounter ).set( "field_name", "FieldNodalTEMP" );
        tParameterList( 6 )( tFieldCounter ).set( "field_entity_type", "NODAL" );
        tParameterList( 6 )( tFieldCounter ).set( "field_type", "FIELD_1" );
        tParameterList( 6 )( tFieldCounter ).set( "field_create_from_file", tFieldRefPath );
    }

    //------------------------------------------------------------------------------

    void
    SOLParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 7 );

        tParameterlist( 0 ).push_back( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-09 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.0 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 1 );

        tParameterlist( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 2 )( 1 ).set( "NLA_rel_res_norm_drop", 10.0 );
        tParameterlist( 2 )( 1 ).set( "NLA_max_iter", 1 );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "STRESS_DOF" );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "VX" );

        tParameterlist( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 )( 3 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
        tParameterlist( 3 )( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
        tParameterlist( 3 )( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2" );
        tParameterlist( 3 )( 3 ).set( "NLA_DofTypes", "UX,UY;STRESS_DOF;VX" );

        tParameterlist( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_solver", 3 );

        tParameterlist( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY;STRESS_DOF;VX" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );
    }

    //------------------------------------------------------------------------------

    void
    MSIParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
    }

    //------------------------------------------------------------------------------

    void
    VISParameterList( moris::Cell< moris::Cell< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", "block_1,block_2" );
        tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY,STRESS1,STRESS2,Temp,EMOD1,EMOD2" );
        tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL,ELEMENTAL,ELEMENTAL" );
        tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQIVMStress1,IQIVMStress2,IQITempField,IQIEmod1,IQIEmod2" );
        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
