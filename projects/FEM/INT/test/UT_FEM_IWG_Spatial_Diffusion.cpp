#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp" //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp" //FEM//INT//src

#include "op_equal_equal.hpp"


moris::Matrix< moris::DDRMat > tConstValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                  moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                  moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tGeoValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aGeometryInterpolator->valx();
}

moris::Matrix< moris::DDRMat > tFIValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aFieldInterpolator( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                               moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aFieldInterpolator( 0 )->N();;
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_SpatialDiff", "[moris],[fem],[IWG_SpatialDiff]" )
{

    // create a spatial diffusion bulk IWG
    //------------------------------------------------------------------------------

    // create an IWG Spatial Difffusion Bulk
    IWG_Isotropic_Spatial_Diffusion_Bulk tIWG;

    // set space dimension
    tIWG.set_space_dim( 3 );

    // set residual dof type
    tIWG.set_residual_dof_type( { MSI::Dof_Type::TEMP } );

    // set active dof type
    tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

    // set active property type
    tIWG.set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator* tGI = new Geometry_Interpolator( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                              { 1.0, 0.0, 0.0 },
                              { 1.0, 1.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 },
                              { 1.0, 0.0, 1.0 },
                              { 1.0, 1.0, 1.0 },
                              { 0.0, 1.0, 1.0 }};

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

    // set the coefficients xHat, tHat
    tGI->set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI->set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tFIRule ( mtk::Geometry_Type::HEX,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR,
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create coefficients
    Matrix< DDRMat > tDOFHat( 8, 1 );
    tDOFHat = {{1.0},{1.0},{1.0},{1.0},{2.0},{2.0},{2.0},{2.0}};

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( tIWG.get_dof_type_list().size() );

    for( uint iDOF = 0; iDOF < tIWG.get_dof_type_list().size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG.get_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               tGI,
                                               tIWG.get_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // define an epsilon environment
    double tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with constant property" )
    {
        // properties
        //------------------------------------------------------------------------------
        // create property coefficients
        Cell< Matrix< DDRMat > > tPropCoeff = { {{1.0}} };

        // create a cell of properties for IWG
        Cell< Property* > tProps( tIWG.get_property_type_list().size() );

        for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
        {
            // create a property
            tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                            Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                            tPropCoeff,
                                            tConstValFunction,
                                            Cell< PropertyFunc > ( 0 ),
                                            tGI );
        }
        // set IWG properties
        tIWG.set_properties( tProps );

        // set IWG field interpolators
        tIWG.set_field_interpolators( tFIs );

        // check evaluation of the residual for IWG Helmholtz Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual
        Cell< Matrix< DDRMat > > tResidual;
        tIWG.compute_residual( tResidual );

        // check evaluation of the jacobian  by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian
        Cell< Cell< Matrix< DDRMat > > > tJacobians;
        tIWG.compute_jacobian( tJacobians );
        //print( tJacobians( 0 )( 0 ),"tJacobians");

        Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
        tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
        //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

        //define a boolean for check
        bool tCheckJacobian = true;

        for ( uint iJac = 0; iJac < tJacobians.size(); iJac++ )
        {
            for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
            {
                for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                {
                    for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                    {
                        tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                    }
                }
            }
        }

        REQUIRE( tCheckJacobian );

        for( Property* tProp : tProps )
        {
            delete tProp;
        }
        tProps.clear();

    }/* END_SECTION */

    SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with property dependent on x" )
        {
            // properties
            //------------------------------------------------------------------------------
            // create property coefficients
            Cell< Matrix< DDRMat > > tPropCoeff = { {{1.0}} };

            // create a cell of properties for IWG
            Cell< Property* > tProps( tIWG.get_property_type_list().size() );

            for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
            {
                // create a property
                tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                                Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                                tPropCoeff,
                                                tGeoValFunction,
                                                Cell< PropertyFunc > ( 0 ),
                                                tGI );
            }
            // set IWG properties
            tIWG.set_properties( tProps );

            // set IWG field interpolators
            tIWG.set_field_interpolators( tFIs );

            // check evaluation of the residual for IWG Helmholtz Bulk ?
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG.compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // evaluate the jacobian
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            tIWG.compute_jacobian( tJacobians );
            //print( tJacobians( 0 )( 0 ),"tJacobians");

            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
            tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
            //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

            //define a boolean for check
            bool tCheckJacobian = true;

            for ( uint iJac = 0; iJac 									< tJacobians.size(); iJac++ )
            {
                for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
                {
                    for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                        }
                    }
                }
            }

            REQUIRE( tCheckJacobian );

            for( Property* tProp : tProps )
            {
                delete tProp;
            }
            tProps.clear();

        }/* END_SECTION */

        SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with property dependent on TEMP" )
        {
            // properties
            //------------------------------------------------------------------------------
            // create property coefficients
            Cell< Matrix< DDRMat > > tPropCoeff = { {{1.0}} };

            // create a cell of properties for IWG
            Cell< Property* > tProps( tIWG.get_property_type_list().size() );

            for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
            {
                // create a property
                tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                                {{ MSI::Dof_Type::TEMP }},
                                                tPropCoeff,
                                                tFIValFunction,
                                                { tFIDerFunction },
                                                tGI );

                // set field interpolators
                tProps( iProp )->set_field_interpolators( tFIs );
            }

            // set IWG properties
            tIWG.set_properties( tProps );

            // set IWG field interpolators
            tIWG.set_field_interpolators( tFIs );

            // check evaluation of the residual for IWG Helmholtz Bulk ?
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG.compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // evaluate the jacobian
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            tIWG.compute_jacobian( tJacobians );
            //print( tJacobians( 0 )( 0 ),"tJacobians");

            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
            tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
            //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

            //define a boolean for check
            bool tCheckJacobian = true;

            for ( uint iJac = 0; iJac < tJacobians.size(); iJac++ )
            {
                for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
                {
                    for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                        }
                    }
                }
            }

            REQUIRE( tCheckJacobian );

            for( Property* tProp : tProps )
            {
                delete tProp;
            }
            tProps.clear();

        }/* END_SECTION */

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

    delete tGI;

}/* END_TEST_CASE */


TEST_CASE( "IWG_SpatialDiff_Dirichlet", "[moris],[fem],[IWG_SpatialDiff_Dirichlet]" )
{

    // create a spatial diffusion bulk IWG
    //------------------------------------------------------------------------------

    // create an IWG Spatial Difffusion Bulk
    IWG_Isotropic_Spatial_Diffusion_Dirichlet tIWG;

    // set space dimension
    tIWG.set_space_dim( 3 );

    // set residual dof type
    tIWG.set_residual_dof_type( { MSI::Dof_Type::TEMP } );

    // set active dof type
    tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

    // set active property type
    tIWG.set_property_type_list( { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET } );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0}, {-0.25}, { 0.75}, { 0.0 }};

    // set the normal
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG.set_normal( tNormal );

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator* tGI = new Geometry_Interpolator( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                              { 1.0, 0.0, 0.0 },
                              { 1.0, 1.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 },
                              { 1.0, 0.0, 1.0 },
                              { 1.0, 1.0, 1.0 },
                              { 0.0, 1.0, 1.0 }};

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

    // set the coefficients xHat, tHat
    tGI->set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI->set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tFIRule ( mtk::Geometry_Type::HEX,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR,
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create coefficients
    Matrix< DDRMat > tDOFHat( 8, 1 );
    tDOFHat = {{1.0},{1.0},{1.0},{1.0},{2.0},{2.0},{2.0},{2.0}};

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( tIWG.get_dof_type_list().size() );

    for( uint iDOF = 0; iDOF < tIWG.get_dof_type_list().size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG.get_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               tGI,
                                               tIWG.get_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // define an epsilon environment
    double tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with constant property" )
    {
        // properties
        //------------------------------------------------------------------------------
        // create property coefficients
        Cell< Cell< Matrix< DDRMat > > > tPropCoeff( 2 );
        tPropCoeff( 0 ) = { {{1.0}} };
        tPropCoeff( 1 ) = { {{1.0}} };

        // create a cell of properties for IWG
        Cell< Property* > tProps( tIWG.get_property_type_list().size() );

        for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
        {
            // create a property
            tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                            Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                            tPropCoeff( iProp ),
                                            tConstValFunction,
                                            Cell< PropertyFunc > ( 0 ),
                                            tGI );
        }
        // set IWG properties
        tIWG.set_properties( tProps );

        // set IWG field interpolators
        tIWG.set_field_interpolators( tFIs );

        // check evaluation of the residual for IWG Helmholtz Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual
        Cell< Matrix< DDRMat > > tResidual;
        tIWG.compute_residual( tResidual );

        // check evaluation of the jacobian  by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian
        Cell< Cell< Matrix< DDRMat > > > tJacobians;
        tIWG.compute_jacobian( tJacobians );
        //print( tJacobians( 0 )( 0 ),"tJacobians");

        Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
        tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
        //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

        //define a boolean for check
        bool tCheckJacobian = true;

        for ( uint iJac = 0; iJac < tJacobians.size(); iJac++ )
        {
            for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
            {
                for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                {
                    for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                    {
                        tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                    }
                }
            }
        }

        REQUIRE( tCheckJacobian );

        for( Property* tProp : tProps )
        {
            delete tProp;
        }
        tProps.clear();

    }/* END_SECTION */

    SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with property dependent on x" )
        {
            // properties
            //------------------------------------------------------------------------------
            // create property coefficients
            Cell< Cell< Matrix< DDRMat > > > tPropCoeff( 2 );
            tPropCoeff( 0 ) = { {{1.0}} };
            tPropCoeff( 1 ) = { {{1.0}} };

            // create a cell of properties for IWG
            Cell< Property* > tProps( tIWG.get_property_type_list().size() );

            for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
            {
                // create a property
                tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                                Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                                tPropCoeff( 0 ),
                                                tGeoValFunction,
                                                Cell< PropertyFunc > ( 0 ),
                                                tGI );

            }
            // set IWG properties
            tIWG.set_properties( tProps );

            // set IWG field interpolators
            tIWG.set_field_interpolators( tFIs );

            // check evaluation of the residual for IWG Helmholtz Bulk ?
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG.compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // evaluate the jacobian
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            tIWG.compute_jacobian( tJacobians );
            //print( tJacobians( 0 )( 0 ),"tJacobians");

            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
            tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
            //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

            //define a boolean for check
            bool tCheckJacobian = true;

            for ( uint iJac = 0; iJac < tJacobians.size(); iJac++ )
            {
                for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
                {
                    for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                        }
                    }
                }
            }

            REQUIRE( tCheckJacobian );

            for( Property* tProp : tProps )
            {
                delete tProp;
            }
            tProps.clear();

        }/* END_SECTION */

        SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with property dependent on TEMP" )
        {
            // properties
            //------------------------------------------------------------------------------
            // create property coefficients
            Cell< Cell< Matrix< DDRMat > > > tPropCoeff( 2 );
            tPropCoeff( 0 ) = { {{1.0}} };
            tPropCoeff( 1 ) = { {{1.0}} };

            // create a cell of properties for IWG
            Cell< Property* > tProps( tIWG.get_property_type_list().size() );

            for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
            {
                // create a property
                tProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                                {{ MSI::Dof_Type::TEMP }},
                                                tPropCoeff( 1 ),
                                                tFIValFunction,
                                                { tFIDerFunction },
                                                tGI );

                // set field interpolators
                tProps( iProp )->set_field_interpolators( tFIs );
            }

            // set IWG properties
            tIWG.set_properties( tProps );

            // set IWG field interpolators
            tIWG.set_field_interpolators( tFIs );

            // check evaluation of the residual for IWG Helmholtz Bulk ?
            //------------------------------------------------------------------------------
            // evaluate the residual
            Cell< Matrix< DDRMat > > tResidual;
            tIWG.compute_residual( tResidual );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // evaluate the jacobian
            Cell< Cell< Matrix< DDRMat > > > tJacobians;
            tIWG.compute_jacobian( tJacobians );
            //print( tJacobians( 0 )( 0 ),"tJacobians");

            Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
            tIWG.compute_jacobian_FD( tJacobiansFD, tPerturbation );
            //print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD");

            //define a boolean for check
            bool tCheckJacobian = true;

            for ( uint iJac = 0; iJac < tJacobians.size(); iJac++ )
            {
                for( uint jJac = 0; jJac < tJacobians( iJac ).size(); jJac++ )
                {
                    for( uint iiJac = 0; iiJac < tJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < tJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( tJacobians( iJac )( jJac )( iiJac, jjJac ) - tJacobiansFD( iJac )( jJac )( iiJac, jjJac ) < tEpsilon );
                        }
                    }
                }
            }

            REQUIRE( tCheckJacobian );

            for( Property* tProp : tProps )
            {
                delete tProp;
            }
            tProps.clear();

        }/* END_SECTION */

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

    delete tGI;

}/* END_TEST_CASE */
