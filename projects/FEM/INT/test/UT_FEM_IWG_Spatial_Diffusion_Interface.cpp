#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Interface.hpp" //FEM//INT//src

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

TEST_CASE( "IWG_SpatialDiff_Interface", "[moris],[fem],[IWG_SpatialDiff_Interface]" )
{

    // create a spatial diffusion bulk IWG
    //------------------------------------------------------------------------------

    // create an IWG Spatial Difffusion Bulk
    IWG_Isotropic_Spatial_Diffusion_Interface tIWG;

    // set space dimension
    tIWG.set_space_dim( 3 );

    // set residual dof type
    tIWG.set_residual_dof_type( { MSI::Dof_Type::TEMP } );

    // set master dof type
    tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }});

    // set slave dof type
    tIWG.set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );

    // set active constitutive type
    tIWG.set_constitutive_type_list( { fem::Constitutive_Type::DIFF_LIN_ISO } );

    // set active constitutive type
    tIWG.set_constitutive_type_list( { fem::Constitutive_Type::DIFF_LIN_ISO }, mtk::Master_Slave::SLAVE );

    // set master property type
    tIWG.set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

    // set slave property type
    tIWG.set_property_type_list( { fem::Property_Type::CONDUCTIVITY }, mtk::Master_Slave::SLAVE );

    // set the normal
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG.set_normal( tNormal );

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
    Cell< Field_Interpolator* > tMasterFIs( tIWG.get_dof_type_list().size() );

    for( uint iDOF = 0; iDOF < tIWG.get_dof_type_list().size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG.get_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tMasterFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                                     tFIRule,
                                                     tGI,
                                                     tIWG.get_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tMasterFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tMasterFIs( iDOF )->set_space_time( tParamPoint );
    }

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tSlaveFIs( tIWG.get_dof_type_list( mtk::Master_Slave::SLAVE ).size() );

    for( uint iDOF = 0; iDOF < tIWG.get_dof_type_list( mtk::Master_Slave::SLAVE ).size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG.get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ).size();

        // create the field interpolator
        tSlaveFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                                    tFIRule,
                                                    tGI,
                                                    tIWG.get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ) );

        // set the coefficients uHat
        tSlaveFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tSlaveFIs( iDOF )->set_space_time( tParamPoint );
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
        Cell< Property* > tMasterProps( tIWG.get_property_type_list().size() );

        for( uint iProp = 0; iProp < tIWG.get_property_type_list().size(); iProp++ )
        {
            // create a property
            tMasterProps( iProp ) = new Property( tIWG.get_property_type_list()( iProp ),
                                            Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                            tPropCoeff,
                                            tConstValFunction,
                                            Cell< PropertyFunc > ( 0 ),
                                            tGI );
        }

        // create a cell of properties for IWG
        Cell< Property* > tSlaveProps( tIWG.get_property_type_list( mtk::Master_Slave::SLAVE ).size() );

        for( uint iProp = 0; iProp < tIWG.get_property_type_list( mtk::Master_Slave::SLAVE ).size(); iProp++ )
        {
            // create a property
        	tSlaveProps( iProp ) = new Property( tIWG.get_property_type_list( mtk::Master_Slave::SLAVE )( iProp ),
                                            Cell< Cell< MSI::Dof_Type > > ( 0 ),
                                            tPropCoeff,
                                            tConstValFunction,
                                            Cell< PropertyFunc > ( 0 ),
                                            tGI );
        }

        // constitutive models
        //------------------------------------------------------------------------------
        // create a cell of properties for IWG
        Cell< Constitutive_Model* > tMasterCMs( tIWG.get_constitutive_type_list().size() );

        // create a constitutive model factory
        fem::CM_Factory tCMFactory;

        // create a constitutive model for each constitutive type
        for( uint iCM = 0; iCM < tIWG.get_constitutive_type_list().size(); iCM++ )
        {
            // create a property
            tMasterCMs( iCM ) = tCMFactory.create_CM( tIWG.get_constitutive_type_list()( iCM ) );

            // set space dim
            tMasterCMs( iCM )->set_space_dim( 3 );

            // set dof types
            tMasterCMs( iCM )->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set property type
            tMasterCMs( iCM )->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

            // set properties
            tMasterCMs( iCM )->set_properties( tMasterProps );

            // set field interpolators
            tMasterCMs( iCM )->set_field_interpolators( tMasterFIs );
        }

        // create a cell of properties for IWG
        Cell< Constitutive_Model* > tSlaveCMs( tIWG.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size() );

        // create a constitutive model for each constitutive type
        for( uint iCM = 0; iCM < tIWG.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size(); iCM++ )
        {
            // create a property
            tSlaveCMs( iCM ) = tCMFactory.create_CM( tIWG.get_constitutive_type_list( mtk::Master_Slave::SLAVE )( iCM ) );

            // set space dim
            tSlaveCMs( iCM )->set_space_dim( 3 );

            // set dof types
            tSlaveCMs( iCM )->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

            // set property type
            tSlaveCMs( iCM )->set_property_type_list( { fem::Property_Type::CONDUCTIVITY } );

            // set properties
            tSlaveCMs( iCM )->set_properties( tSlaveProps );

            // set field interpolators
            tSlaveCMs( iCM )->set_field_interpolators( tSlaveFIs );
        }

        // set IWG field interpolators
        tIWG.set_constitutive_models( tMasterCMs );
        tIWG.set_constitutive_models( tSlaveCMs, mtk::Master_Slave::SLAVE );

        // set IWG properties
        tIWG.set_properties( tMasterProps );
        tIWG.set_properties( tSlaveProps, mtk::Master_Slave::SLAVE );

        // set IWG field interpolators
        tIWG.set_field_interpolators( tMasterFIs );
        tIWG.set_field_interpolators( tSlaveFIs, mtk::Master_Slave::SLAVE );

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

        Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
        tIWG.compute_jacobian_FD_double( tJacobiansFD, tPerturbation );

//        print( tJacobians( 0 )( 0 ),"tJacobians00");
//        print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD00");
//
//        print( tJacobians( 0 )( 1 ),"tJacobians01");
//        print( tJacobiansFD( 0 )( 1 ),"tJacobiansFD01");
//
//        print( tJacobians( 1 )( 0 ),"tJacobians10");
//        print( tJacobiansFD( 1 )( 0 ),"tJacobiansFD10");
//
//        print( tJacobians( 1 )( 1 ),"tJacobians11");
//        print( tJacobiansFD( 1 )( 1 ),"tJacobiansFD11");

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

        for( Property* tProp : tMasterProps )
        {
            delete tProp;
        }
        tMasterProps.clear();

        for( Property* tProp : tSlaveProps )
        {
            delete tProp;
        }
        tSlaveProps.clear();

        for( Constitutive_Model* tCM : tMasterCMs )
        {
            delete tCM;
        }
        tMasterCMs.clear();

        for( Constitutive_Model* tCM : tSlaveCMs )
        {
            delete tCM;
        }
        tSlaveCMs.clear();

    }/* END_SECTION */

    // clean up
    for( Field_Interpolator* tFI : tMasterFIs )
    {
        delete tFI;
    }
    tMasterFIs.clear();

    for( Field_Interpolator* tFI : tSlaveFIs )
    {
        delete tFI;
    }
    tSlaveFIs.clear();

    delete tGI;

}/* END_TEST_CASE */
