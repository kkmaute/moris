/*
 * cl_Equation_Object_Pdv.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#define protected public
#define private   public
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Interpolation_Element.hpp"
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src

#undef protected
#undef private

#include "MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Vertex_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cell_Proxy.hpp"
#include "MSI_Test_Proxy/cl_MTK_Cluster_Proxy.hpp"

#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Geometry_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src


namespace moris
{
    namespace MSI
    {
    moris::Matrix< moris::DDRMat > tConstValFunction_FDTest( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                    moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                                    moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                                    moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
    {
        return aParameters( 0 );
    }

    TEST_CASE("Eqn_Obj_pdv","[MSI],[Eqn_Obj_pdv]")
    {
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 2;

        // Create generic Node Object
        fem::Node_Base * Node1 = new Node_Proxy( 0 );
        fem::Node_Base * Node2 = new Node_Proxy( 1 );
        fem::Node_Base * Node3 = new Node_Proxy( 2 );
        fem::Node_Base * Node4 = new Node_Proxy( 3 );
        //---------------------------------------------------------------------------------

        moris::uint tNumNodes = 4;

        // Create List with node pointer corresponding to generic equation object
        moris::Cell< moris::Cell< fem::Node_Base * > > tNodeIds_1( 1 );
        tNodeIds_1( 0 ).resize( tNumNodes );
        tNodeIds_1( 0 )( 0 ) = Node1;
        tNodeIds_1( 0 )( 1 ) = Node2;
        tNodeIds_1( 0 )( 2 ) = Node3;
        tNodeIds_1( 0 )( 3 ) = Node4;

        mtk::Vertex * tVertex_0 = new mtk::Vertex_Proxy( 0 ) ;
        mtk::Vertex * tVertex_1 = new mtk::Vertex_Proxy( 1 ) ;
        mtk::Vertex * tVertex_2 = new mtk::Vertex_Proxy( 2 ) ;
        mtk::Vertex * tVertex_3 = new mtk::Vertex_Proxy( 3 ) ;
        mtk::Vertex * tVertex_4 = new mtk::Vertex_Proxy( 4 ) ;
        mtk::Vertex * tVertex_5 = new mtk::Vertex_Proxy( 5 ) ;

        mtk::Cell * tCell_1 = new mtk::Cell_Proxy( 0, { tVertex_0, tVertex_4, tVertex_5 });
        mtk::Cell * tCell_2 = new mtk::Cell_Proxy( 1, { tVertex_4, tVertex_1, tVertex_2 });
        mtk::Cell * tCell_3 = new mtk::Cell_Proxy( 2, { tVertex_4, tVertex_2, tVertex_5 });
        mtk::Cell * tCell_4 = new mtk::Cell_Proxy( 3, { tVertex_5, tVertex_2, tVertex_3 });

        Matrix< DDRMat > tLocalCoords( 6, 2 );
        tLocalCoords( 0, 0 ) = -1;  tLocalCoords( 0, 1 ) = -1;
        tLocalCoords( 1, 0 ) =  1;  tLocalCoords( 1, 1 ) = -1;
        tLocalCoords( 2, 0 ) =  1;  tLocalCoords( 2, 1 ) =  1;
        tLocalCoords( 3, 0 ) = -1;  tLocalCoords( 3, 1 ) =  1;
        tLocalCoords( 4, 0 ) =  0;  tLocalCoords( 4, 1 ) = -1;
        tLocalCoords( 5, 0 ) = -1;  tLocalCoords( 5, 1 ) =  0;

        mtk::Cluster * tCluster = new mtk::Cluster_Proxy( { tCell_1 }, { tCell_2, tCell_3, tCell_4 }, tLocalCoords );

        Design_Variable_Interface * tDesignVariableInterface = new Design_Variable_Interface_Proxy();

        MSI::Equation_Set * tSet = new fem::Set();

        tSet->set_Dv_interface( tDesignVariableInterface );

        // Create generic equation objects
        MSI::Equation_Object * EquObj = new fem::Interpolation_Element();

        EquObj->mEquationBlock = tSet;
        reinterpret_cast< fem::Interpolation_Element *> ( EquObj )->mSet = reinterpret_cast< fem::Set *> (tSet);

        tSet->mMasterDofTypes = { { MSI::Dof_Type::TEMP } };


        std::shared_ptr< fem::Cluster > tFemCluster = std::make_shared< fem::Cluster >( fem::Element_Type::BULK,
                                                                                        tCluster,
                                                                                        reinterpret_cast< fem::Set *> (tSet),
                                                                                        EquObj );

        reinterpret_cast< fem::Interpolation_Element * >( EquObj )->mFemCluster = {tFemCluster};


        // create the properties
        std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
        tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
        tPropMasterConductivity->set_val_function( tConstValFunction_FDTest );

        std::shared_ptr< fem::Property > tPropMasterTempLoad = std::make_shared< fem::Property > ();
        tPropMasterTempLoad->set_parameters( { {{ 1.0 }} } );
        tPropMasterTempLoad->set_val_function( tConstValFunction_FDTest );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
        tCMMasterDiffLinIso->set_space_dim( 2 );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
        tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
        tIWG->set_constitutive_model( tCMMasterDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
        tIWG->set_property( tPropMasterTempLoad, "Load", mtk::Master_Slave::MASTER );

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        fem::Interpolation_Rule tIPGIRule( mtk::Geometry_Type::QUAD,
                                      fem::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      fem::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        fem::Geometry_Interpolator tIPGI( tIPGIRule );

        // create a space geometry interpolation rule
        fem::Interpolation_Rule tIGGIRule( mtk::Geometry_Type::TRI,
                                           fem::Interpolation_Type::LAGRANGE,
                                           mtk::Interpolation_Order::LINEAR,
                                           fem::Interpolation_Type::LAGRANGE,
                                           mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        fem::Geometry_Interpolator tIGGI( tIGGIRule );

        // field interpolators
        //------------------------------------------------------------------------------
        //create a space time interpolation rule
        fem::Interpolation_Rule tFIRule ( mtk::Geometry_Type::QUAD,
                                          fem::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR,
                                          fem::Interpolation_Type::CONSTANT,
                                          mtk::Interpolation_Order::CONSTANT );

        // create a cell of field interpolators for IWG
        Cell< fem::Field_Interpolator* > tFIs( 1 );

        // create the field interpolator
        tFIs( 0 ) = new fem::Field_Interpolator( 1, tFIRule, &tIPGI, { MSI::Dof_Type::TEMP } );

        // set a fem set pointer
        tIWG->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
        std::cout<<tIWG->mSet<<" 1"<<std::endl;

        // set size for the set EqnObjDofTypeList
        tIWG->mSet->mEqnObjDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

        // set master dof type list
        tIWG->mSet->mMasterDofTypes = { { MSI::Dof_Type::TEMP } };

        // set size and populate the set dof type map
        tIWG->mSet->mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tIWG->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

        // set size and populate the set master dof type map
        tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
        tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

        // set size and fill the set residual assembly map
        tIWG->mSet->mResDofAssemblyMap.resize( 1 );
        tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 7 } };

        // set size and init the set residual and jacobian
        tIWG->mSet->mResidual.set_size( 8, 1, 0.0 );

        // set requested residual dof type flag to true
        tIWG->mResidualDofTypeRequested = true;

        // build global dof type list
        tIWG->get_global_dof_type_list();

        // populate the requested master dof type
        tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

        // create a field interpolator manager
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
        fem::Field_Interpolator_Manager tFIManager( tDummy, tSet );

        // populate the field interpolator manager
        tFIManager.mFI = tFIs;
        tFIManager.mIPGeometryInterpolator = &tIPGI;
        tFIManager.mIGGeometryInterpolator = &tIGGI;

        tIWG->mSet->mMasterFIManager = &tFIManager;

        // set IWG field interpolator manager
        tIWG->set_field_interpolator_manager( &tFIManager );

        EquObj->mPdofValues = { {0, 0, 0, 0} };

        EquObj->compute_dRdp();



//        // Create the pdof hosts of this equation object
//        moris::Cell < Pdof_Host * > tPdofHostList;
//        tPdofHostList.resize( 3, nullptr );
//        moris::uint tNumMaxPdofTypes = 1;
//
//        Matrix< DDSMat > tDofTypeIndexMap(4, 1, -1);
//        tDofTypeIndexMap(3, 0) = 0;
//
//        Matrix< DDUMat > tTimePerDofType(4, 1, 1);
//
//        Equation_Set tEqnBlock;
//        tEqnBlock.mEqnObjDofTypeList.resize( 1, MSI::Dof_Type::TEMP );
//        EquObj.mEquationBlock = &tEqnBlock;
//
//        EquObj.create_my_pdof_hosts( tNumMaxPdofTypes, tDofTypeIndexMap, tTimePerDofType, tPdofHostList );
//
//        // Check if right pdof host was created in given pdof host list
//        CHECK( equal_to( tPdofHostList( 0 )->mNodeID, 0 ) );
//        REQUIRE( tPdofHostList( 1 ) == NULL );
//        CHECK( equal_to( tPdofHostList( 2 )->mNodeID, 2 ) );
//
//        // Check equation objects internal pdof host list
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 ).size(), 2 ) );
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 )( 0 )->mNodeID, 0 ) );
//        CHECK( equal_to( EquObj.mMyPdofHosts( 0 )( 1 )->mNodeID, 2 ) );
        delete Node1;
        delete Node2;
//        delete tPdofHostList(0);
//        delete tPdofHostList(1);
//        delete tPdofHostList(2);
    }

    }
}


