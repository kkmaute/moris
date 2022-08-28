/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Set.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "linalg_typedefs.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src

#define protected public
#define private   public
#include "cl_FEM_Set.hpp"                //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp" //FEM/INT/src
#undef protected
#undef private

namespace moris
{
    namespace fem
    {

void tConstValFunction_UTFEMSET
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

        // This test check all the member functions of the FEM Set in the particular case where
        // there is only a master.
        TEST_CASE( "Set", "[moris],[fem],[FEMSet_Master]" )
        {
//            //create a set
//            Set tSet;
//
//            // set spatial dimension
//            tSet.mSpaceDim = 3;
//
//            // list of IWG types
//            moris::Cell< fem::IWG_Type >  tIWGTypeList= { fem::IWG_Type::SPATIALDIFF_BULK ,
//                                                          fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE,
//                                                          fem::IWG_Type::HELMHOLTZ,
//                                                          fem::IWG_Type::LSNORMAL };
//
//            // number of IWGs to be created
//            uint tNumOfIWGs = tIWGTypeList.size();
//
//            // list of residual dof type
//            moris::Cell< moris::Cell< MSI::Dof_Type > > aResidualDofType( tNumOfIWGs );
//            aResidualDofType( 0 ) = { MSI::Dof_Type::TEMP };
//            aResidualDofType( 1 ) = { MSI::Dof_Type::TEMP };
//            aResidualDofType( 2 ) = { MSI::Dof_Type::VX };
//            aResidualDofType( 3 ) = { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ };
//
//            // list of active dof type
//            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > aMasterDofTypes( tNumOfIWGs );
//            aMasterDofTypes( 0 ) = {{ MSI::Dof_Type::TEMP }};
//            aMasterDofTypes( 1 ) = {{ MSI::Dof_Type::TEMP }};
//            aMasterDofTypes( 2 ) = {{ MSI::Dof_Type::VX }};
//            aMasterDofTypes( 3 ) = {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ },
//                                    { MSI::Dof_Type::LS1}};
//
//            // list of active property type
//            moris::Cell< moris::Cell< fem::Property_Type > > aMasterPropTypes( tNumOfIWGs );
//            aMasterPropTypes( 0 ) = { fem::Property_Type::CONDUCTIVITY };
//            aMasterPropTypes( 1 ) = { fem::Property_Type::CONDUCTIVITY,
//                                      fem::Property_Type::TEMP_DIRICHLET };
//
//            // list of active constitutive type
//            moris::Cell< moris::Cell< fem::Constitutive_Type > > aMasterConstitutiveTypes( tNumOfIWGs );
//            aMasterConstitutiveTypes( 0 ) = { fem::Constitutive_Type::DIFF_LIN_ISO };
//            aMasterConstitutiveTypes( 1 ) = { fem::Constitutive_Type::DIFF_LIN_ISO  };
//
//            // a factory to create the IWGs
//            fem::IWG_Factory tIWGFactory;
//
//            // create a cell of IWGs for the problem considered
//            moris::Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );
//
//           // loop over the IWG types
//           for( uint i = 0; i < tNumOfIWGs; i++)
//           {
//               // create an IWG with the factory for the ith IWG type
//               tIWGs( i ) = tIWGFactory.create_IWG( tIWGTypeList( i ) );
//
//               // set residual dof type
//               tIWGs( i )->set_residual_dof_type( aResidualDofType( i ) );
//
//               // set active dof type
//               tIWGs( i )->set_dof_type_list( aMasterDofTypes( i ) );
//
//               // set active property type
//               tIWGs( i )->set_property_type_list( aMasterPropTypes( i ) );
//
//               // set active constitutive type
//               tIWGs( i )->set_constitutive_type_list( aMasterConstitutiveTypes( i ) );
//           }
//
//           // pass in the cell of IWG pointers to the element block
//            tSet.mIWGs = tIWGs;
//
//            //std::cout<<"Test create_constitutive_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create unique list of constitutive type
//                tSet.create_constitutive_type_list();
//
//                // check mMasterConstitutiveTypes size
//                CHECK( equal_to( tSet.mMasterConstitutiveTypes.size(), 1 ) );
//
//                // check mMasterConstitutiveTypes content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterConstitutiveTypes( 0 ) ), 1 ) ); // DIFF_LIN_ISO
//
//            //std::cout<<"Test create_constitutive_type_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a constitutive type map for the set
//                tSet.create_constitutive_type_map();
//
//                // check mMasterConstitutiveTypeMap size
//                CHECK( equal_to( tSet.mMasterConstitutiveTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterConstitutiveTypeMap.n_rows(), 2 ) );
//
//                // check mMasterConstitutiveTypeMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterConstitutiveTypeMap( 0, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterConstitutiveTypeMap( 1, 0 ) ),  0 ) );
//
//            //std::cout<<"Test create_constitutive_models"<<std::endl;
//            //------------------------------------------------------------------------------
//
//                Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 1 );
//                tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//                tConstitutiveUserDefinedInfo( 0 )( 0 ) = Constitutive_User_Defined_Info( { fem::Constitutive_Type::DIFF_LIN_ISO },
//                                                                                         {{ MSI::Dof_Type::TEMP }},
//                                                                                         { fem::Property_Type::CONDUCTIVITY } );
//
//                // create the properties for the set
//                tSet.create_constitutive_models( tConstitutiveUserDefinedInfo );
//
//                // check mMasterProperties size
//                CHECK( equal_to( tSet.mMasterCM.size(), 1 ) );
//
//                // check mMasterProperties content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterCM( 0 )->get_constitutive_type() ), 1 ) );
//
//                tSet.set_IWG_constitutive_models();
//
//            //std::cout<<"Test create_property_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create unique list of property type
//                tSet.create_property_type_list();
//
//                // check mPropertyTypeList size
//                CHECK( equal_to( tSet.mMasterPropTypes.size(), 2 ) );
//
//                // check mInterpDofAssemblyMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypes( 0 ) ), 1 ) ); // CONDUCTIVITY
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypes( 1 ) ), 3 ) ); // TEMP_DIRICHLET
//
//            //std::cout<<"Test create_property_type_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a property type map for the set
//                tSet.create_property_type_map();
//
//                // check mMasterPropTypeMap size
//                CHECK( equal_to( tSet.mMasterPropTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterPropTypeMap.n_rows(), 4 ) );
//
//                // check mMasterPropTypeMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 0, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 1, 0 ) ),  0 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 3, 0 ) ),  1 ) );
//
//            //std::cout<<"Test create_properties"<<std::endl;
//            //------------------------------------------------------------------------------
//
//                Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 1 );
//                tPropertyUserDefinedInfo( 0 ).resize( 2 );
//                tPropertyUserDefinedInfo( 0 )( 0 ) = Property_User_Defined_Info( { fem::Property_Type::CONDUCTIVITY },
//                                                                                 {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }},
//                                                                                 {{{ 1.0 }}},
//                                                                                 tConstValFunction_UTFEMSET,
//                                                                                 Cell< PropertyFunc >( 0 ) );
//                tPropertyUserDefinedInfo( 0 )( 1 ) = Property_User_Defined_Info( { fem::Property_Type::TEMP_DIRICHLET },
//                                                                                 Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                                 {{{ 5.0 }}},
//                                                                                 tConstValFunction_UTFEMSET,
//                                                                                 Cell< PropertyFunc >( 0 ) );
//
//                // create the properties for the set
//                tSet.create_properties( tPropertyUserDefinedInfo );
//
//                // check mMasterProperties size
//                CHECK( equal_to( tSet.mMasterProperties.size(), 2 ) );
//
//                // check mMasterProperties content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterProperties( 0 )->get_property_type() ), 1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterProperties( 1 )->get_property_type() ), 3 ) );
//
//            //std::cout<<"Test set_IWG_properties"<<std::endl;
//            //------------------------------------------------------------------------------
//                // set properties for each IWG
//                tSet.set_IWG_properties();
//
//                // check each IWG received right number of properties
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_properties().size(), 1 ) );
//                CHECK( equal_to( tSet.mIWGs( 1 )->get_properties().size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGs( 2 )->get_properties().size(), 0 ) );
//                CHECK( equal_to( tSet.mIWGs( 3 )->get_properties().size(), 0 ) );
//
//                // check each IWG received the right property type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_properties()( 0 )->get_property_type() ), 3 ) ); // CONDUCTIVITY
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 1 )->get_properties()( 0 )->get_property_type() ), 3 ) ); // CONDUCTIVITY
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 1 )->get_properties()( 1 )->get_property_type() ), 1 ) ); // TEMP_DIRICHLET
//
//            //std::cout<<"Test create_unique_dof_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // call create_uniaue_dof_type_lists
//                tSet.create_unique_dof_type_list();
//
//                // check the size of the list
//                CHECK( equal_to( tSet.mEqnObjDofTypeList.size(), 7 ) );
//
//                // check the content of the list
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 0 ) ),  0 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 1 ) ),  3 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 2 ) ),  6 ) ); //LS1
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 3 ) ),  8 ) ); //NLSX
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 4 ) ),  9 ) ); //NLSY
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 5 ) ), 10 ) ); //NLSZ
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 6 ) ), 11 ) ); //VX
//
//            //std::cout<<"Test create_dof_and_dv_type_lists"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create dof type list
//                tSet.create_dof_and_dv_type_lists();
//
//                // check the size of mInterpDofTypeList
//                CHECK( equal_to( tSet.mMasterDofTypes.size(), 5 ) );
//
//                // check the content of mInterpDofTypeList
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 0 )( 0 ) ),  3 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 1 )( 0 ) ),  0 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 2 )( 0 ) ), 11 ) ); //VX
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 3 )( 0 ) ),  8 ) ); //NLSX
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 3 )( 1 ) ),  9 ) ); //NLSY
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 3 )( 2 ) ), 10 ) ); //NLSZ
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 4 )( 0 ) ),  6 ) ); //LS1
//
//            //std::cout<<"Test create_dof_and_dv_type_maps"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a dof type map
//                tSet.create_dof_and_dv_type_maps();
//
//                // check mMasterDofTypeMap size
//                CHECK( equal_to( tSet.mMasterDofTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterDofTypeMap.n_rows(), 12 ) );
//
//                // check the content of mInterpDofTypeMap
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  0, 0 ) ),  1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  1, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  3, 0 ) ),  0 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  4, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  5, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  6, 0 ) ),  4 ) ); //LS1
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  7, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  8, 0 ) ),  3 ) ); //NLSX, NLSY, NLSZ
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  9, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap( 10, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap( 11, 0 ) ),  2 ) ); //VX
//
//            //std::cout<<"Mimic create_field_interpolators"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a cell of field interpolator pointers---------------------------------
//                moris::Cell< Field_Interpolator* > tFieldInterpolators( tSet.get_number_of_field_interpolators(), nullptr );
//
//                // set the number of coefficients for each field interpolator
//                tFieldInterpolators( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP } );
//                tFieldInterpolators( 1 ) = new Field_Interpolator( 1, { MSI::Dof_Type::UX } );
//                tFieldInterpolators( 2 ) = new Field_Interpolator( 1, { MSI::Dof_Type::VX } );
//                tFieldInterpolators( 3 ) = new Field_Interpolator( 3, { MSI::Dof_Type::NLSX,
//                                                                        MSI::Dof_Type::NLSY,
//                                                                        MSI::Dof_Type::NLSZ } );
//                tFieldInterpolators( 4 ) = new Field_Interpolator( 1, { MSI::Dof_Type::LS1 } );
//
//                // pass in the cell of FI pointers to the element block
//                tSet.mMasterFI = tFieldInterpolators;
//
//            //std::cout<<"Test create_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create dof assembly map
//                tSet.create_dof_assembly_map();
//
//                // check mInterpDofAssemblyMap size
//                CHECK( equal_to( tSet.mDofAssemblyMap.n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mDofAssemblyMap.n_rows(), 5 ) );
//
//                // check mInterpDofAssemblyMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 0, 0 ) ), 0 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 0, 1 ) ), 0 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 1, 0 ) ), 1 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 1, 1 ) ), 1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 2, 0 ) ), 2 ) ); //VX
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 2, 1 ) ), 2 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 3, 0 ) ), 3 ) ); //NLSX, NLSY, NLSZ
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 3, 1 ) ), 5 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 4, 0 ) ), 6 ) ); //LS1
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 4, 1 ) ), 6 ) );
//
//            //std::cout<<"Test create_IWG_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // set the IWG field inerpolators
//                tSet.set_IWG_field_interpolators();
//
//                // check each IWG received right number of field interpolators
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_dof_field_interpolators().size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGs( 1 )->get_dof_field_interpolators().size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGs( 2 )->get_dof_field_interpolators().size(), 1 ) );
//                CHECK( equal_to( tSet.mIWGs( 3 )->get_dof_field_interpolators().size(), 2 ) );
//
//                // check each IWG received the right dof type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ),  3 ) ); // TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ),  0 ) ); // UX
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 1 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ),  3 ) ); // TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 1 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ),  0 ) ); // UX
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 2 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ), 11 ) ); // VX
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 3 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ),  8 ) ); // NLSX, NLSY, NLSZ
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 3 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ),  6 ) ); // LS1
//
//            //std::cout<<"Test create_IWG_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create IWG dof assembly map
//                tSet.create_IWG_dof_assembly_map();
//
//                // check res dof assembly maps size
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap.size(), 4 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 ).n_rows(), 1 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 ).n_rows(), 1 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 2 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 2 ).n_rows(), 1 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 3 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 3 ).n_rows(), 1 ) );
//
//                // check resDofAssembly matrix per IWG
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 0, 1 ), 0 ) );
//
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 )( 0, 1 ), 0 ) );
//
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 2 )( 0, 0 ), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 2 )( 0, 1 ), 2 ) );
//
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 3 )( 0, 0 ), 3 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 3 )( 0, 1 ), 5 ) );
//
//                // check jac dof assembly maps size
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap.size(), 4 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 ).n_rows(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 ).n_rows(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 2 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 2 ).n_rows(), 1 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 ).n_rows(), 2 ) );
//
//                // check dofAssembly matrix per IWG
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 0, 1 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 1, 0 ), 1 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 1, 1 ), 1 ) );
//
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 0, 1 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 1, 0 ), 1 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 1, 1 ), 1 ) );
//
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 2 )( 0, 0 ), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 2 )( 0, 1 ), 2 ) );
//
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 )( 0, 0 ), 3 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 )( 0, 1 ), 5 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 )( 1, 0 ), 6 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 3 )( 1, 1 ), 6 ) );
//
//            //std::cout<<"Test create_IWG_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create IWG dof assembly map
//                tSet.set_properties_field_interpolators();
//
//                // check each IWG received right number of field interpolators
//                CHECK( equal_to( tSet.mMasterProperties( 0 )->get_dof_field_interpolators().size(), 0 ) );
//                CHECK( equal_to( tSet.mMasterProperties( 1 )->get_dof_field_interpolators().size(), 2 ) );

        }/* TEST_CASE */

        // This test check all the member functions of the FEM Set in the case where
        // there is a master and a slave.
        TEST_CASE( "FEMSet_Slave", "[moris],[fem],[FEMSet_Slave]" )
        {
//            //create a set
//            Set tSet;
//
//            // list of IWG types
//            moris::Cell< fem::IWG_Type > tIWGTypeList= { fem::IWG_Type::SPATIALDIFF_GHOST,
//                                                         fem::IWG_Type::HELMHOLTZ };
//
//            // number of IWGs to be created
//            uint tNumOfIWGs = tIWGTypeList.size();
//
//            // list of residual dof type
//            moris::Cell< moris::Cell< MSI::Dof_Type > > aResidualDofType( tNumOfIWGs );
//            aResidualDofType( 0 ) = { MSI::Dof_Type::TEMP };
//            aResidualDofType( 1 ) = { MSI::Dof_Type::VX };
//
//            // list of active dof type
//            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > aMasterDofTypes( tNumOfIWGs );
//            aMasterDofTypes( 0 ) = {{ MSI::Dof_Type::TEMP }};
//            aMasterDofTypes( 1 ) = {{ MSI::Dof_Type::VX }};
//            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > aSlaveDofTypes( tNumOfIWGs );
//            aSlaveDofTypes( 0 ) = {{ MSI::Dof_Type::TEMP }};
//
//            // list of active property type
//            moris::Cell< moris::Cell< fem::Property_Type > > aMasterPropTypes( tNumOfIWGs );
//            aMasterPropTypes( 0 ) = { fem::Property_Type::CONDUCTIVITY };
//            moris::Cell< moris::Cell< fem::Property_Type > > aSlavePropTypes( tNumOfIWGs );
//            aSlavePropTypes( 0 ) = { fem::Property_Type::CONDUCTIVITY };
//
//            moris::Cell< moris::Cell< fem::Constitutive_Type > > aMasterConstitutiveTypes( tNumOfIWGs );
//
//            // create a cell of IWGs for the problem considered
//            moris::Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );
//
//            // a factory to create the IWGs
//            fem::IWG_Factory tIWGFactory;
//
//            // loop over the IWG types
//            for( uint i = 0; i < tNumOfIWGs; i++)
//            {
//                // create an IWG with the factory for the ith IWG type
//                tIWGs( i ) = tIWGFactory.create_IWG( tIWGTypeList( i ) );
//
//                // set residual dof type
//                tIWGs( i )->set_residual_dof_type( aResidualDofType( i ) );
//
//                // set active master and slave dof type
//                tIWGs( i )->set_dof_type_list( aMasterDofTypes( i ) );
//                tIWGs( i )->set_dof_type_list( aSlaveDofTypes( i ), mtk::Master_Slave::SLAVE );
//
//                // set active master and slave property type
//                tIWGs( i )->set_property_type_list( aMasterPropTypes( i ) );
//                tIWGs( i )->set_property_type_list( aSlavePropTypes( i ), mtk::Master_Slave::SLAVE );
//
//                // set active constitutive type
//                tIWGs( i )->set_constitutive_type_list( aMasterConstitutiveTypes( i ) );
//                tIWGs( i )->set_constitutive_type_list( aMasterConstitutiveTypes( i ), mtk::Master_Slave::SLAVE );
//            }
//
//            // pass in the cell of IWG pointers to the element block
//            tSet.mIWGs = tIWGs;
//
//            //std::cout<<"Test create_constitutive_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create unique list of constitutive type
//                tSet.create_constitutive_type_list();
//
//                // check mMasterConstitutiveTypes size
//                CHECK( equal_to( tSet.mMasterConstitutiveTypes.size(), 0 ) );
//                CHECK( equal_to( tSet.mSlaveConstitutiveTypes.size(), 0 ) );
//
//            //std::cout<<"Test create_constitutive_type_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a constitutive type map for the set
//                tSet.create_constitutive_type_map();
//
//                // check mMasterConstitutiveTypeMap size
//                CHECK( equal_to( tSet.mMasterConstitutiveTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterConstitutiveTypeMap.n_rows(), 0 ) );
//
//            //std::cout<<"Test create_constitutive_models"<<std::endl;
//            //------------------------------------------------------------------------------
//
//                Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 2 );
//
//                // create the properties for the set
//                tSet.create_constitutive_models( tConstitutiveUserDefinedInfo );
//
//                // check mMasterProperties size
//                CHECK( equal_to( tSet.mMasterCM.size(), 0 ) );
//                CHECK( equal_to( tSet.mSlaveCM.size(), 0 ) );
//
//                tSet.set_IWG_constitutive_models();
//
//            //std::cout<<"Test create_property_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create unique list of property type
//                tSet.create_property_type_list();
//
//                // check mMasterPropTypes size
//                CHECK( equal_to( tSet.mMasterPropTypes.size(), 1 ) );
//
//                // check mMasterPropTypes content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypes( 0 ) ), 3 ) );
//
//                // check mSlavePropTypes size
//                CHECK( equal_to( tSet.mSlavePropTypes.size(), 1 ) );
//
//                // check mMasterPropTypes content
//                CHECK( equal_to( static_cast< uint >( tSet.mSlavePropTypes( 0 ) ), 3 ) );
//
//            //std::cout<<"Test create_property_type_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a property type map for the set
//                tSet.create_property_type_map();
//
//                // check mMasterPropTypeMap size
//                CHECK( equal_to( tSet.mMasterPropTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterPropTypeMap.n_rows(), 4 ) );
//
//                // check mMasterPropTypeMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 0, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 1, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterPropTypeMap( 3, 0 ) ),  0 ) );
//
//                // check mSlavePropTypeMap size
//                CHECK( equal_to( tSet.mSlavePropTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mSlavePropTypeMap.n_rows(), 4 ) );
//
//                // check mSlavePropTypeMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mSlavePropTypeMap( 0, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlavePropTypeMap( 1, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlavePropTypeMap( 2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlavePropTypeMap( 3, 0 ) ),  0 ) );
//
//            //std::cout<<"Test create_properties"<<std::endl;
//            //------------------------------------------------------------------------------
//
//                // create a property user defined info container
//                Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 2 );
//
//                // fill the master property user defined info container
//                tPropertyUserDefinedInfo( 0 ).resize( 1 );
//                tPropertyUserDefinedInfo( 0 )( 0 ) = Property_User_Defined_Info( { fem::Property_Type::CONDUCTIVITY },
//                                                                                 {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }},
//                                                                                 {{{ 1.0 }}},
//                                                                                 tConstValFunction_UTFEMSET,
//                                                                                 { tConstValFunction_UTFEMSET, tConstValFunction_UTFEMSET } );
//                // fill the slave property user defined info container
//                tPropertyUserDefinedInfo( 1 ).resize( 1 );
//                tPropertyUserDefinedInfo( 1 )( 0 ) = Property_User_Defined_Info( { fem::Property_Type::CONDUCTIVITY },
//                                                                                 {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }},
//                                                                                 {{{ 1.0 }}},
//                                                                                 tConstValFunction_UTFEMSET,
//                                                                                 { tConstValFunction_UTFEMSET, tConstValFunction_UTFEMSET } );
//
//                // create the properties for the set
//                tSet.create_properties( tPropertyUserDefinedInfo );
//
//                // check mMasterProperties size
//                CHECK( equal_to( tSet.mMasterProperties.size(), 1 ) );
//
//                // check mMasterProperties content
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterProperties( 0 )->get_property_type() ), 3 ) );
//
//                // check mSlaveProperties size
//                CHECK( equal_to( tSet.mSlaveProperties.size(), 1 ) );
//
//                // check mMasterProperties content
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveProperties( 0 )->get_property_type() ), 3 ) );
//
//            //std::cout<<"Test set_IWG_properties"<<std::endl;
//            //------------------------------------------------------------------------------
//                // set properties for each IWG
//                tSet.set_IWG_properties();
//
//                // check each IWG received right number of master properties
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_properties().size(), 1 ) );
//                CHECK( equal_to( tSet.mIWGs( 1 )->get_properties().size(), 0 ) );
//
//                // check each IWG received the right master property type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_properties()( 0 )->get_property_type() ), 3 ) ); // CONDUCTIVITY
//
//                // check each IWG received right number of slave properties
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_properties( mtk::Master_Slave::SLAVE ).size(), 1 ) );
//
//                // check each IWG received the right slave property type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_properties( mtk::Master_Slave::SLAVE )( 0 )->get_property_type() ), 3 ) ); // CONDUCTIVITY
//
//            //std::cout<<"Test create_unique_dof_type_list"<<std::endl;
//            //------------------------------------------------------------------------------
//                // call create_uniaue_dof_type_lists
//                tSet.create_unique_dof_type_list();
//
//                // check the size of the list
//                CHECK( equal_to( tSet.mEqnObjDofTypeList.size(), 3 ) );
//
//                // check the content of the list
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 0 ) ),  0 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 1 ) ),  3 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mEqnObjDofTypeList( 2 ) ), 11 ) ); //VX
//
//            //std::cout<<"Test create_dof_and_dv_type_lists"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create dof type list
//                tSet.create_dof_and_dv_type_lists();
//
//                // check the size of mMasterDofTypes
//                CHECK( equal_to( tSet.mMasterDofTypes.size(), 3 ) );
//
//                // check the content of mMasterDofTypes
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 0 )( 0 ) ),  3 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 1 )( 0 ) ),  0 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypes( 2 )( 0 ) ), 11 ) ); //VX
//
//                // check the size of mSlaveDofTypes
//                CHECK( equal_to( tSet.mSlaveDofTypes.size(), 2 ) );
//
//                // check the content of mSlaveDofTypes
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypes( 0 )( 0 ) ),  3 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypes( 1 )( 0 ) ),  0 ) ); //UX
//
//            //std::cout<<"Test create_dof_and_dv_type_maps"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a dof type map
//                tSet.create_dof_and_dv_type_maps();
//
//                // check mMasterDofTypeMap size
//                CHECK( equal_to( tSet.mMasterDofTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mMasterDofTypeMap.n_rows(), 12 ) );
//
//                // check the content of mInterpDofTypeMap
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  0, 0 ) ),  1 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  1, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  3, 0 ) ),  0 ) ); //TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  4, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  5, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  6, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  7, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  8, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap(  9, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap( 10, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterDofTypeMap( 11, 0 ) ),  2 ) ); //VX
//
//                // check mSlaveDofTypeMap size
//                CHECK( equal_to( tSet.mSlaveDofTypeMap.n_cols(), 1 ) );
//                CHECK( equal_to( tSet.mSlaveDofTypeMap.n_rows(), 4 ) );
//
//                // check the content of mSlaveDofTypeMap
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypeMap(  0, 0 ) ),  1 ) ); //UX
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypeMap(  1, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypeMap(  2, 0 ) ), -1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveDofTypeMap(  3, 0 ) ),  0 ) ); //TEMP
//
//
//            //std::cout<<"Mimic create_field_interpolators"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create a cell of master field interpolator pointers
//                moris::Cell< Field_Interpolator* > tMasterFI( tSet.get_number_of_field_interpolators(), nullptr );
//
//                // set the number of coefficients for each field interpolator
//                tMasterFI( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP } );
//                tMasterFI( 1 ) = new Field_Interpolator( 1, { MSI::Dof_Type::UX } );
//                tMasterFI( 2 ) = new Field_Interpolator( 1, { MSI::Dof_Type::VX } );
//
//                // pass in the cell of FI pointers to the element block                print(tSet.mIWGDofAssemblyMap( 1 ),"mIWGDofAssemblyMap( 1 )");
//
//                tSet.mMasterFI = tMasterFI;
//
//                // create a cell of master field interpolator pointers
//                moris::Cell< Field_Interpolator* > tSlaveFI( tSet.get_number_of_field_interpolators( mtk::Master_Slave::SLAVE ), nullptr );
//
//                // set the number of coefficients for each field interpolator
//                tSlaveFI( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP } );
//                tSlaveFI( 1 ) = new Field_Interpolator( 1, { MSI::Dof_Type::UX } );
//
//                // pass in the cell of FI pointers to the element block
//                tSet.mSlaveFI = tSlaveFI;
//
//            //std::cout<<"Test create_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create dof assembly map
//                tSet.create_dof_assembly_map();
//
//                // check mInterpDofAssemblyMap size
//                CHECK( equal_to( tSet.mDofAssemblyMap.n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mDofAssemblyMap.n_rows(), 5 ) );
//
//                // check mInterpDofAssemblyMap content
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 0, 0 ) ), 0 ) ); //TEMP MASTER
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 0, 1 ) ), 0 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 1, 0 ) ), 1 ) ); //UX MASTER
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 1, 1 ) ), 1 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 2, 0 ) ), 2 ) ); //VX MASTER
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 2, 1 ) ), 2 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 3, 0 ) ), 3 ) ); //TEMP SLAVE
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 3, 1 ) ), 3 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 4, 0 ) ), 4 ) ); //UX SLAVE
//                CHECK( equal_to( static_cast< uint >( tSet.mDofAssemblyMap( 4, 1 ) ), 4 ) );
//
//            //std::cout<<"Test set_IWG_field_interpolators"<<std::endl;
//            //------------------------------------------------------------------------------
//                // set the IWG field inerpolators
//                tSet.set_IWG_field_interpolators();
//
//                // check each IWG received right number of master field interpolators
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_dof_field_interpolators().size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGs( 1 )->get_dof_field_interpolators().size(), 1 ) );
//
//                // check each IWG received the right master dof type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ),  3 ) ); // TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ),  0 ) ); // UX
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 1 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ), 11 ) ); // VX
//
//                // check each IWG received right number of slave field interpolators
//                CHECK( equal_to( tSet.mIWGs( 0 )->get_dof_field_interpolators( mtk::Master_Slave::SLAVE ).size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGs( 1 )->get_dof_field_interpolators( mtk::Master_Slave::SLAVE ).size(), 0 ) );
//
//                // check each IWG received the right slave dof type
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators( mtk::Master_Slave::SLAVE )( 0 )->get_dof_type()( 0 ) ),  3 ) ); // TEMP
//                CHECK( equal_to( static_cast< uint >( tSet.mIWGs( 0 )->get_dof_field_interpolators( mtk::Master_Slave::SLAVE )( 1 )->get_dof_type()( 0 ) ),  0 ) ); // UX
//
//            //std::cout<<"Test create_IWG_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create IWG dof assembly map
//                tSet.create_IWG_dof_assembly_map();
//
//                // check mIWGResDofAssemblyMap size
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap.size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 ).n_rows(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 ).n_rows(), 1 ) );
//
//                // check dofAssembly matrix per IWG
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 0, 1 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 1, 0 ), 3 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 0 )( 1, 1 ), 3 ) );
//
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 )( 0, 0 ), 2 ) );
//                CHECK( equal_to( tSet.mIWGResDofAssemblyMap( 1 )( 0, 1 ), 2 ) );
//
//                // check mIWGResDofAssemblyMap size
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap.size(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 ).n_rows(), 4 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 ).n_cols(), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 ).n_rows(), 1 ) );
//
//                // check dofAssembly matrix per IWG
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 0, 0 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 0, 1 ), 0 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 1, 0 ), 1 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 1, 1 ), 1 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 2, 0 ), 3 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 2, 1 ), 3 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 3, 0 ), 4 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 0 )( 3, 1 ), 4 ) );
//
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 0, 0 ), 2 ) );
//                CHECK( equal_to( tSet.mIWGJacDofAssemblyMap( 1 )( 0, 1 ), 2 ) );
//
//            //std::cout<<"Test create_IWG_dof_assembly_map"<<std::endl;
//            //------------------------------------------------------------------------------
//                // create IWG dof assembly map
//                tSet.set_properties_field_interpolators();
//
//                // check each IWG received right number of master field interpolators
//                CHECK( equal_to( tSet.mMasterProperties( 0 )->get_dof_field_interpolators().size(), 2 ) );
//
//                // check each IWG received right master dof types
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterProperties( 0 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ), 3 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mMasterProperties( 0 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ), 0 ) );
//
//                // check each IWG received right number of slave field interpolators
//                CHECK( equal_to( tSet.mSlaveProperties( 0 )->get_dof_field_interpolators().size(), 2 ) );
//
//                // check each IWG received right slave dof types
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveProperties( 0 )->get_dof_field_interpolators()( 0 )->get_dof_type()( 0 ) ), 3 ) );
//                CHECK( equal_to( static_cast< uint >( tSet.mSlaveProperties( 0 )->get_dof_field_interpolators()( 1 )->get_dof_type()( 0 ) ), 0 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */

