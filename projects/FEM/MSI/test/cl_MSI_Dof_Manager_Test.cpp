/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Dof_Manager_Test.cpp
 *
 */

#include "catch.hpp"
#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#define protected public
#define private   public
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#undef protected
#undef private

#include "fn_PRM_MSI_Parameters.hpp"
namespace moris
{
    namespace MSI
    {
    TEST_CASE("Dof_Manager_Max_Pdof_Host","[MSI],[Dof_max_pdof_hosts]")
    {
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 1;

        fem::Node_Base * Node1;
        fem::Node_Base * Node2;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsId1( 2, 1 );
        Matrix< IdMat > tAdofsId2( 2, 1 );

        tAdofsId1( 0, 0 ) = 0;
        tAdofsId1( 1, 0 ) = 1;
        tAdofsId2( 0, 0 ) = 0;
        tAdofsId2( 1, 0 ) = 1;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsInd1( 2, 1 );
        Matrix< IdMat > tAdofsInd2( 2, 1 );

        tAdofsInd1( 0, 0 ) = 0;
        tAdofsInd1( 1, 0 ) = 1;
        tAdofsInd2( 0, 0 ) = 0;
        tAdofsInd2( 1, 0 ) = 1;

        // Create generic T-matrices
        Matrix< DDRMat > tMatrix1( 2, 1 );
        Matrix< DDRMat > tMatrix2( 2, 1 );

        // Create generic T-matrices
        tMatrix1( 0, 0 ) = 1.0;
        tMatrix1( 1, 0 ) = 1.0;
        tMatrix2( 0, 0 ) = 1.0;
        tMatrix2( 1, 0 ) = -2.0;

        // Create generic adof owning processor
        Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
        Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

        tAdofOwningProcessor1( 0, 0 ) = 0;
        tAdofOwningProcessor1( 1, 0 ) = 0;
        tAdofOwningProcessor2( 0, 0 ) = 0;
        tAdofOwningProcessor2( 1, 0 ) = 0;

        // Create generic Node Object
        Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
        Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );

        moris::uint tNumEquationObjects = 2;

        moris::uint tNumNodes = 2;

        moris::Cell < Equation_Object* >tListEqnObj( tNumEquationObjects, nullptr );

        // Create List with node pointern correponding to generic equation object
        moris::Cell< moris::Cell< fem::Node_Base * > > tNodeIds_1( 1 );
        tNodeIds_1( 0 ).resize( tNumNodes );
        tNodeIds_1( 0 )( 0 ) = Node1;
        tNodeIds_1( 0 )( 1 ) = Node2;

        moris::Cell< moris::Cell< fem::Node_Base * > > tNodeIds_2( 1 );
        tNodeIds_2( 0 ).resize( tNumNodes );
        tNodeIds_2( 0 )( 0 ) = Node1;
        tNodeIds_2( 0 )( 1 ) = Node2;

        // Create generic equation objects
        Equation_Object EquObj_1( tNodeIds_1 );
        Equation_Object EquObj_2( tNodeIds_2 );

        // Create List with equation objects
        tListEqnObj( 0 ) = & EquObj_1;
        tListEqnObj( 1 ) = & EquObj_2;

        Matrix< IdMat > tCommTable( 1, 1, 0 );

        Dof_Manager tDofMgn;

        CHECK( equal_to( tDofMgn.initialize_max_number_of_possible_pdof_hosts( tListEqnObj ), 2 ) );

        delete tNodeIds_1( 0 )( 0 );
        delete tNodeIds_1( 0 )( 1 );
    }

    TEST_CASE("Dof_Manager_Pdof_Host_Time_Level","[MSI],[Dof_time_level]")
    {
        // Create dof maager
        Dof_Manager tDofMgn;

        // Set pdof type list
        tDofMgn.mPdofTypeList.resize( 2 );
        tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;
        tDofMgn.mPdofTypeList( 1 ) = Dof_Type::UX;

        // create pdof host list and fill it with needed information
        tDofMgn.mPdofHostList.resize( 2 );
        tDofMgn.mPdofHostList( 0 ) = new Pdof_Host();
        tDofMgn.mPdofHostList( 1 ) = new Pdof_Host();

        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 2 );
        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 3 );
        tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType.resize( 2 );
        tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 ).resize( 2 );

        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 1 )) = new Pdof;
        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;
        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 1 )) = new Pdof;
        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 2 )) = new Pdof;
        (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
        (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 1 )) = new Pdof;

        // Call the function which shall be tested
        tDofMgn.initialize_pdof_host_time_level_list();

        // Check length and entries of the resulting vector
        CHECK( equal_to( tDofMgn.mPdofHostTimeLevelList.length(), 2 ) );

        CHECK( equal_to( tDofMgn.mPdofHostTimeLevelList( 0, 0 ), 2 ) );
        CHECK( equal_to( tDofMgn.mPdofHostTimeLevelList( 1, 0 ), 3 ) );
    }

    TEST_CASE("Dof_Mgn_ini_pdof_host_list","[MSI],[Dof_ini_pdof_host_list]")
    {
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 1;

        fem::Node_Base * Node1;
        fem::Node_Base * Node2;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsId1( 2, 1 );
        Matrix< IdMat > tAdofsId2( 2, 1 );

        tAdofsId1( 0, 0 ) = 0;
        tAdofsId1( 1, 0 ) = 1;
        tAdofsId2( 0, 0 ) = 0;
        tAdofsId2( 1, 0 ) = 1;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsInd1( 2, 1 );
        Matrix< IdMat > tAdofsInd2( 2, 1 );

        tAdofsInd1( 0, 0 ) = 0;
        tAdofsInd1( 1, 0 ) = 1;
        tAdofsInd2( 0, 0 ) = 0;
        tAdofsInd2( 1, 0 ) = 1;

        // Create generic T-matrices
        Matrix< DDRMat > tMatrix1( 2, 1 );
        Matrix< DDRMat > tMatrix2( 2, 1 );

        // Create generic T-matrices
        tMatrix1( 0, 0 ) = 1.0;
        tMatrix1( 1, 0 ) = 1.0;
        tMatrix2( 0, 0 ) = 1.0;
        tMatrix2( 1, 0 ) = -2.0;

        // Create generic adof owning processor
        Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
        Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

        tAdofOwningProcessor1( 0, 0 ) = 0;
        tAdofOwningProcessor1( 1, 0 ) = 0;
        tAdofOwningProcessor2( 0, 0 ) = 0;
        tAdofOwningProcessor2( 1, 0 ) = 0;

        // Create generic Node Object
        Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
        Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );

        moris::uint tNumEquationObjects = 2;

        moris::uint tNumNodes = 2;

        moris::Cell < Equation_Object* >tListEqnObj;
        tListEqnObj.resize( tNumEquationObjects, nullptr );

        // Create List with node pointern correponding to generic equation object
        moris::Cell< moris::Cell< fem::Node_Base * > > tNodeIds_1( 1 );
        tNodeIds_1( 0 ).resize( tNumNodes );
        tNodeIds_1( 0 )( 0 ) = Node1;
        tNodeIds_1( 0 )( 1 ) = Node2;

        moris::Cell< moris::Cell< fem::Node_Base * > > tNodeIds_2( 1 );
        tNodeIds_2( 0 ).resize( tNumNodes );
        tNodeIds_2( 0 )( 0 ) = Node1;
        tNodeIds_2( 0 )( 1 ) = Node2;

        // Create generic equation objects
        Equation_Object EquObj_1( tNodeIds_1 );
        Equation_Object EquObj_2( tNodeIds_2 );

        Equation_Set tEqnBlock;
        tEqnBlock.mUniqueDofTypeList.resize( 1, MSI::Dof_Type::TEMP );
        tEqnBlock.mUniqueDofTypeListMasterSlave.resize( 1 );
        tEqnBlock.mUniqueDofTypeListMasterSlave(0).resize( 1, MSI::Dof_Type::TEMP );
        EquObj_1.mEquationSet = &tEqnBlock;
        EquObj_2.mEquationSet = &tEqnBlock;

        // Create List with equation objects
        tListEqnObj( 0 ) = & EquObj_1;
        tListEqnObj( 1 ) = & EquObj_2;

        Dof_Manager tDofMgn;

        tDofMgn.mPdofTypeList.resize(1, Dof_Type::TEMP);
        tDofMgn.mPdofTypeMap.set_size( 4, 1 );
        tDofMgn.mPdofTypeMap( 3, 0 ) = 0;

        tDofMgn.mTimePerDofType.set_size( 1, 1, 1 );

        tDofMgn.initialize_pdof_host_list( tListEqnObj );

        // Check size of pdof host list
        CHECK( equal_to( tDofMgn.mPdofHostList.size(), 2 ) );

        CHECK( equal_to( ((tDofMgn.mPdofHostList(0))->mNodeObj)->get_id(), 0 ) );
        CHECK( equal_to( ((tDofMgn.mPdofHostList(1))->mNodeObj)->get_id(), 1 ) );

        delete Node1;
        delete Node2;
    }

    TEST_CASE("Dof_Mgn_create_adofs","[MSI],[Dof_create_adofs]")
    {
        if( par_size() == 1 )
        {
            // Create node obj
            moris::uint tNodeId1 = 0;
            moris::uint tNodeId2 = 1;

            fem::Node_Base * Node1;
            fem::Node_Base * Node2;

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsId1( 2, 1 );
            Matrix< IdMat > tAdofsId2( 2, 1 );

            tAdofsId1( 0, 0 ) = 0;
            tAdofsId1( 1, 0 ) = 5;
            tAdofsId2( 0, 0 ) = 3;
            tAdofsId2( 1, 0 ) = 0;

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsInd1( 2, 1 );
            Matrix< IdMat > tAdofsInd2( 2, 1 );

            tAdofsInd1( 0, 0 ) = 0;
            tAdofsInd1( 1, 0 ) = 5;
            tAdofsInd2( 0, 0 ) = 3;
            tAdofsInd2( 1, 0 ) = 0;

            // Create generic T-matrices
            Matrix< DDRMat > tMatrix1( 2, 1 );
            Matrix< DDRMat > tMatrix2( 2, 1 );

            // Create generic T-matrices
            tMatrix1( 0, 0 ) = 1.0;
            tMatrix1( 1, 0 ) = -4.0;
            tMatrix2( 0, 0 ) = 2.0;
            tMatrix2( 1, 0 ) = -2.0;

            // Create generic adof owning processor
            Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
            Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

            tAdofOwningProcessor1( 0, 0 ) = 0;
            tAdofOwningProcessor1( 1, 0 ) = 0;
            tAdofOwningProcessor2( 0, 0 ) = 0;
            tAdofOwningProcessor2( 1, 0 ) = 0;

            // Create generic Node Object
            Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
            Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );

            // Create dof manager and hardcode initial values
            Dof_Manager tDofMgn;

            tDofMgn.mNumMaxAdofs = -1;

            tDofMgn.mUseHMR = true;

            tDofMgn.mPdofTypeList.resize( 2 );
            tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;
            tDofMgn.mPdofTypeList( 1 ) = Dof_Type::UX;

            tDofMgn.mPdofHostList.resize( 2 );
            tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );
            tDofMgn.mPdofHostList( 1 ) = new Pdof_Host( 2, Node2 );

            tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
            tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
            tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );
            tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType.resize( 2 );
            tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 ).resize( 1 );

            (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
            (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;
            (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;

            (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;
            (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 ))->mDofTypeIndex = 1;
            (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;

            tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId1;
            tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds = tAdofsId1;
            tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId2;

            moris::Cell < Equation_Object* >tListEqnObj;
            moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
            tDofMgn.mModelSolverInterface = new Model_Solver_Interface( tMSIParameters, tListEqnObj );
            tDofMgn.mModelSolverInterface->mDofMgn = tDofMgn;
            // end hardcoding stuff

            // Create adofs and build adof lists
            tDofMgn.create_adofs();

            CHECK( equal_to( tDofMgn.mAdofList.size(), 5 ) );
            CHECK( equal_to( tDofMgn.mAdofList( 0 )->mAdofId, 0 ) );
            CHECK( equal_to( tDofMgn.mAdofList( 4 )->mAdofId, 4 ) );

            CHECK( equal_to( tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds( 0 ), 0 ) );
            CHECK( equal_to( tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds( 1 ), 2 ) );
            CHECK( equal_to( tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds( 0 ), 3 ) );
            CHECK( equal_to( tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds( 1 ), 4 ) );
            CHECK( equal_to( tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds( 0 ), 1 ) );
            CHECK( equal_to( tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds( 1 ), 0 ) );

            delete Node1;
            delete Node2;
        }
    }

    TEST_CASE("Dof_Mgn_Set_T_Matrix","[MSI],[Dof_set_t_matrix]")
    {
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 1;

        fem::Node_Base * Node1;
        fem::Node_Base * Node2;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsId1( 2, 1 );
        Matrix< IdMat > tAdofsId2( 2, 1 );

        tAdofsId1( 0, 0 ) = 0;
        tAdofsId1( 1, 0 ) = 5;
        tAdofsId2( 0, 0 ) = 3;
        tAdofsId2( 1, 0 ) = 0;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsInd1( 2, 1 );
        Matrix< IdMat > tAdofsInd2( 2, 1 );

        tAdofsInd1( 0, 0 ) = 0;
        tAdofsInd1( 1, 0 ) = 5;
        tAdofsInd2( 0, 0 ) = 3;
        tAdofsInd2( 1, 0 ) = 0;

        // Create generic T-matrices
        Matrix< DDRMat > tMatrix1( 2, 1 );
        Matrix< DDRMat > tMatrix2( 2, 1 );

        // Create generic T-matrices
        tMatrix1( 0, 0 ) = 1.0;
        tMatrix1( 1, 0 ) = -4.0;
        tMatrix2( 0, 0 ) = 2.0;
        tMatrix2( 1, 0 ) = -2.0;

        // Create generic adof owning processor
        Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
        Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

        tAdofOwningProcessor1( 0, 0 ) = 0;
        tAdofOwningProcessor1( 1, 0 ) = 0;
        tAdofOwningProcessor2( 0, 0 ) = 0;
        tAdofOwningProcessor2( 1, 0 ) = 0;

        // Create generic Node Object
        Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1,  tMatrix1, tAdofOwningProcessor1 );
        Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );

        // Create dof manager and hardcode initial values
        Dof_Manager tDofMgn;

        tDofMgn.mPdofTypeList.resize( 2 );
        tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;
        tDofMgn.mPdofTypeList( 1 ) = Dof_Type::UX;

        tDofMgn.mNumMaxAdofs = -1;

        tDofMgn.mUseHMR = true;

        tDofMgn.mPdofHostList.resize( 2 );
        tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );
        tDofMgn.mPdofHostList( 1 ) = new Pdof_Host( 2, Node2 );

        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
        tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );
        tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType.resize( 1 );
        tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 ).resize( 1 );

        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
        (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;
        (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;

        moris::Cell < Equation_Object* >tListEqnObj;
        moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
        tDofMgn.mModelSolverInterface = new Model_Solver_Interface( tMSIParameters, tListEqnObj );
        tDofMgn.mModelSolverInterface->mDofMgn = tDofMgn;
        // end hardcoding stuff

        // Create adofs and build adof lists
        tDofMgn.set_pdof_t_matrix();

        CHECK( equal_to( (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mTmatrix)( 0, 0 ),  1 ) );
        CHECK( equal_to( (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mTmatrix)( 1, 0 ), -4 ) );
        CHECK( equal_to( (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mTmatrix)( 0, 0 ),  1 ) );
        CHECK( equal_to( (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mTmatrix)( 1, 0 ), -4 ) );
        CHECK( equal_to( (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mTmatrix)( 0, 0 ),  2 ) );
        CHECK( equal_to( (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mTmatrix)( 1, 0 ), -2 ) );

        delete Node1;
        delete Node2;
    }

    TEST_CASE("Dof_Mgn_create_unique_dof_type_list","[MSI],[Dof_create_dof_and_dv_type_lists][MSI_parallel]")
     {
         // Create generic equation objects
         Equation_Set EquObj_1;
         Equation_Set EquObj_2;

         moris::Cell < Equation_Set* >tListEqnObj;

         // Determine process rank
         size_t tRank = par_rank();

         // Hardcode input test values
         switch( tRank )
             {
             case 0:
                 EquObj_1.mUniqueDofTypeList.resize( 1 );
                 EquObj_2.mUniqueDofTypeList.resize( 2 );
                 EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                 EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                 EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                 tListEqnObj.resize( 2, nullptr );
                 tListEqnObj( 0 ) = & EquObj_1;
                 tListEqnObj( 1 ) = & EquObj_2;
               break;
             case 1:
                 EquObj_1.mUniqueDofTypeList.resize( 2 );
                 EquObj_2.mUniqueDofTypeList.resize( 2 );
                 EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                 EquObj_1.mUniqueDofTypeList( 1 ) = Dof_Type::UX;
                 EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                 EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                 tListEqnObj.resize( 2, nullptr );
                 tListEqnObj( 0 ) = & EquObj_1;
                 tListEqnObj( 1 ) = & EquObj_2;
               break;
             case 2:
                 EquObj_1.mUniqueDofTypeList.resize( 3 );
                 EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                 EquObj_1.mUniqueDofTypeList( 1 ) = Dof_Type::TEMP;
                 EquObj_1.mUniqueDofTypeList( 2 ) = Dof_Type::UZ;
                 tListEqnObj.resize( 1, nullptr );
                 tListEqnObj( 0 ) = & EquObj_1;
               break;
             case 3:
                 EquObj_1.mUniqueDofTypeList.resize( 1 );
                 EquObj_2.mUniqueDofTypeList.resize( 2 );
                 EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                 EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                 EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                 tListEqnObj.resize( 2, nullptr );
                 tListEqnObj( 0 ) = & EquObj_1;
                 tListEqnObj( 1 ) = & EquObj_2;
               break;
              }

         // Create dof manager
         Dof_Manager tDofMgn;

         // Call initialize pdof type list function
         tDofMgn.initialize_pdof_type_list( tListEqnObj );

         // Check pdof type list
         CHECK( equal_to( static_cast<int>( tDofMgn.mPdofTypeList( 0 ) ), 0 ) );
         CHECK( equal_to( static_cast<int>( tDofMgn.mPdofTypeList( 1 ) ), 2 ) );
         CHECK( equal_to( static_cast<int>( tDofMgn.mPdofTypeList( 2 ) ), 3 ) );

     }

    TEST_CASE("Dof_Mgn_create_unique_dof_type_map_matrix","[MSI],[Dof_create_dof_type_map][MSI_parallel]")
    {
        // Create generic equation objects
        Equation_Set EquObj_1;
        Equation_Set EquObj_2;

        moris::Cell < Equation_Set* >tListEqnObj;

        // Determine process rank
        size_t tRank = par_rank();

        // Hardcode input test values
        switch( tRank )
            {
            case 0:
                EquObj_1.mUniqueDofTypeList.resize( 1 );
                EquObj_2.mUniqueDofTypeList.resize( 2 );
                EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                tListEqnObj.resize( 2, nullptr );
                tListEqnObj( 0 ) = & EquObj_1;
                tListEqnObj( 1 ) = & EquObj_2;
              break;
            case 1:
                EquObj_1.mUniqueDofTypeList.resize( 2 );
                EquObj_2.mUniqueDofTypeList.resize( 2 );
                EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                EquObj_1.mUniqueDofTypeList( 1 ) = Dof_Type::UX;
                EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                tListEqnObj.resize( 2, nullptr );
                tListEqnObj( 0 ) = & EquObj_1;
                tListEqnObj( 1 ) = & EquObj_2;
              break;
            case 2:
                EquObj_1.mUniqueDofTypeList.resize( 3 );
                EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                EquObj_1.mUniqueDofTypeList( 1 ) = Dof_Type::TEMP;
                EquObj_1.mUniqueDofTypeList( 2 ) = Dof_Type::UZ;
                tListEqnObj.resize( 1, nullptr );
                tListEqnObj( 0 ) = & EquObj_1;
              break;
            case 3:
                EquObj_1.mUniqueDofTypeList.resize( 1 );
                EquObj_2.mUniqueDofTypeList.resize( 2 );
                EquObj_1.mUniqueDofTypeList( 0 ) = Dof_Type::TEMP;
                EquObj_2.mUniqueDofTypeList( 0 ) = Dof_Type::UX;
                EquObj_2.mUniqueDofTypeList( 1 ) = Dof_Type::UZ;
                tListEqnObj.resize( 2, nullptr );
                tListEqnObj( 0 ) = & EquObj_1;
                tListEqnObj( 1 ) = & EquObj_2;
              break;
             }

        // Create dof manager
        Dof_Manager tDofMgn;

        // Call initialize pdof type list function
        tDofMgn.initialize_pdof_type_list( tListEqnObj );

        // Call initialize pdof type list function
        tDofMgn.create_dof_type_map();

        // Check pdof type list
        CHECK( equal_to( tDofMgn.mPdofTypeMap( 0, 0 ), 0 ) );
        CHECK( equal_to( tDofMgn.mPdofTypeMap( 1, 0 ), -1 ) );
        CHECK( equal_to( tDofMgn.mPdofTypeMap( 2, 0 ), 1 ) );
        CHECK( equal_to( tDofMgn.mPdofTypeMap( 3, 0 ), 2 ) );
    }

    TEST_CASE("Dof_Mgn_create_adofs_parallell_1","[MSI],[Dof_create_adofs_parallel_1][MSI_parallel]")
    {
        size_t tSize = par_size();

        if( tSize == 2 )
        {
            // Create node obj
            moris::uint tNodeId1 = 0;
            moris::uint tNodeId2 = 1;

            fem::Node_Base * Node1;
            fem::Node_Base * Node2;

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsId1( 2, 1 );
            Matrix< IdMat > tAdofsId2( 2, 1 );

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsInd1( 2, 1 );
            Matrix< IdMat > tAdofsInd2( 2, 1 );

            // Create generic T-matrices
            Matrix< DDRMat > tMatrix1( 2, 1 );
            Matrix< DDRMat > tMatrix2( 2, 1 );

            // Create generic adof owning processor
            Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
            Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

            // Determine process rank
            size_t tRank = par_rank();

            // Hardcode input test values
            switch( tRank )
            {
            case 0:
                tAdofsId1( 0, 0 ) = 0;
                tAdofsId1( 1, 0 ) = 5;
                tAdofsId2( 0, 0 ) = 3;
                tAdofsId2( 1, 0 ) = 0;

                tAdofsInd1( 0, 0 ) = 0;
                tAdofsInd1( 1, 0 ) = 5;
                tAdofsInd2( 0, 0 ) = 3;
                tAdofsInd2( 1, 0 ) = 0;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = -4.0;
                tMatrix2( 0, 0 ) = 2.0;
                tMatrix2( 1, 0 ) = -2.0;

                tAdofOwningProcessor1( 0, 0 ) = 0;
                tAdofOwningProcessor1( 1, 0 ) = 0;
                tAdofOwningProcessor2( 0, 0 ) = 1;
                tAdofOwningProcessor2( 1, 0 ) = 0;

                // Create generic Node Object
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );
              break;
            case 1:
                tAdofsId1( 0, 0 ) = 3;
                tAdofsId1( 1, 0 ) = 5;

                tAdofsInd1( 0, 0 ) = 0;
                tAdofsInd1( 1, 0 ) = 1;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = 3.0;

                tAdofOwningProcessor1( 0, 0 ) = 1;
                tAdofOwningProcessor1( 1, 0 ) = 0;

                // Create generic Node Object
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
              break;
            default:
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                break;
            }

            // Create dof manager and hardcode initial values
            Dof_Manager tDofMgn;

            // create a map object
            moris::map< moris::moris_id, moris::moris_index > tMap;

            tDofMgn.mNumMaxAdofs = -1;

            tDofMgn.mUseHMR = true;

            tDofMgn.mPdofTypeList.resize( 2 );
            tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;
            tDofMgn.mPdofTypeList( 1 ) = Dof_Type::UX;

            switch( tRank )
            {
            case 0:
                tDofMgn.mPdofHostList.resize( 2 );
                tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );
                tDofMgn.mPdofHostList( 1 ) = new Pdof_Host( 2, Node2 );

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 ).resize( 1 );

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 ))->mDofTypeIndex = 1;
                (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId2;

                tDofMgn.mCommTable.set_size( 2, 1, 0);
                tDofMgn.mCommTable( 1, 0 ) = 1;

                tMap[ 0 ] = 0;
                tMap[ 5 ] = 5;
                tMap[ 3 ] = 3;

              break;
            case 1:
                tDofMgn.mPdofHostList.resize( 1 );
                tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 ))->mDofTypeIndex = 1;

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds = tAdofsId1;

                tDofMgn.mCommTable.set_size( 2, 1, 1);
                tDofMgn.mCommTable( 1, 0 ) = 0;

                tMap[ 3 ] = 0;
                tMap[ 5 ] = 1;

              break;
            }

            // set map of dof manager
            tDofMgn.set_adof_map( &tMap );

            moris::Cell < Equation_Object* >tListEqnObj;
            moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
            tDofMgn.mModelSolverInterface = new Model_Solver_Interface( tMSIParameters, tListEqnObj );
            tDofMgn.mModelSolverInterface->mDofMgn = tDofMgn;
            // end hardcoding stuff

            // Create adofs and build adof lists
            tDofMgn.create_adofs();

            if ( par_rank() == 0 )
            {
                CHECK( equal_to( tDofMgn.mAdofList.size(), 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 0 )->mAdofId, 0 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 1 )->mAdofId, 4 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 2 )->mAdofId, 1 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 3 )->mAdofId, 2 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 4 )->mAdofId, 3 ) );
            }
            if ( par_rank() == 1 )
            {
                CHECK( equal_to( tDofMgn.mAdofList.size(), 4 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 0 )->mAdofId, 4 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 1 )->mAdofId, 1 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 2 )->mAdofId, 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 3 )->mAdofId, 3 ) );
            }
            delete Node1;
            delete Node2;
        }
    }

    TEST_CASE("Dof_Mgn_create_adofs_parallell_2","[MSI],[Dof_create_adofs_parallel_2][MSI_parallel]")
    {

        size_t tSize = par_size();
        if ( tSize == 2 )
        {
            // Create node obj
            moris::uint tNodeId1 = 0;
            moris::uint tNodeId2 = 1;

            fem::Node_Base * Node1;
            fem::Node_Base * Node2;

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsId1( 2, 1 );
            Matrix< IdMat > tAdofsId2( 2, 1 );

            // Create generic adofs to this nodes pdof
            Matrix< IdMat > tAdofsInd1( 2, 1 );
            Matrix< IdMat > tAdofsInd2( 2, 1 );

            // Create generic T-matrices
            Matrix< DDRMat > tMatrix1( 2, 1 );
            Matrix< DDRMat > tMatrix2( 2, 1 );

            // Create generic adof owning processor
            Matrix< IdMat > tAdofOwningProcessor1( 2, 1 );
            Matrix< IdMat > tAdofOwningProcessor2( 2, 1 );

            // Determine process rank
            size_t tRank = par_rank();

            // Hardcode input test values
            switch( tRank )
            {
            case 0:
                tAdofsId1( 0, 0 ) = 0;
                tAdofsId1( 1, 0 ) = 5;
                tAdofsId2( 0, 0 ) = 4;
                tAdofsId2( 1, 0 ) = 0;

                tAdofsInd1( 0, 0 ) = 0;
                tAdofsInd1( 1, 0 ) = 5;
                tAdofsInd2( 0, 0 ) = 4;
                tAdofsInd2( 1, 0 ) = 0;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = -4.0;
                tMatrix2( 0, 0 ) = 2.0;
                tMatrix2( 1, 0 ) = -2.0;

                tAdofOwningProcessor1( 0, 0 ) = 0;
                tAdofOwningProcessor1( 1, 0 ) = 0;
                tAdofOwningProcessor2( 0, 0 ) = 1;
                tAdofOwningProcessor2( 1, 0 ) = 0;

                // Create generic Node Object
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId2, tAdofsId2, tAdofsInd2, tMatrix2, tAdofOwningProcessor2 );
              break;
            case 1:
                tAdofsId1( 0, 0 ) = 3;
                tAdofsId1( 1, 0 ) = 5;

                tAdofsInd1( 0, 0 ) = 0;
                tAdofsInd1( 1, 0 ) = 2;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = 3.0;

                tAdofOwningProcessor1( 0, 0 ) = 1;
                tAdofOwningProcessor1( 1, 0 ) = 0;
                tAdofOwningProcessor2( 0, 0 ) = 1;
                tAdofOwningProcessor2( 1, 0 ) = 0;

                // Create generic Node Object
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
              break;
            default:
                Node1 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Proxy( tNodeId1, tAdofsId1, tAdofsInd1, tMatrix1, tAdofOwningProcessor1 );
                break;
            }

            // Create dof manager and hardcode initial values
            Dof_Manager tDofMgn;

            // create a map object
            moris::map< moris::moris_id, moris::moris_index > tMap;

            tDofMgn.mNumMaxAdofs = -1;

            tDofMgn.mUseHMR = true;

            tDofMgn.mPdofTypeList.resize( 2 );
            tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;
            tDofMgn.mPdofTypeList( 1 ) = Dof_Type::UX;

            switch( tRank )
            {
            case 0:
                tDofMgn.mPdofHostList.resize( 2 );
                tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );
                tDofMgn.mPdofHostList( 1 ) = new Pdof_Host( 2, Node2 );

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 ).resize( 1 );

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 ))->mDofTypeIndex = 1;
                (tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 1 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId2;

                tDofMgn.mCommTable.set_size( 2, 1, 0);
                tDofMgn.mCommTable( 1, 0 ) = 1;

                tMap[ 0 ] = 0;
                tMap[ 5 ] = 5;
                tMap[ 4 ] = 4;
                tMap[ 3 ] = 3;

              break;
            case 1:
                tDofMgn.mPdofHostList.resize( 1 );
                tDofMgn.mPdofHostList( 0 ) = new Pdof_Host( 2, Node1 );

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType.resize( 2 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 ).resize( 1 );
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 ).resize( 1 );

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )) = new Pdof;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )) = new Pdof;

                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 ))->mDofTypeIndex = 0;
                (tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 ))->mDofTypeIndex = 1;

                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 0 )( 0 )->mAdofIds = tAdofsId1;
                tDofMgn.mPdofHostList( 0 )->mListOfPdofTimePerType( 1 )( 0 )->mAdofIds = tAdofsId1;

                tDofMgn.mCommTable.set_size( 2, 1, 1);
                tDofMgn.mCommTable( 1, 0 ) = 0;

                tMap[ 3 ] = 0;
                tMap[ 5 ] = 2;
                tMap[ 4 ] = 1;

              break;
            }

            // set map of dof manager
            tDofMgn.set_adof_map( &tMap );

            moris::Cell < Equation_Object* >tListEqnObj;
            moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();
            tDofMgn.mModelSolverInterface = new Model_Solver_Interface( tMSIParameters, tListEqnObj );
            tDofMgn.mModelSolverInterface->mDofMgn = tDofMgn;
            // end hardcoding stuff

            // Create adofs and build adof lists
            tDofMgn.create_adofs();

            if ( par_rank() == 0 )
            {
                CHECK( equal_to( tDofMgn.mAdofList.size(), 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 0 )->mAdofId, 0 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 1 )->mAdofId, 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 2 )->mAdofId, 1 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 3 )->mAdofId, 2 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 4 )->mAdofId, 3 ) );
            }
            if ( par_rank() == 1 )
            {
                CHECK( equal_to( tDofMgn.mAdofList.size(), 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 0 )->mAdofId, 4 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 1 )->mAdofId, 5 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 2 )->mAdofId, 1 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 3 )->mAdofId, 6 ) );
                CHECK( equal_to( tDofMgn.mAdofList( 4 )->mAdofId, 3 ) );
            }
            delete Node1;
            delete Node2;
        }
    }

    }
}

