/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Pdof_Host_Test.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"

#define protected public
#define private   public
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#undef protected
#undef private

#include "fn_PRM_MSI_Parameters.hpp"

namespace moris
{
    namespace MSI
    {
    TEST_CASE("Pdof_host_set_dof_type","[MSI],[Pdof_host_set_dof_type]")
    {
        // Create node obj
        moris::uint tNodeId = 4;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsListId( 2, 1 );

        tAdofsListId( 0, 0 ) = 0;
        tAdofsListId( 1, 0 ) = 1;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsListInd( 2, 1 );

        tAdofsListInd( 0, 0 ) = 0;
        tAdofsListInd( 1, 0 ) = 1;

        // Create generic T-matrices
        Matrix< DDRMat > tMatrix( 2, 1 );

        // Create generic T-matrices
        tMatrix( 0, 0 ) = 1.0;
        tMatrix( 1, 0 ) = -2.0;

        // Create generic adof owning processor
        Matrix< IdMat > tAdofOwningProcessor( 2, 1 );

        tAdofOwningProcessor( 0, 0 ) = 0;
        tAdofOwningProcessor( 1, 0 ) = 0;

        // Create generic Node Object
        fem::Node_Base * tNode;
        tNode = new Node_Proxy( tNodeId, tAdofsListId, tAdofsListInd, tMatrix, tAdofOwningProcessor );

        // Create Pdof Host
        Pdof_Host tPdofHost( 1, tNode );

        //Check noodeId of the created pdof host
        CHECK( equal_to( tPdofHost.mNodeID, 4 ) );

        // Create pdof type enum and time step moris Mat
        enum Dof_Type tDofType = Dof_Type::TEMP;
        Matrix< DDUMat >  tTimeSteps(1, 1, 1);
        moris::uint tNumMaxPdofTypes = 1;

        Matrix< DDSMat >tDofTypeIndexMap( 4, 1, -1);
        tDofTypeIndexMap(3, 0) = 0;

        Matrix< DDUMat >tTimePerDofType( 4, 1, 1);

        //  Set pdof type and timestep
        tPdofHost.set_pdof_type( tDofType, tTimePerDofType, tNumMaxPdofTypes, tDofTypeIndexMap );

        // Check size of type and time list
        CHECK( equal_to( tPdofHost.mListOfPdofTimePerType.size(), 1 ) );
        CHECK( equal_to( tPdofHost.mListOfPdofTimePerType( 0 ).size(), 1 ) );

        // check time step indx of this pdof
        CHECK( equal_to( (tPdofHost.mListOfPdofTimePerType( 0 )( 0 ))->mTimeStepIndex, 0 ) );

        delete tNode;
    }

    TEST_CASE("Pdof_Host_Get_Adofs","[MSI],[Pdof_host_get_adofs]")
    {
        // Create node obj
        moris::uint tNodeId = 4;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsListId( 2, 1 );

        tAdofsListId( 0, 0 ) = 0;
        tAdofsListId( 1, 0 ) = 2;

        // Create generic adofs to this nodes pdof
        Matrix< IdMat > tAdofsListInd( 2, 1 );

        tAdofsListInd( 0, 0 ) = 0;
        tAdofsListInd( 1, 0 ) = 2;

        // Create generic T-matrices
        Matrix< DDRMat > tMatrix( 2, 1 );

        // Create generic T-matrices
        tMatrix( 0, 0 ) = 1.0;
        tMatrix( 1, 0 ) = -2.0;

        // Create generic adof owning processor
        Matrix< IdMat > tAdofOwningProcessor( 2, 1 );

        tAdofOwningProcessor( 0, 0 ) = 0;
        tAdofOwningProcessor( 1, 0 ) = 0;

        // Create generic Node Object
        fem::Node_Base * tNode;
        tNode = new Node_Proxy( tNodeId, tAdofsListId, tAdofsListInd, tMatrix, tAdofOwningProcessor );

        // Create Pdof Host
        Pdof_Host tPdofHost( 1, tNode );

        //Check noodeId of the created pdof host
        CHECK( equal_to( tPdofHost.mNodeID, 4 ) );

        // Create pdof type enum and time step moris Mat
        enum Dof_Type tDofType = Dof_Type::TEMP;
        Matrix< DDUMat >  tTimeSteps(1, 1, 1);
        moris::uint tNumMaxPdofTypes = 1;

        Matrix< DDSMat >tDofTypeIndexMap(4, 1, -1);
        tDofTypeIndexMap(3, 0) = 0;

        Matrix< DDUMat >tTimePerDofType(4, 1, 1);

        // Set pdof type and timestep
        tPdofHost.set_pdof_type( tDofType, tTimePerDofType, tNumMaxPdofTypes, tDofTypeIndexMap );

        // Create external adof list
        Vector< Vector< Adof * > > tAdofList;
        tAdofList.resize( 1 );
        tAdofList( 0 ).resize( 5 );

        Matrix< DDUMat > tTimeLevelOffsets( 1, 1, 0);

        moris::ParameterList tMSIParameters = prm::create_msi_parameter_list();

        Vector< Equation_Object* >tListEqnObj;
        Model_Solver_Interface  tMSI( tMSIParameters, tListEqnObj );

        Dof_Manager tDofMgn;
        tDofMgn.mPdofTypeList.resize( 1 );
        tDofMgn.mPdofTypeList( 0 ) = Dof_Type::TEMP;

        tMSI.mDofMgn = tDofMgn;

        tPdofHost.get_adofs( tTimeLevelOffsets, tAdofList, &tMSI, true );

        // Check if adofs are set to right spot
        REQUIRE( tAdofList( 0 )( 0 ) != NULL );
        REQUIRE( tAdofList( 0 )( 2 ) != NULL );
        REQUIRE( tAdofList( 0 )( 1 ) == NULL );

        delete tNode;
    }

    }
}

