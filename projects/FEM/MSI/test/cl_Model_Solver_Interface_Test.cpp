/*
 * cl_Dof_Manager_Test.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifdef MORIS_HAVE_PARALLEL
 #include "Epetra_MpiComm.h"
 #include <mpi.h>
#endif

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "typedefs.hpp"
#include "cl_Mat.hpp"
#include "cl_Communication_Tools.hpp"

#define protected public
#define private   public
#include "cl_Equation_Object.hpp"
#include "cl_Node_Obj.hpp"
#include "cl_Model_Solver_Interface.hpp"
#include "cl_Dof_Manager.hpp"
#include "cl_Pdof_Host.hpp"
#undef protected
#undef private

namespace moris
{
    namespace MSI
    {
    TEST_CASE("MSI_Test","[MSI],[MSI_Test]")
    {
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 1;

        Node_Obj * Node1;
        Node_Obj * Node2;

        // Create generic adofs to this nodes pdof
        moris::Mat< moris::sint> tAdofs1( 2, 1 );
        moris::Mat< moris::sint> tAdofs2( 2, 1 );

        tAdofs1( 0, 0 ) = 0;
        tAdofs1( 1, 0 ) = 1;
        tAdofs2( 0, 0 ) = 1;
        tAdofs2( 1, 0 ) = 0;

        // Create generic T-matrices
        moris::Mat< moris::real> tMatrix1( 2, 1 );
        moris::Mat< moris::real> tMatrix2( 2, 1 );

        // Create generic T-matrices
        tMatrix1( 0, 0 ) = 1.0;
        tMatrix1( 1, 0 ) = 1.0;
        tMatrix2( 0, 0 ) = 1.0;
        tMatrix2( 1, 0 ) = 2.0;

        // Create generic Node Object
        Node1 = new Node_Obj( tNodeId1, tAdofs1, tMatrix1 );
        Node2 = new Node_Obj( tNodeId2, tAdofs2, tMatrix2 );

        moris::uint tNumEquationObjects = 2;

        moris::uint tNumNodes = 2;

        moris::Cell < Equation_Object* >tListEqnObj;
        tListEqnObj.resize( tNumEquationObjects, nullptr );

        // Create List with node pointern correponding to generic equation object
        moris::Cell< Node_Obj* > tNodeIds_1( tNumNodes );
        tNodeIds_1( 0 ) = Node1;
        tNodeIds_1( 1 ) = Node2;

        moris::Cell< Node_Obj* > tNodeIds_2( tNumNodes );
        tNodeIds_2( 0 ) = Node1;
        tNodeIds_2( 1 ) = Node2;

        // Create generic equation objects
        Equation_Object EquObj_1( tNodeIds_1 );
        Equation_Object EquObj_2( tNodeIds_2 );

        EquObj_1.mJacobian.set_size( 2, 2, 0.0);
        EquObj_2.mJacobian.set_size( 2, 2, 0.0);
        EquObj_1.mResidual.set_size( 2, 1, 0.0);
        EquObj_2.mResidual.set_size( 2, 1, 0.0);

        EquObj_1.mJacobian( 0, 0 ) = 1;
        EquObj_1.mJacobian( 0, 1 ) = 2;
        EquObj_2.mJacobian( 1, 0 ) = 1;
        EquObj_2.mJacobian( 1, 1 ) = -3;

        EquObj_1.mResidual( 0, 0 ) = 5;

        // Create List with equation objects
        tListEqnObj( 0 ) = & EquObj_1;
        tListEqnObj( 1 ) = & EquObj_2;

        Model_Solver_Interface tMSI( tNumEquationObjects, tListEqnObj );

        moris::Mat< moris::real > tSolution;
        tMSI.solve_system( tSolution );

        CHECK( equal_to( tSolution( 0, 0 ), 3 ) );
        CHECK( equal_to( tSolution( 1, 0 ), 1 ) );
    }

    TEST_CASE("MSI_Test_parallel","[MSI],[MSI_Test_parallel]")
    {

   // Determine process rank
       size_t rank = par_rank();
       size_t size = par_size();

       if (size == 4)
       {
       // Set input integer and pointer
       uint tNumMyDofs = 0;
       Mat < int > tMyGlobalElements;
       Mat < uint > tMyConstraintDofs;

       // Define input test values
       switch( rank )
           {
           case 0:
             tNumMyDofs = 8;
             tMyGlobalElements.resize( tNumMyDofs, 1 );
             tMyConstraintDofs.resize( 2, 1 );
             tMyGlobalElements(0,0) = 0;    tMyGlobalElements(1,0) = 1;  tMyGlobalElements(2,0) = 8;    tMyGlobalElements(3,0) = 9;    tMyGlobalElements(4,0) = 16;    tMyGlobalElements(5,0) = 17;    tMyGlobalElements(6,0) = 14;    tMyGlobalElements(7,0) = 15;
             tMyConstraintDofs(0,0) = 0;    tMyConstraintDofs(1,0) = 1;
             break;
           case 1:
             tNumMyDofs = 4;
             tMyGlobalElements.resize( tNumMyDofs, 1 );
             tMyConstraintDofs.resize( 1, 1 );
             tMyGlobalElements(0,0) = 2;   tMyGlobalElements(1,0) = 3;    tMyGlobalElements(2,0) = 10;    tMyGlobalElements(3,0) = 11;
             tMyConstraintDofs(0,0) = 3;
             break;
           case 2:
             tNumMyDofs = 2;
             tMyGlobalElements.resize( tNumMyDofs, 1 );
             tMyGlobalElements(0,0) = 4;    tMyGlobalElements(1,0) = 5;
             break;
           case 3:
             tNumMyDofs = 4;
             tMyGlobalElements.resize( tNumMyDofs, 1 );
             tMyGlobalElements(0,0) = 12;    tMyGlobalElements(1,0) = 13;    tMyGlobalElements(2,0) = 6;    tMyGlobalElements(3,0) = 7;
             break;
            }
        // Create node obj
        moris::uint tNodeId1 = 0;
        moris::uint tNodeId2 = 1;

        Node_Obj * Node1;
        Node_Obj * Node2;

        // Create generic adofs to this nodes pdof
        moris::Mat< moris::sint> tAdofs1( 2, 1 );
        moris::Mat< moris::sint> tAdofs2( 2, 1 );



    }

    }
    }
}


