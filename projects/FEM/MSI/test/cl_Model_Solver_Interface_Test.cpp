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

        tMSI.solve_system();
    }

    }
}


