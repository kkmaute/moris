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


#include "cl_Solver_Factory.hpp"
#include "cl_Solver_Input.hpp"

#define protected public
#define private   public
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Node_Obj.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#undef protected
#undef private

namespace moris
{
    namespace MSI
    {
    TEST_CASE("MSI_Test","[MSI],[MSI_Test]")
    {
        if ( par_size() == 1 )
        {
            // Create node obj
            moris::uint tNodeId1 = 0;
            moris::uint tNodeId2 = 1;

            mtk::Vertex * Node1;
            mtk::Vertex * Node2;

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

            // Create generic adof owning processor
            moris::Mat< moris::uint> tAdofOwningProcessor1( 2, 1 );
            moris::Mat< moris::uint> tAdofOwningProcessor2( 2, 1 );

            tAdofOwningProcessor1( 0, 0 ) = 0;
            tAdofOwningProcessor1( 1, 0 ) = 0;
            tAdofOwningProcessor2( 0, 0 ) = 0;
            tAdofOwningProcessor2( 1, 0 ) = 0;

            // Create generic Node Object
            Node1 = new Node_Obj( tNodeId1, tAdofs1, tMatrix1, tAdofOwningProcessor1 );
            Node2 = new Node_Obj( tNodeId2, tAdofs2, tMatrix2, tAdofOwningProcessor2 );

            moris::uint tNumEquationObjects = 2;

            moris::uint tNumNodes = 2;

            moris::Cell < Equation_Object* >tListEqnObj( tNumEquationObjects, nullptr );

            // Create List with node pointern correponding to generic equation object
            moris::Cell< mtk::Vertex* > tNodeIds_1( tNumNodes );
            tNodeIds_1( 0 ) = Node1;
            tNodeIds_1( 1 ) = Node2;

            moris::Cell< mtk::Vertex* > tNodeIds_2( tNumNodes );
            tNodeIds_2( 0 ) = Node1;
            tNodeIds_2( 1 ) = Node2;

            // Create generic equation objects
            Equation_Object EquObj_1( tNodeIds_1 );
            Equation_Object EquObj_2( tNodeIds_2 );

            EquObj_1.mEqnObjDofTypeList.resize( 1, Dof_Type::TEMP);
            EquObj_2.mEqnObjDofTypeList.resize( 1, Dof_Type::TEMP);

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

            moris::Mat< moris::uint > tCommTable( 1, 1, 0 );

            Model_Solver_Interface tMSI( tListEqnObj, tCommTable );

            // create solver input object
            moris::MSI::MSI_Solver_Interface *  tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( &tMSI, tMSI.get_dof_manager() );

            // create solver factory
            moris::Solver_Factory  tSolFactory;

            // create solver object
            std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

            tLin->solve_linear_system();

            moris::Mat< moris::real > tSolution;
            tLin->get_solution( tSolution );

            CHECK( equal_to( tSolution( 0, 0 ), -2 ) );
            CHECK( equal_to( tSolution( 1, 0 ), 5 ) );

            delete Node1;
            delete Node2;
            delete tSolverInput;
        }
    }

    TEST_CASE("MSI_Test_parallel","[MSI],[MSI_Test_parallel][MSI_parallel]")
    {
        size_t tSize = par_size();

        if ( tSize == 2 )
        {
            // Create node obj
            moris::uint tNodeId1 = 0;
            moris::uint tNodeId2 = 1;
            moris::uint tNodeId3 = 2;
            moris::uint tNodeId4 = 3;

            mtk::Vertex * Node1;
            mtk::Vertex * Node2;

            // Create generic adofs to this nodes pdof
            moris::Mat< moris::sint> tAdofs1( 2, 1 );
            moris::Mat< moris::sint> tAdofs2( 2, 1 );

            // Create generic T-matrices
            moris::Mat< moris::real> tMatrix1( 2, 1 );
            moris::Mat< moris::real> tMatrix2( 2, 1 );

            // Create generic adof owning processor
            moris::Mat< moris::uint> tAdofOwningProcessor1( 2, 1 );
            moris::Mat< moris::uint> tAdofOwningProcessor2( 2, 1 );

            // Determine process rank
            size_t tRank = par_rank();
            size_t tSize = par_size();

            moris::Mat< moris::uint > tCommTable( 2, 1 );
            moris::uint tNumEquationObjects;
            moris::uint tNumNodes;
            moris::Cell < Equation_Object* >tListEqnObj;
            moris::Cell< mtk::Vertex* > tNodeIds_1;
            moris::Cell< mtk::Vertex* > tNodeIds_2;

            // Hardcode input test values
            switch( tRank )
            {
            case 0:
                tAdofs1( 0, 0 ) = 0;
                tAdofs1( 1, 0 ) = 1;
                tAdofs2( 0, 0 ) = 1;
                tAdofs2( 1, 0 ) = 2;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = 1.0;
                tMatrix2( 0, 0 ) = 1.0;
                tMatrix2( 1, 0 ) = 2.0;

                tAdofOwningProcessor1( 0, 0 ) = 0;
                tAdofOwningProcessor1( 1, 0 ) = 0;
                tAdofOwningProcessor2( 0, 0 ) = 0;
                tAdofOwningProcessor2( 1, 0 ) = 1;

                // Create generic Node Object
                Node1 = new Node_Obj( tNodeId1, tAdofs1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Obj( tNodeId2, tAdofs2, tMatrix2, tAdofOwningProcessor2 );

                tCommTable( 0, 0 ) = 0;
                tCommTable( 1, 0 ) = 1;

                tNumEquationObjects = 2;
                tNumNodes = 2;
                tListEqnObj.resize( tNumEquationObjects, nullptr );

                // Create List with node pointern correponding to generic equation object
                tNodeIds_1.resize( tNumNodes );
                tNodeIds_1( 0 ) = Node1;
                tNodeIds_1( 1 ) = Node2;

                tNodeIds_2.resize( tNumNodes );
                tNodeIds_2( 0 ) = Node1;
                tNodeIds_2( 1 ) = Node2;

              break;
            case 1:
                tAdofs1( 0, 0 ) = 3;
                tAdofs1( 1, 0 ) = 1;
                tAdofs2( 0, 0 ) = 3;
                tAdofs2( 1, 0 ) = 0;

                tMatrix1( 0, 0 ) = 1.0;
                tMatrix1( 1, 0 ) = 1.0;
                tMatrix2( 0, 0 ) = 1.0;
                tMatrix2( 1, 0 ) = 2.0;

                tAdofOwningProcessor1( 0, 0 ) = 1;
                tAdofOwningProcessor1( 1, 0 ) = 0;
                tAdofOwningProcessor2( 0, 0 ) = 1;
                tAdofOwningProcessor2( 1, 0 ) = 0;

                // Create generic Node Object
                Node1 = new Node_Obj( tNodeId3, tAdofs1, tMatrix1, tAdofOwningProcessor1 );
                Node2 = new Node_Obj( tNodeId4, tAdofs2, tMatrix2, tAdofOwningProcessor2 );

                tCommTable( 0, 0 ) = 1;
                tCommTable( 1, 0 ) = 0;

                tNumEquationObjects = 2;
                tNumNodes = 2;
                tListEqnObj.resize( tNumEquationObjects, nullptr );

                // Create List with node pointern correponding to generic equation object
                tNodeIds_1.resize( tNumNodes );
                tNodeIds_1( 0 ) = Node1;
                tNodeIds_1( 1 ) = Node2;

                tNodeIds_2.resize( tNumNodes );
                tNodeIds_2( 0 ) = Node1;
                tNodeIds_2( 1 ) = Node2;
              break;
            }

            // Create generic equation objects
            Equation_Object EquObj_1( tNodeIds_1 );
            Equation_Object EquObj_2( tNodeIds_2 );

            EquObj_1.mEqnObjDofTypeList.resize( 1, Dof_Type::TEMP );
            EquObj_2.mEqnObjDofTypeList.resize( 1, Dof_Type::TEMP );

            EquObj_1.mJacobian.set_size( 2, 2, 0.0);
            EquObj_2.mJacobian.set_size( 2, 2, 0.0);
            EquObj_1.mResidual.set_size( 2, 1, 0.0);
            EquObj_2.mResidual.set_size( 2, 1, 0.0);

            EquObj_1.mJacobian( 0, 0 ) = 1;
            EquObj_1.mJacobian( 0, 1 ) = 2;
            EquObj_2.mJacobian( 1, 0 ) = 1;
            EquObj_2.mJacobian( 1, 1 ) = -3;

            EquObj_1.mResidual( 0, 0 ) = 2;

            EquObj_2.mResidual( 1, 0 ) =2;

            // Create List with equation objects
            tListEqnObj( 0 ) = & EquObj_1;
            tListEqnObj( 1 ) = & EquObj_2;

            Model_Solver_Interface tMSI( tListEqnObj, tCommTable );

            // create solver input object
            moris::MSI::MSI_Solver_Interface * tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( &tMSI, tMSI.get_dof_manager() );

            // create solver factory
            moris::Solver_Factory  tSolFactory;

            // create solver object
            std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput );

            tLin->solve_linear_system();

            moris::Mat< moris::real > tSolution;
            tLin->get_solution( tSolution );

            if ( par_rank() == 0 )
            {
                CHECK( equal_to( tSolution( 0, 0 ), 0 ) );
                CHECK( equal_to( tSolution( 1, 0 ), 2 ) );
            }
            else if ( par_rank() == 1 )
            {
                CHECK( equal_to( tSolution( 0, 0 ), -1 ) );
                CHECK( equal_to( tSolution( 1, 0 ), 0 ) );
            }
            delete Node1;
            delete Node2;
        }
    }
    }
}


