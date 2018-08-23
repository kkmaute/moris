/*
 * cl_Dof_Manager.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "cl_MSI_Model_Solver_Interface.hpp"

//namespace moris
//{
//    namespace MSI
//    {
//        Model_Solver_Interface::Model_Solver_Interface( moris::Linear_Solver * aLin,
//                                                        moris::Solver_Input  * aInput,
//                                                               Sparse_Matrix * aMat,
//                                                        moris::Dist_Vector   * aVectorRHS )
//        {
//            // Get local number of elements
//            moris::uint numLocElements = aInput->get_num_my_elements();
//
//            // Loop over all local elements to build matrix graph
//            for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
//            {
//                Mat< int > tElementTopology;
//                aInput->get_element_topology(Ii, tElementTopology );
//
//                aMat->build_graph( tElementTopology.length(), tElementTopology );
//            }
//            // global assembly to switch entries to the right proceccor
//            aMat->matrix_global_asembly();
//
//            // Loop over all local elements to fill matrix and RHS
//            for (moris::uint Ii=0; Ii< numLocElements; Ii++)
//            {
//                moris::Mat< int > tElementTopology;
//                aInput->get_element_topology( Ii, tElementTopology );
//
//                moris::Mat< moris::real > tElementMatrix;
//                aInput->get_element_matrix( Ii, tElementMatrix );
//
//                moris::Mat< moris::real > tElementRHS;
//                aInput->get_element_rhs( Ii, tElementRHS );
//
//                // Fill element in distributed matrix
//                aMat->fill_matrix( tElementTopology.length(),
//                                         tElementMatrix,
//                                         tElementTopology );
//
//                // Fill elementRHS in distributed RHS
//                aVectorRHS->sum_into_global_values( tElementTopology.length(),
//                                                    tElementTopology,
//                                                    tElementRHS );
//            }
//            // global assembly to switch entries to the right proceccor
//            aVectorRHS->vector_global_asembly();
//            aMat->matrix_global_asembly();
//
//            // build linear system on solver class
//            aLin->build_linear_system();
//        }
//    }
//}
