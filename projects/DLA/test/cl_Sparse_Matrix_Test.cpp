/*
 * cl_Sparse_Matrix_Test.cpp
 *
 *  Created on: Mar 19, 2018
 *      Author: schmidt
 */

#ifdef MORIS_HAVE_PARALLEL
 #include "Epetra_MpiComm.h"
 #include <mpi.h>
#endif

#include "catch.hpp"

#include "fn_equal_to.hpp" // ALG/src

#include "typedefs.hpp" // COR/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp" // COM/src/
#include "cl_Matrix_Vector_Factory.hpp" // DLA/src/
#include "cl_Solver_Input_Test.hpp" // DLA/src/
#include "cl_Vector.hpp" // DLA/src/
#include "cl_Sparse_Matrix.hpp" // DLA/src/

namespace moris
{
TEST_CASE("Sparse Mat","[Sparse Mat],[DistLinAlg]")
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    if (size == 4)
    {
    // Build Input Class
    Solver_Input* tSolverInput = new Solver_Input_Test( );

    // Build matrix factory
    Matrix_Vector_Factory      tMatFactory;

    // Build map
    Map_Class * tLocalMap = tMatFactory.create_map( tSolverInput->get_num_my_dofs(),
                                                    tSolverInput->get_my_local_global_map(),
                                                    tSolverInput->get_constr_dof(),
                                                    tSolverInput->get_my_local_global_map() );

    // Create pointer to sparse matrix
    Sparse_Matrix * tMat = tMatFactory.create_matrix( tSolverInput, tLocalMap );

    // Build sparse matrix graph
    for ( moris::uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        tMat->build_graph( tElementTopology.length(),  tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    // Fill element matrices into global matrix
    for ( uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        tSolverInput->get_element_matrix( Ii, tElementMatrix );

        tMat->fill_matrix( tElementTopology.length(), tElementMatrix, tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    //tMat->print_matrix_to_screen();

    // Set up output matrix
    sint tGlobalRow = 8;    sint tLength = 13;    sint tNumEntries = 5;
    moris::Matrix< DDRMat > tValues ( tLength, 1, 0.0 );

    // Get matrix values
    sint err = tMat->get_matrix()->ExtractGlobalRowCopy( tGlobalRow, tLength, tNumEntries, tValues.data() );

    // Compare to true values.
    if ( rank == 0 )
    {
        CHECK( equal_to( tValues( 0, 0 ), 24) );
        CHECK( equal_to( tValues( 4, 0 ), -6) );
        CHECK( equal_to( tValues( 8, 0 ), -3) );
    }
    delete ( tSolverInput );
    delete ( tLocalMap );
    delete ( tMat );

    }
}

TEST_CASE("Scale Sparse Mat","[Scale Sparse Mat],[DistLinAlg]")
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    if (size == 4)
    {
    // Build Input Class
    Solver_Input* tSolverInput = new Solver_Input_Test( );

    // Build matrix factory
    Matrix_Vector_Factory      tMatFactory;

    // Build map
    Map_Class * tMap = tMatFactory.create_map( tSolverInput->get_num_my_dofs(),
                                               tSolverInput->get_my_local_global_map(),
                                               tSolverInput->get_constr_dof(),
                                               tSolverInput->get_my_local_global_map() );

    // build distributed vector
    Dist_Vector * tVectorScale = tMatFactory.create_vector( tSolverInput, tMap, VectorType::FREE );

    // Create pointer to sparse matrix
    Sparse_Matrix * tMat = tMatFactory.create_matrix( tSolverInput, tMap );

    // Build sparse matrix graph
    for ( moris::uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        tMat->build_graph( tElementTopology.length(),  tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    // Fill element matrices into global matrix
    for ( uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        tSolverInput->get_element_matrix( Ii, tElementMatrix );

        tMat->fill_matrix( tElementTopology.length(), tElementMatrix, tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    //tMat->print_matrix_to_screen();

    // Put Values in scaling vector
    tVectorScale->vec_put_scalar( 2.25 );

    tMat->sparse_mat_left_scale( *tVectorScale );

    // Set up output matrix
    moris::sint tGlobalRow = 8;    moris::sint tLength = 13;    moris::sint tNumEntries = 5;
    moris::Matrix< DDRMat > tValues (tLength, 1, 0.0);

    // Get matrix values
    sint err = tMat->get_matrix()->ExtractGlobalRowCopy(tGlobalRow, tLength, tNumEntries, tValues.data() );

    // Compare to true values.
    if (rank == 0)
    {
        CHECK(equal_to(tValues(0,0), 54));
        CHECK(equal_to(tValues(4,0), -13.5));
        CHECK(equal_to(tValues(8,0), -6.75));
    }
    delete ( tSolverInput );
    delete ( tMap );
    delete ( tVectorScale );
    delete ( tMat );
    }
}

TEST_CASE("Diagonal Sparse Mat","[Diagonal Sparse Mat],[DistLinAlg]")
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    if (size == 4)
    {
    // Build Input Class
    Solver_Input* tSolverInput = new Solver_Input_Test( );

    // Build matrix factory
    Matrix_Vector_Factory      tMatFactory;

    // Build map
    Map_Class * tMap = tMatFactory.create_map( tSolverInput->get_num_my_dofs(),
                                               tSolverInput->get_my_local_global_map(),
                                               tSolverInput->get_constr_dof(),
                                               tSolverInput->get_my_local_global_map() );

    // build distributed vector
    Dist_Vector * tVectorDiagonal = tMatFactory.create_vector( tSolverInput, tMap, VectorType::FREE );

    // Create pointer to sparse matrix
    Sparse_Matrix * tMat = tMatFactory.create_matrix( tSolverInput, tMap );

    // Build sparse matrix graph
    for ( moris::uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        tMat->build_graph( tElementTopology.length(),  tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    // Fill element matrices into global matrix
    for ( uint Ii=0; Ii< tSolverInput->get_num_my_elements(); Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        tSolverInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        tSolverInput->get_element_matrix( Ii, tElementMatrix );

        tMat->fill_matrix( tElementTopology.length(), tElementMatrix, tElementTopology );
    }

    // Call Global Asemby to ship information between processes
    tMat->matrix_global_asembly();

    //tMat->print_matrix_to_screen();

    // extract diagonal value into Dist Vector
    tMat->get_diagonal( *tVectorDiagonal );

    moris::Matrix< DDRMat > tDiagonal ( 15, 1, 0.0 );

    // needed as offset parameter for Epetra. =0
    sint tMyLDA = 0;

    // Get solution and output it in moris::Mat LHSValues
    tVectorDiagonal->get_vector()->ExtractCopy( tDiagonal.data(), tMyLDA );

    // Compare to true values.
    if (rank == 0)
    {
        CHECK(equal_to(tDiagonal(0,0), 24));
        CHECK(equal_to(tDiagonal(1,0), 24));
        CHECK(equal_to(tDiagonal(2,0), 48));
        CHECK(equal_to(tDiagonal(3,0), 48));
        CHECK(equal_to(tDiagonal(4,0), 24));
        CHECK(equal_to(tDiagonal(5,0), 24));
    }
    if (rank == 3)
    {
        CHECK(equal_to(tDiagonal(0,0), 24));
        CHECK(equal_to(tDiagonal(1,0), 24));
        CHECK(equal_to(tDiagonal(2,0), 12));
        CHECK(equal_to(tDiagonal(3,0), 12));
    }

    // Set all entries of vector to 100.33
    tVectorDiagonal->vec_put_scalar( 100.33 );

    // Replace Matrix diagonal values with tVectorDiagonal values
    tMat->replace_diagonal_values( *tVectorDiagonal );

    // Set up output matrix
    moris::sint tGlobalRow = 8;    moris::sint tLength = 13;    moris::sint tNumEntries = 5;
    moris::Matrix< DDRMat > tValues (tLength, 1, 0.0);

    // Get matrix values
    sint err = tMat->get_matrix()->ExtractGlobalRowCopy(tGlobalRow, tLength, tNumEntries, tValues.data() );

    // Compare to true values.
    if (rank == 0)
    {
        CHECK(equal_to(tValues(0,0), 100.33));
        CHECK(equal_to(tValues(4,0), -6));
        CHECK(equal_to(tValues(8,0), -3));
    }
    delete ( tSolverInput );
    delete ( tMap );
    delete ( tVectorDiagonal );
    delete ( tMat );
    }
}
}
