/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Communicator_Epetra.cpp
 *
 */

#include "cl_Communicator_Epetra.hpp"

using namespace moris;

 Communicator_Epetra::Communicator_Epetra() : mEpetraComm(NULL)
{
     //Fixme replace after next commit
     MPI_Comm world=MPI_COMM_WORLD;
//     MPI_Comm tSolv;
//
//     int tRank;
//     //sint tColor=2;
//
//     MPI_Comm_rank(world, &tRank);
//     MPI_Comm_split(world, tRank%1, 0, &tSolv);
     //MPI_Comm_split(world, tColor, 0, &tSolv);

#ifdef MORIS_HAVE_PARALLEL
    mEpetraComm = new Epetra_MpiComm(world);
#else
    mEpetraComm = new Epetra_SerialComm();
#endif
}
// ----------------------------------------------------------------------------

Communicator_Epetra::~Communicator_Epetra()
{
    delete( mEpetraComm );
}

