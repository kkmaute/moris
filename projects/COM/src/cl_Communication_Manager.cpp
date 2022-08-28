/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Communication_Manager.cpp
 *
 */

#include "cl_Communication_Manager.hpp"

namespace moris
{
    Comm_Manager::Comm_Manager(
            int *argc,
            char ***argv)
    {
        MPI_Init(argc, argv);

        // Set MPI_COMM_WORLD as the default communicator
        mComm = Cell<MPI_Comm>({MPI_COMM_WORLD});

        // Give it a global name
        std::string tGlobCommName = "moris_glob";
        mCommName = Cell<std::string>({tGlobCommName});

        // save path to moris executable
        mMorisExec = std::string( *argv[ 0 ] );
    }

    //--------------------------------------------------------------------------------

    Comm_Manager::Comm_Manager(MPI_Comm & aComm)
    {
        // Set your communicator to the mComm as global
        mComm = Cell<MPI_Comm>({aComm});

        // Give it a name
        std::string tGlobCommName = "moris_glob";

        // Set name
        mCommName = Cell<std::string>({tGlobCommName});
    }

    //--------------------------------------------------------------------------------

    Comm_Manager::~Comm_Manager()
    {
    }

    //--------------------------------------------------------------------------------

    void
    Comm_Manager::initialize(
            int *argc,
            char ***argv)
    {
        MPI_Init(argc, argv);

        // Set MPI_COMM_WORLD as the default communicator
        mComm = Cell<MPI_Comm>({MPI_COMM_WORLD});

        // Give it a global name
        std::string tGlobCommName = "moris_glob";
        mCommName = Cell<std::string>({tGlobCommName});

        // save path to moris executable
        mMorisExec = std::string( *argv[ 0 ] );
    }

    //--------------------------------------------------------------------------------

    MPI_Comm
    Comm_Manager::get_global_comm()
    {
        return mComm(0);
    }

    //--------------------------------------------------------------------------------

    MPI_Comm
    Comm_Manager::get_comm(size_t aCommIndex)
    {
        MORIS_ASSERT( aCommIndex<mComm.size()-1,
                "Communicator Index out of bounds, did you add the communicator?" );

        return mComm( aCommIndex );
    }

    //--------------------------------------------------------------------------------

    MPI_Comm
    Comm_Manager::get_comm()
    {
        return mComm( mActiveCommunicator );
    }

    //--------------------------------------------------------------------------------

    size_t
    Comm_Manager::add_communicator(
            MPI_Comm          & aNewComm,
            const std::string & aCommName)
    {
        // Use the current size as the next index
        size_t aCommIndex = mComm.size();
        mActiveCommunicator = aCommIndex;

        // Append the communicator to the list
        mComm.push_back( aNewComm );
        mCommName.push_back( aCommName );

        // return the index
        return aCommIndex;
    }

    //--------------------------------------------------------------------------------

    void Comm_Manager::remove_communicator(size_t aCommIndex)
    {
        MORIS_ERROR((aCommIndex < mComm.size()), "Tried to remove a communicator index that doesn't exist.");
        mComm.erase(aCommIndex);
    }

    //--------------------------------------------------------------------------------

    void
    Comm_Manager::finalize()
    {
        MPI_Finalize();
    }

    //--------------------------------------------------------------------------------

    const std::string &
    Comm_Manager::get_exec_path()
    {
        return mMorisExec;
    }
}

