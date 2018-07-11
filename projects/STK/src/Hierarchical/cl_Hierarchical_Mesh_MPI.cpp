/*
 * cl_HIERARCHICAL_MESH_MPI.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_MPI.hpp" // STK/src/Hierarchical
using namespace moris;

// FIXME put these to functions into moris::ComTools

void Hierarchical_Mesh_MPI::broadcast_bitset_logical_and(
        BoostBitset & aBitset)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1) // call only in parallel
    {
        if( (aBitset.size() %2) == 0)
            aBitset.resize(aBitset.size()+1);
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1;
        std::string tBuffer;
        tBitset1 = &aBitset;
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0) // first proc does the calculation
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);

                // since each only deactivates its own basis, and we always refine but never coarsen
                // this bitset is a logical and
                tBitset1 = (tBitset1 & tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }

        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        aBitset = tBitset3;
    }
}

void Hierarchical_Mesh_MPI::broadcast_bitset_logical_or(
        BoostBitset & aBitset)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        if( (aBitset.size() %2) == 0)
            aBitset.resize(aBitset.size()+1);
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1;
        std::string tBuffer;
        tBitset1 = &aBitset;
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0)
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);

                // since each only refines its own basis, and we always refine but never coarsen
                // this bitset is a logical or

                tBitset1 = (tBitset1 | tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }

        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        aBitset = tBitset3;
    }
}

void Hierarchical_Mesh_MPI::gather_value_and_bcast_max(
        uint & aMessage)
{
    int tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        int aMessageInt = (int)aMessage;
        int *recvbuf;
        recvbuf = (int *)malloc(tProcSize*1*sizeof(int));
        int tMaxValuee = 0;
        if(tProcRank == 0)
        {
            MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
            tMaxValuee = *std::max_element(recvbuf, recvbuf+tProcSize);
            aMessage = (uint)tMaxValuee;
        }
        if(tProcRank > 0)
        {
            MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
        }
        broadcast(aMessage);
    }
}
