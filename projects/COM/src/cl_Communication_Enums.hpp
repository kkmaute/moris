/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Communication_Enums.hpp
 *
 */

#ifndef SRC_COMM_CL_COMMUNCATION_ENUMS_HPP_
#define SRC_COMM_CL_COMMUNCATION_ENUMS_HPP_

enum class CommunicationType
{
    BCAST,           // Broadcast enum
    GATHER,          // Gather enum
    SCATTER,         // Scatter enum
    RECV,            // Blocking receive enum
    SEND,            // Blocking send enum
    IRECV,           // Nonblocking receive enum
    ISEND,           // Nonblocking send enum
    ID_ALLOCATE,     // ENTITY ALLOCATION
    END_ENUM         // This should be the last enum (used to check the length of the in some functions)
};

#endif /* SRC_COMM_CL_COMMUNCATION_ENUMS_HPP_ */

