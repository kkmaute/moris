/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner.cpp
 *
 */

#include "cl_DLA_Preconditioner.hpp"
using namespace moris::dla;


//-------------------------------------------------------------------------------

Preconditioner::Preconditioner()
{
}

//-------------------------------------------------------------------------------

Preconditioner::Preconditioner(
        moris::Parameter_List* aParameterlist,
        Linear_Problem*       aLinearSystem )
{
    this->initialize(
            aParameterlist,
            aLinearSystem );
}

//-------------------------------------------------------------------------------

Preconditioner::~Preconditioner()
{
}

//-------------------------------------------------------------------------------



