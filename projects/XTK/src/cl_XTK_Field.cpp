/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Field.cpp
 *
 */

#include "cl_XTK_Field.hpp"

namespace xtk
{
    Field::Field( std::string  aFieldLabel,
            moris::moris_index aFieldPhase )
            : mFieldLabel( aFieldLabel )
            , mFieldPhase( aFieldPhase )
            , mFieldData( 0, 0 )
    {
        // do nothing besides initializing the member data
    }

}    // namespace xtk
