/*
 * cl_XTK_Field.cpp
 *
 *  Created on: Sep 6, 2019
 *      Author: doble
 */

#include "cl_XTK_Field.hpp"

namespace xtk
{
Field::Field(std::string aFieldLabel,
             moris::moris_index aFieldPhase):
                     mFieldLabel(aFieldLabel),
                     mFieldPhase(aFieldPhase),
                     mFieldData(0,0)
{

}

}
