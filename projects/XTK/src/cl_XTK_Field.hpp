/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Field.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

using namespace moris;
namespace moris::xtk
{
    class Field
    {
      public:
        Field(){};
        Field(
                std::string        aFieldLabel,
                moris::moris_index aFieldPhase );

        std::string        mFieldLabel;
        moris::moris_index mFieldPhase;
        Matrix< DDRMat >   mFieldData; /*Structure Node (0), Cell(1)*/
    };                                 // end: class xtk::Field
}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_ */
