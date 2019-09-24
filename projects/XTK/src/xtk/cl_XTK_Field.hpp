/*
 * cl_XTK_Field.hpp
 *
 *  Created on: Sep 6, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"


using namespace moris;
namespace xtk
{
class Field
{
public:
    Field(){};
    Field(std::string aFieldLabel,
          moris::moris_index aFieldPhase);

    std::string                  mFieldLabel;
    moris::moris_index           mFieldPhase;
    moris::Matrix<moris::DDRMat> mFieldData;   /*Structure Node (0), Cell(1)*/
};
}

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_FIELD_HPP_ */
