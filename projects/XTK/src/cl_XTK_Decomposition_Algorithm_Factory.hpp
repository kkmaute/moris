/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Algorithm_Factory.hpp
 *
 */

#ifndef SRC_cl_XTK_Decomposition_Algorithm_Factory
#define SRC_cl_XTK_Decomposition_Algorithm_Factory

#include "cl_XTK_Enums.hpp"
#include <memory>
#include "cl_Param_List.hpp"

namespace xtk
{
class Decomposition_Algorithm;

std::shared_ptr<Decomposition_Algorithm>
create_decomposition_algorithm(
  enum Subdivision_Method aSubdivisionMethod,
  moris::ParameterList&   aParameterList);

}// namespace xtk

#endif /* cl_XTK_Decomposition_Algorithm_Factory.hpp */
