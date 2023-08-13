/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_mult.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_OP_ELEMWISE_MULT_HPP_
#define PROJECTS_LINALG_SRC_OP_ELEMWISE_MULT_HPP_

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_elemwise_mult_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_elemwise_mult_Arma.hpp"
#endif

//FIXME: Christian, All functions need a description and function call with two moris matrix arguments
#endif /* PROJECTS_LINALG_SRC_OP_ELEMWISE_MULT_HPP_ */

