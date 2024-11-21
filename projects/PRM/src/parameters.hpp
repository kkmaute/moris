/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * parameters.hpp
 *
 */

#pragma once

// include all parameter function headers defined in this directory
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_MIG_Parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_STK_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_WRK_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_Module_Parameter_Lists.hpp"

// Typedefs for ease of typing
namespace moris
{
    typedef OPT_Submodule OPT;
    typedef HMR_Submodule HMR;
    typedef GEN_Submodule GEN;
    typedef FEM_Submodule FEM;
    typedef SOL_Submodule SOL;
}
