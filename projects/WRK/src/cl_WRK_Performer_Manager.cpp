/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Performer_Manager.cpp
 *
 */

#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_HMR.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_WRK_Performer_Manager.hpp"

#include <utility>

#include "cl_Library_IO.hpp"

namespace moris::wrk
{

    Performer_Manager::Performer_Manager( std::shared_ptr< Library_IO > aLibrary )
            : mLibrary( std::move( aLibrary ) )
    {
    }

    //------------------------------------------------------------------------------

    Performer_Manager::~Performer_Manager()
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris::wrk
