/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow.hpp
 *
 */

#pragma once

#include "cl_OPT_Criteria_Interface.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Communication_Tools.hpp"

namespace xtk
{
    class Model;
}

namespace moris
{
    class Library_IO;
    //------------------------------------------------------------------------------
    namespace hmr
    {
        class HMR;
    }
    namespace mtk
    {
        class Mesh_Manager;
        class Field;
    }    // namespace mtk
    namespace ge
    {
        class GEN_Geometry_Engine;
    }
    namespace mdl
    {
        class Model;
    }
    namespace opt
    {
        class Manager;
    }

    namespace wrk
    {
        class Performer_Manager;
        //------------------------------------------------------------------------------

        class Workflow : public opt::Criteria_Interface
        {
          protected:
            wrk::Performer_Manager* mPerformerManager;

            Vector< mtk::Field* >                    Fields;
            moris::map< std::string, moris_index > tFieldNameTiIndexMap;

            uint mNumCriteria = MORIS_UINT_MAX;

            bool tIsFirstOptSolve = true;

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Workflow( wrk::Performer_Manager* aPerformerManager )
                    : mPerformerManager( aPerformerManager ){};

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Workflow(){};

        };
        //------------------------------------------------------------------------------
    }    // namespace wrk
} /* namespace moris */
