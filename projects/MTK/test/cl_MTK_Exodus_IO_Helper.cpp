/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Exodus_IO_Helper.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Exodus_IO_Helper.hpp"

#include "paths.hpp"

namespace moris
{
    namespace mtk
    {
        TEST_CASE("Parallel Exodus IO Helper","[EXO_IO_2P]")
        {
            if(par_size() == 2)
            {
                // TODO: use an existing exo file
                std::string tPrefix = moris::get_base_moris_dir();
                std::string tFileOutput = tPrefix + "projects/MTK/test/Test_Files/mtk_2_proc_test.exo";

                // File with element cmap
                std::string tFileOutputwElemCmap = tPrefix + "projects/MTK/test/Test_Files/mtk_2_proc_test_elem_cmap.exo";

                // initialize exodus io help
                Exodus_IO_Helper tExoIO(tFileOutput);

                // element cmap data;
                if(par_rank() == 0)
                {
                    Matrix<IdMat> tElementIdsOnBoundary ={{1}};
                    Matrix<IdMat> tSideOrdinalOnBoundary ={{5}};
                    Matrix<IdMat> tSideSharedProc ={{1}};
                    // Create a new exodus file with element cmaps
                    tExoIO.create_new_exo_with_elem_cmaps_from_existing_exo(tFileOutputwElemCmap,tElementIdsOnBoundary,tSideOrdinalOnBoundary,tSideSharedProc);
                }
                else if(par_rank() == 1)
                {
                    Matrix<IdMat> tElementIdsOnBoundary ={{1}};
                    Matrix<IdMat> tSideOrdinalOnBoundary ={{5}};
                    Matrix<IdMat> tSideSharedProc ={{0}};
                    // Create a new exodus file with element cmaps
                    tExoIO.create_new_exo_with_elem_cmaps_from_existing_exo(tFileOutputwElemCmap,tElementIdsOnBoundary,tSideOrdinalOnBoundary,tSideSharedProc);
                }
            }
        }
    }
}

