/*
 * cl_MTK_Exodus_IO_Helper.cpp
 *
 *  Created on: Nov 28, 2018
 *      Author: doble
 */
#include "catch.hpp"

#include "cl_MTK_Exodus_IO_Helper.hpp"

namespace moris
{
namespace mtk
{

TEST_CASE("Parallel Exodus IO Helper","[EXO_IO_2P]")
                {
    if(par_size() == 2)
    {
        // TODO: use an existing exo file
        std::string tPrefix = std::getenv("MORISROOT");
        std::string tFileOutput = tPrefix + "projects/MTK/test/Test_Files/mtk_2_proc_test.exo";

        // File with element cmap
        std::string tFileOutputwElemCmap = tPrefix + "projects/MTK/test/Test_Files/mtk_2_proc_test_elem_cmap.exo";

        // initialize exodus io help
        Exodus_IO_Helper tExoIO(tFileOutput.c_str());

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



