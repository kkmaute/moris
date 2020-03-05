/*
 * cl_XTK_Multigrid.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: Schmidt
 */

#include "cl_XTK_Multigrid.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>

#include "cl_Communication_Tools.hpp"
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_sort.hpp"
#include "typedefs.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "xtk_typedefs.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Model.hpp"

namespace xtk
{

    Multigrid::Multigrid( xtk::Model * aXTKModelPtr ) : mXTKModelPtr( aXTKModelPtr )
    {
//        mXTKModelPtr
    }

//------------------------------------------------------------------------------

    void Multigrid::create_fine_to_coarse_relationship()
    {
        moris::mtk::Interpolation_Mesh & tInterpolationMesh = mXTKModelPtr->get_background_mesh().get_mesh_data();

        uint tNumBasis = tInterpolationMesh.get_num_basis( 0 );

        for ( uint Ik = 0; Ik < tNumBasis; Ik ++ )
        {
            // get num coarse basis for this basis
            uint tNumCoarseBasis = tInterpolationMesh.get_num_coarse_basis_of_basis( 0, Ik );

            std::cout<<tNumCoarseBasis<<" tNumCoarseBasis"<<std::endl;

            for ( uint Ii = 0; Ii < tNumCoarseBasis; Ii ++ )
            {
                moris_index tCoarseBasisIndex = tInterpolationMesh.get_coarse_basis_index_of_basis( 0, Ik, Ii );

                std::cout<<tCoarseBasisIndex<<" tCoarseBasisIndex"<<std::endl;
            }

        }



        std::cout<<tNumBasis<<std::endl;






//        const Cell<Matrix<IndexMat>> & tEnrichedCoeffsToBackrountCoeffs = mXTKModelPtr->mEnrichedIntegMesh( 0 )
//                                               ->get_enriched_coefficients_to_background_coefficients();                  //FIXME

    }

}


