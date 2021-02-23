/*
 * cl_MTK_Field_Proxy.hpp
 *
 *  Created on: Jan 20, 2021
 *      Author: schmidt
 */
#ifndef SRC_MTK_FIELD_PROXY_HPP_
#define SRC_MTK_FIELD_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        class Field_Analytic : public Field
        {
            private:

            public :

                // ----------------------------------------------------------------------------------------------

                Field_Analytic(
                        mtk::Mesh_Pair * aMeshPairs,
                        uint const     & aDiscretizationMeshIndex = 0 );

                // ----------------------------------------------------------------------------------------------

                ~Field_Analytic();

                // ----------------------------------------------------------------------------------------------

                template<typename T>
                void evaluate_scalar_function( T aLambda )
                {
                    Interpolation_Mesh* tInterpolationMesh =
                            mMeshPair->mInterpolationMesh;

                    // get number of nodes on block
                    uint tNumberOfVertices = tInterpolationMesh->get_num_nodes();

                    // set size of node values
                    mNodalValues.set_size( tNumberOfVertices, 1 );

                    // loop over all vertices
                    for( uint Ik = 0; Ik < tNumberOfVertices; ++Ik )
                    {
                        // evaluate function at vertex coordinates
                        mNodalValues( Ik ) = aLambda( tInterpolationMesh->get_mtk_vertex( Ik ).get_coords() );
                    }
                }
        };
    }
}
#endif /* SRC_MTK_FIELD_PROXY_HPP_ */
