/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field_Proxy.hpp
 *
 */

#ifndef SRC_MTK_FIELD_PROXY_HPP_
#define SRC_MTK_FIELD_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::mtk
{
    class Field_Proxy : public Field
    {
      private:

      public:
        // ----------------------------------------------------------------------------------------------

        Field_Proxy(
                const mtk::Mesh_Pair& aMeshPairs,
                uint const &          aDiscretizationMeshIndex = 0 );

        // ----------------------------------------------------------------------------------------------

        ~Field_Proxy() override;

        // ----------------------------------------------------------------------------------------------

        template< typename T >
        void evaluate_scalar_function( T aLambda )
        {
            Interpolation_Mesh* tInterpolationMesh =
                    mMeshPair.get_interpolation_mesh();

            // get number of nodes on block
            uint tNumberOfVertices = tInterpolationMesh->get_num_nodes();

            // set size of node values
            mValues.set_size( tNumberOfVertices, 1 );

            // loop over all vertices
            for ( uint Ik = 0; Ik < tNumberOfVertices; ++Ik )
            {
                // evaluate function at vertex coordinates
                mValues( Ik ) = aLambda( tInterpolationMesh->get_mtk_vertex( Ik ).get_coords() );
            }
        }
    };
    }

#endif /* SRC_MTK_FIELD_PROXY_HPP_ */

