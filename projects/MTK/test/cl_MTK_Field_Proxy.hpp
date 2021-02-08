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
        class Field_Proxy : public Field
        {
        private:
            Mesh_Pair mMeshPair;
            Matrix< DDRMat > mNodalValues;

        public :

            Field_Proxy(
                    Mesh_Pair aMeshPair,
                    uint      aDiscretizationMeshIndex = 0 );

            ~Field_Proxy();

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the mesh pair associated with this field.
             *
             * @return Mesh pair
             */
            Mesh_Pair get_mesh_pair();

            // ----------------------------------------------------------------------------------------------

            //            void save_field_to_hdf5( const std::string & aFilePath, const bool aCreateNewFile=true );
            //
            //            void save_node_values_to_hdf5( const std::string & aFilePath, const bool aCreateNewFile=true );
            //
            //            void load_field_from_hdf5( const std::string & aFilePath,
            //                    const uint          aBSplineOrder=0 );
            //
            //            void save_bspline_coeffs_to_binary( const std::string & aFilePath );
            //
            //            void save_node_values_to_binary( const std::string & aFilePath );

            template<typename T>
            void evaluate_scalar_function( T aLambda )
            {
                Interpolation_Mesh* tInterpolationMesh = mMeshPair.mInterpolationMesh;

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
