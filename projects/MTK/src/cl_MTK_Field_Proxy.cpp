
#include "cl_MTK_Field_Proxy.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Proxy::Field_Proxy(
                std::shared_ptr<mtk::Mesh_Manager>   aMeshManager,
                uint const                         & aMeshIndex,
                uint const                         & aDiscretizationMeshIndex )
        {
            mMeshManager = aMeshManager;
            mMeshIndex = aMeshIndex;
            mDiscretizationMeshIndex = aDiscretizationMeshIndex;
        }

        // ----------------------------------------------------------------------------------------------

        Field_Proxy::~Field_Proxy()
        {

        }

        // ----------------------------------------------------------------------------------------------

        uint Field_Proxy::get_discretization_order() const
        {
            return mMeshManager->
                    get_interpolation_mesh( mMeshIndex )->
                    get_HMR_lagrange_mesh()->
                    get_bspline_order( mDiscretizationMeshIndex );
        }

        // ----------------------------------------------------------------------------------------------

        Matrix< DDRMat > & Field_Proxy::get_node_values()
        {
            return mNodalValues;
        }

        // ----------------------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Proxy::get_node_values() const
        {
            return mNodalValues;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > & Field_Proxy::get_coefficients()
        {
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Proxy::get_coefficients() const
        {
            return mCoefficients;
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Proxy::evaluate_node_values()
        {
            Interpolation_Mesh* tInterpolationMesh =
                    mMeshManager->get_interpolation_mesh( mMeshIndex );

            // get number of nodes on block
            uint tNumberOfNodes= tInterpolationMesh->get_num_nodes();

            // set size of node values
            mNodalValues.set_size( tNumberOfNodes, 1 );

            for( uint Ik = 0; Ik < tNumberOfNodes; ++Ik )
            {
                // get pointer to node
                auto tNode = &tInterpolationMesh->get_mtk_vertex( Ik );

                // get PDOFs from node
                auto tBSplines = tNode->
                        get_interpolation( mDiscretizationMeshIndex )->
                        get_coefficients();

                // get T-Matrix
                const Matrix< DDRMat > & tTMatrix = *tNode->
                        get_interpolation( mDiscretizationMeshIndex )->
                        get_weights();

                // get number of coefficients
                uint tNumberOfCoeffs = tTMatrix.length();

                MORIS_ASSERT( tNumberOfCoeffs > 0, "No coefficients defined for node" ) ;

                // fill coeffs vector
                Matrix< DDRMat > tCoeffs( tNumberOfCoeffs, 1 );
                for( uint Ii = 0; Ii < tNumberOfCoeffs; ++Ii )
                {
                    tCoeffs( Ii ) = mCoefficients( tBSplines( Ii )->get_index() );
                }

                // write value into solution
                mNodalValues( Ik ) = moris::dot( tTMatrix, tCoeffs );
            }
        }

        // ----------------------------------------------------------------------------------------------

    }
}

