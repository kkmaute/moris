/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */
#include "op_times.hpp" //LNA/src
#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Interface.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        HMR::HMR ( const Parameters * aParameters ) :
                mParameters( aParameters )
        {

            // create factory object
            Factory tFactory;

            // create background mesh object
            mBackgroundMesh = tFactory.create_background_mesh( mParameters );

            // update element table
            mBackgroundMesh->collect_active_elements();

            // initialize mesh objects
            this->create_meshes();

            // initialize T-Matrix objects
            this->init_t_matrices();
        }

// -----------------------------------------------------------------------------

        HMR::~HMR()
        {
            // delete T-Matrix objects
            this->delete_t_matrices();

            // delete B-Spline and lagrange meshes
            this->delete_meshes();

            // delete Background Mesh
            if( mBackgroundMesh )
            {
                delete mBackgroundMesh;
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::delete_meshes()
        {

            // delete all pointers
            for( auto tMesh : mBSplineMeshes )
            {
                // test if mesh exists
                if( tMesh != NULL )
                {
                    // delete this mesh
                    delete tMesh;
                }
            }

            // delete all pointers
            for( auto tMesh : mLagrangeMeshes )
            {
                // test if mesh exists
                if( tMesh != NULL )
                {
                    // delete this mesh
                    delete tMesh;
                }
            }

        }
// -----------------------------------------------------------------------------

        void
        HMR::create_meshes()
        {
            // delete existing meshes
            this->delete_meshes();

            // create factory object
            Factory tFactory;

            // get max interpolation degree
            uint tMaxOrder = mParameters->get_max_polynomial();

            // assign memory for B-Spline meshes
            mBSplineMeshes.resize ( tMaxOrder, nullptr );

            // assign memory for Lagrange meshes
            mLagrangeMeshes.resize ( tMaxOrder, nullptr );

            // loop over all meshes
            for( uint k=0; k<tMaxOrder; ++k )
            {
                mBSplineMeshes( k )
                        = tFactory.create_bspline_mesh(
                                mParameters, mBackgroundMesh, k+1 );

                mLagrangeMeshes( k )
                        = tFactory.create_lagrange_mesh(
                                mParameters, mBackgroundMesh, k+1 );
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::update_meshes()
        {
            // update all B-Spline meshes
            for( auto tMesh : mBSplineMeshes )
            {
                // test if mesh exists
                if( tMesh != NULL )
                {
                    // synchronize mesh with background mesh
                    tMesh->update_mesh();
                }
            }

            // get counter for meshes
            uint tMeshCount = 0;

            // update all Lagrange meshes and link elements to their
            // B-Spline twins
            for( auto tMesh : mLagrangeMeshes )
            {
                // test if mesh exists
                if( tMesh != NULL )
                {
                    // synchronize mesh with background mesh
                    tMesh->update_mesh();

                    // link twins
                    if ( mBSplineMeshes( tMeshCount ) != NULL )
                    {
                        tMesh->link_twins( mBSplineMeshes( tMeshCount ) );
                    }
                }

                // increment counter
                ++tMeshCount;
            }
        }

// -----------------------------------------------------------------------------

        Interface
        HMR::create_interface()
        {
            return Interface( *this );
        }

// -----------------------------------------------------------------------------

        void
        HMR::init_t_matrices()
        {
            // get number of meshes
            uint tNumberOfMeshes = mParameters->get_max_polynomial();

            // allocate T-Matrix cell
            mTMatrix.resize( tNumberOfMeshes, nullptr );

            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                // test if both meshes exist
                if ( mBSplineMeshes( k ) != NULL
                    && mLagrangeMeshes( k ) != NULL )
                {
                    // initialize T-Matrix object
                    mTMatrix( k ) = new T_Matrix( mParameters,
                            mBSplineMeshes( k ),
                            mLagrangeMeshes( k ) );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::delete_t_matrices()
        {
            // delete pointers of calculation objects
            for( auto tTMatrix :  mTMatrix )
            {
                if ( tTMatrix != NULL )
                {
                    delete tTMatrix;
                }
            }
        }
// -----------------------------------------------------------------------------

        void
        HMR::calculate_t_matrices()
        {
            // mesh counter
            uint tMeshCount = 0;
            // loop over all Lagrange meshes
            for( auto tMesh : mLagrangeMeshes )
            {
                // unflag all nodes on this mesh
                tMesh->unflag_all_basis();

                // unflag all nodes on B-Spline mesh
                mBSplineMeshes( tMeshCount )->unflag_all_basis();

                // test if mesh exists
                if( tMesh != NULL )
                {
                    // get number of elements on background mesh
                    auto tNumberOfElements = mBackgroundMesh
                            ->get_number_of_active_elements_on_proc();

                    // calculate trabsposed Lagrange T-Matrix
                    Mat< real > tL( mTMatrix( tMeshCount )->get_lagrange_matrix() );

                    auto tNumberOfNodes = tMesh->get_number_of_basis_per_element();

                    // loop over all active elements
                    for( uint e=0; e<tNumberOfElements; ++e )
                    {
                        // get pointer to element on background mesh
                        auto tElement = mBackgroundMesh->get_element( e );

                        // test if element is flagged
                        if ( tElement->get_t_matrix_flag() )
                        {
                            // get index of element
                            auto tMemoryIndex = tElement->get_memory_index() ;

                            // calculate the B-Spline T-Matrix
                            Mat< real > tB;
                            Cell< Basis* > tDOFs;

                            mTMatrix( tMeshCount )->calculate_truncated_t_matrix(
                                    tMemoryIndex,
                                    tB,
                                    tDOFs );

                            // transposed T-Matrix
                            Mat< real > tT = tL * tB;

                            // number of columns in T-Matrix
                            uint tNCols = tT.n_cols();

                            // get pointer to Lagrange Element
                            auto tLagrangeElement
                                = tMesh->get_element_by_memory_index( tMemoryIndex );

                            // epsilon to count T-Matrix
                            real tEpsilon = 1e-12;

                            // loop over all nodes of this basis
                            for( uint k = 0; k<tNumberOfNodes; ++k  )
                            {
                                // pointer to node
                                auto tNode = tLagrangeElement->get_basis( k );

                                // test if node is flagged
                                if ( ! tNode->is_flagged() )
                                {
                                    // initialize counter
                                    uint tCount = 0;

                                    // count number of nonzero entries
                                    for( uint i=0; i<tNCols; ++i )
                                    {
                                        if ( std::abs( tT( k, i ) ) > tEpsilon )
                                        {
                                            // increment counter
                                            ++tCount;
                                        }
                                    }

                                    // reserve DOF cell
                                    Cell< mtk::Vertex* > tNodeDOFs( tCount, nullptr );

                                    // reserve matrix with coefficients
                                    Mat< real > tCoefficients( tCount, 1 );

                                    // reset counter
                                    tCount = 0;

                                    // loop over all nonzero entries
                                    for( uint i=0; i<tNCols; ++i )
                                    {
                                        if ( std::abs( tT( k, i ) ) > tEpsilon )
                                        {
                                            // copy entry of T-Matrix
                                            tCoefficients( tCount ) = tT( k, i );

                                            // copy pointer of dof and convert to mtk::Vertex
                                            tNodeDOFs( tCount ) = tDOFs( i );

                                            // flag this DOF
                                            tDOFs( i )->flag();

                                            // increment counter
                                            ++tCount;
                                        }
                                    }

                                    // store the coefficients
                                    tNode->set_t_matrix( tCoefficients );

                                    // store pointers to the DOFs
                                    tNode->set_dofs( tNodeDOFs );

                                    // flag this node as processed
                                    tNode->flag();
                                }
                            }
                        }
                    }

                    // calculate indices for flagged basis
                    mBSplineMeshes( tMeshCount )->calculate_basis_indices();
                }

                // increment mesh counter
                ++tMeshCount;
            }
        }

// -----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
