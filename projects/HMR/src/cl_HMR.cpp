/*
 * cl_HMR.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */
#include "op_times.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Interface.hpp" //HMR/src
#include "cl_HMR_MTK.hpp" //HMR/src

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

            // delete B-Spline and Lagrange meshes
            this->delete_meshes();

            // delete Background Mesh
            if( mBackgroundMesh != NULL )
            {
                delete mBackgroundMesh;
            }

            for( auto tField: mFields )
            {
                delete tField;
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
        HMR::synchronize_t_matrix_flags()
        {
            // synchronize flags for T-Matrices on aura
            mBackgroundMesh->synchronize_t_matrix_flags();

            // get proc neighbors
            auto tMyProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            uint tNumberOfNeighbors = tMyProcNeighbors.length();

            uint tMyRank = par_rank();

            for( auto tMesh: mBSplineMeshes )
            {
                // ask mesh about number of basis per element
                auto tNumberOfBasisPerElement
                    = tMesh->get_number_of_basis_per_element();

                // loop over all procs
                for ( uint p=0; p<tNumberOfNeighbors; ++p )
                {
                    // only do this if there is a neighbor
                    if(        tMyProcNeighbors( p ) != gNoProcNeighbor
                            && tMyProcNeighbors( p ) != tMyRank )
                    {
                        Cell< Background_Element_Base* > tElements;

                        // get active elements from aura
                        mBackgroundMesh->collect_active_elements_from_aura(
                                p, 0, tElements );

                        // loop over all elements from aura
                        for( auto tElement : tElements )
                        {
                            // test if element is flagged
                            if( tElement->get_t_matrix_flag() )
                            {
                                // get pointer to B-Spline Element
                                auto tBElement = tMesh->get_element_by_memory_index(
                                        tElement->get_memory_index() );

                                // loop over all basis of B-Spline element
                                for( uint k=0; k<tNumberOfBasisPerElement; ++k )
                                {
                                    // get pointer to basis
                                    auto tBasis = tBElement->get_basis( k );

                                    // test if basis is owned by current proc
                                    if( tBasis->get_owner() == tMyRank )
                                    {
                                        // find first element that is owned
                                        auto tNumberOfElements = tBasis->get_element_counter();

                                        for( uint i=0; i<tNumberOfElements; ++i )
                                        {
                                            auto tOtherElement = tBasis->get_element( i );

                                            // test if this element is owned
                                            if ( tOtherElement->get_owner() == tMyRank )
                                            {
                                                // flag T-Matrix of this element
                                                tOtherElement->get_background_element()->set_t_matrix_flag();

                                                // exit loop
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        HMR::finalize()
        {
            // mesh counter
            uint tMeshCount = 0;

            // synchronize flahs for T-Matrices
            this->synchronize_t_matrix_flags();

            // loop over all Lagrange meshes
            for( auto tMesh : mLagrangeMeshes )
            {
                // unflag all nodes on this mesh
                tMesh->unflag_all_basis();

                // unflag all nodes on B-Spline mesh
                mBSplineMeshes( tMeshCount )->unflag_all_basis();

                // get number of elements on background mesh
                auto tNumberOfElements = mBackgroundMesh
                      ->get_number_of_active_elements_on_proc();

                // calculate transposed Lagrange T-Matrix
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

                        mTMatrix( tMeshCount )->calculate_t_matrix(
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

                        // loop over all nodes of this element
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

                // increment mesh counter
                ++tMeshCount;
            }

            // create the communication table for this mesh
            this->create_communication_table();
        }

// -----------------------------------------------------------------------------

        void
        HMR::create_communication_table()
        {
            uint tParSize = par_size();
            uint tMyRank  = par_rank();

            if( tParSize > 1 )
            {
                // in a first step, we identify all processers this proc wants
                // to talk to

                // this is a Bool-like matrix
                Mat< uint > tColumn( tParSize, 1, 0 );

                // test owners of B-Splines
                for( auto tMesh: mBSplineMeshes )
                {
                    // get number of active B-Splines
                    auto tNumberOfBSplines = tMesh->get_number_of_active_basis_on_proc();

                    // loop over all active basis on this mesh
                    for( uint k=0; k<tNumberOfBSplines; ++k )
                    {
                        // get pointer to basis
                        auto tBasis = tMesh->get_active_basis( k );

                        // test if flag of basis is set
                        if ( tBasis->is_flagged() )
                        {
                            // set flag for this proc
                            tColumn( tBasis->get_owner() ) = 1;
                        }
                    }
                }

                // remove self from row
                tColumn( tMyRank ) = 0;

                // communication table
                Mat< uint > tCommTable;

                // matrices to send
                Cell< Mat< uint > > tSend;

                // matrices to receive
                Cell< Mat< uint > > tRecv;

                if( tMyRank != 0 )
                {
                    // create communication table with one entry
                    tCommTable.set_size( 1, 1, 0 );
                    tSend.resize( 1, tColumn );
                }
                else
                {
                    // create comm matrix
                    tCommTable.set_size( tParSize, 1, 0 );

                    // communicate with all other procs
                    for( uint k=1; k<tParSize; ++k )
                    {
                        tCommTable( k ) = k;
                    }

                    // nothing to send
                    Mat< uint > tEmpty;
                    tSend.resize( tParSize, tEmpty );
                }

                // exchange matrices
                communicate_mats( tCommTable, tSend, tRecv );

                // process information on master proc
                if ( tMyRank == 0 )
                {
                    // create communication matrix
                    Mat< uint > tCommMatrix( tParSize, tParSize, 0 );

                    // process first row
                    tRecv( 0 ) = tColumn;

                    // loop over all procs and create comm matrix
                    for( uint j=0; j<tParSize; ++j )
                    {
                        for( uint i=0; i<tParSize; ++i )
                        {
                            if ( tRecv( j )( i, 0 ) != 0 )
                            {
                                tCommMatrix( i, j ) = 1;
                                tCommMatrix( j, i ) = 1;
                            }
                        }
                    }

                    // remove diagonal
                    for( uint i=0; i<tParSize; ++i )
                    {
                        tCommMatrix( i, i ) = 0;
                    }

                    // create sending list
                    Mat< uint > tEmpty;
                    tSend.resize( tParSize, tEmpty );

                    for( uint j=0; j<tParSize; ++j )
                    {
                        // count nonzero entries
                        uint tCount = 0;
                        for( uint i=0; i<tParSize; ++i )
                        {
                            if ( tCommMatrix( i, j ) != 0 )
                            {
                                ++tCount;
                            }
                        }

                        // assign memory
                        tSend( j ).set_size( tCount, 1, 0 );

                        // reset counter
                        tCount = 0;

                        // write values into matrix
                        for( uint i=0; i<tParSize; ++i )
                        {
                            if ( tCommMatrix( i, j ) != 0 )
                            {
                                tSend( j )( tCount++ ) = i;
                            }
                        }
                    }


                }

                // exchange matrices
                communicate_mats( tCommTable, tSend, tRecv );

                if ( tMyRank == 0 )
                {
                    mCommunicationTable = tSend( 0 );
                }
                else
                {
                    mCommunicationTable = tRecv( 0 );
                }
            }
            else
            {
                // output is empty
                mCommunicationTable.set_size( 0, 1 );
            }

        }
// -----------------------------------------------------------------------------

        void
        HMR::add_field( const std::string & aLabel,
                     const uint & aOrder,
                     const Mat<real> & aValues )
        {
            // create new field
            Field * tField = new Field(
                    mParameters,
                    aLabel,
                    mBackgroundMesh,
                    mBSplineMeshes( aOrder-1 ),
                    mLagrangeMeshes( aOrder -1 ) );
            // copy values
            tField->set_lagrange_values( aValues );

            // add to database
            mFields.push_back( tField );
        }

// -----------------------------------------------------------------------------

        /**
         * this function is for testing purpose only. Data is always copied.
         * This is not an efficient way to do things!
         */
        void
        HMR::save_to_exodus( const uint & aOrder, const std::string & aPath )
        {
            // create MTK object
            MTK * tMTK = mLagrangeMeshes( aOrder -1 )->create_mtk_object();

            // append fiends
            for( auto tField : mFields )
            {
                if( tField->get_order() == aOrder )
                {
                    tMTK->add_node_data(
                            tField->get_label(),
                            tField->get_data() );
                }
            }

            // save MTK to exodus
            tMTK->save_to_file( aPath );

            // delete file
            delete tMTK;
        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
