/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Basis.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BASIS_HPP_
#define SRC_HMR_CL_HMR_BASIS_HPP_

#include "cl_HMR_Edge.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Facet.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "HMR_Globals.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Cell.hpp" //CNT/src
#include "cl_Matrix.hpp"

#include "cl_MTK_Vertex.hpp" //MTK/src

namespace moris::hmr
{
    //Cell< mtk::Vertex* > gEmptyVertexCell;
    //------------------------------------------------------------------------------
    /**
     * \brief base class for templated Lagrange Nodes and B-Splines
     */
    class Basis : public mtk::Vertex
    {
            //------------------------------------------------------------------------------

        protected:

            //------------- -----------------------------------------------------------------

            //! Level on which basis is defined
            const uint       mLevel;

            //! owner of basis
            moris_id         mOwner = gNoProcOwner;

            //! processors whom share basis (in ascending proc rank order)
            Matrix< IdMat > mSharingProcs;

            //! counts how many elements are connected to this basis
            uint             mNumberOfConnectedElements = 0;

            //! multi purpose flag
            bool             mFlag = false;

            //! ID of basis on global domain
            luint            mDomainID = gNoEntityID;

            //! global index in whole domain
            luint            mDomainIndex = gNoEntityID;

            //! index in local memory
            luint            mMemoryIndex = gNoEntityID;

            //! index on local proc for MTK
            luint            mLocalIndex = gNoEntityID;

            //! flag telling if node is used by owned elements
            bool             mUsedFlag = false;

            //! flag telling if node is used by owned and shared elements
            bool             mUsedOwnedAndSharedFlag = false;

            //  array containing connected elements
            moris::Cell< Element * > mElements;

            //! counts how many facets are connected to this basis
            uint             mNumberOfConnectedFacets = 0;

            //  array containing connected facets
            Facet**          mFacets = nullptr;

            //! counts how many edges are connected to this basis
            uint             mNumberOfConnectedEdges = 0;

            //  array containing connected facets
            Edge**           mEdges = nullptr;

            // -----------------------------------------------------------------------------

        public:

            // -----------------------------------------------------------------------------

            /**
             * default basis constructor
             *
             * @param[in]   aLevel        level on which basis exists
             * @param[in]   aOwner        owner of basis
             */

            Basis( uint aLevel,
                    uint aOwner ) : mLevel( aLevel ),
                    mOwner( aOwner ),
                    mSharingProcs(0,0,MORIS_INDEX_MAX)
        {
        }

            // -----------------------------------------------------------------------------

            /**
             * Virtual destructor. Does nothing.
             */
            virtual ~Basis(){};

            // -----------------------------------------------------------------------------

            /**
             * MTK Interface: returns owner of basis
             *
             * @return uint    ID of proc that owns this basis
             */
            moris_id get_owner() const
            {
                return mOwner;
            }
            // -----------------------------------------------------------------------------

            /**
             * MTK Interface: returns all procs which share this basis
             *
             * @return const Matrix<IdMat> &    IDs of proc that share this basis
             */
            const Matrix<IdMat> & get_node_sharing() const
            {
                return mSharingProcs;
            }
            // -----------------------------------------------------------------------------

            /**
             * Add node sharing for processor
             *
             * @return void
             */
            void add_node_sharing( moris_id aSharedProcRank )
            {
                uint tNumShared = mSharingProcs.n_rows();
                mSharingProcs.resize( tNumShared + 1, 1 );
                mSharingProcs( tNumShared, 0 ) = aSharedProcRank;
            }
            // -----------------------------------------------------------------------------
            /**
             * Return whether node is shared
             *
             * @return bool, true = has node sharing, false = does not have node sharing
             */
            bool has_node_sharing()
            {
                if( mSharingProcs.numel()>0 )
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            // -----------------------------------------------------------------------------
            /**
             * MTK Interface: returns a domain wide id of the vertex
             */
            moris_id get_id() const
            {
                // fixme: add +1 and check against MTK output
                return mDomainIndex + 1 ; // < -- this is correct
                // HMR's domain index is MTK's domain id +1

                //return mDomainID;
            }

            // -----------------------------------------------------------------------------

            /**
             * MTK Interface: returns a local proc index of the vertex
             */
            virtual moris_index get_index() const
            {
                return mLocalIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the flag of the basis to true
             *
             * @return void
             */
            void flag()
            {
                mFlag = true;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the active flag of the basis to false
             *
             * @return void
             */
            void unflag()
            {
                mFlag = false;
            }

            //------------------------------------------------------------------------------

            /**
             * tests if the basis is flagged
             *
             * @return bool
             */
            bool is_flagged() const
            {
                return mFlag;
            }

            // -----------------------------------------------------------------------------

            /**
             * returns the level of a basis
             *
             * @return   uint level of basis
             */
            uint get_level() const
            {
                return mLevel;
            }

            // -----------------------------------------------------------------------------

            /**
             * Sets the proc owner of basis to specified value
             *
             * @param[in]  aOwner   ID of proc that will own this basis
             *
             * @return void
             */
            void set_owner( uint aOwner )
            {
                mOwner = aOwner;
            }

            // -----------------------------------------------------------------------------

            /**
             * increment the element counter
             *
             * @return void
             */
            void increment_element_counter()
            {
                ++mNumberOfConnectedElements;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the element counter to zero
             *
             * @return void
             */
            void reset_element_counter()
            {
                mNumberOfConnectedElements = 0;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the value of the element counter
             */
            auto get_element_counter() const
            -> decltype ( mNumberOfConnectedElements )
            {
                return mNumberOfConnectedElements;
            }

            // -----------------------------------------------------------------------------

            /**
             * increment the element counter
             *
             * @return void
             */
            void increment_facet_counter()
            {
                ++mNumberOfConnectedFacets;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the element counter to zero
             *
             * @return void
             */
            void reset_facet_counter()
            {
                mNumberOfConnectedFacets = 0;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the value of the element counter
             */
            auto get_facet_counter() const
            -> decltype ( mNumberOfConnectedFacets )
            {
                return mNumberOfConnectedFacets;
            }

            //------------------------------------------------------------------------------

            /**
             * increment the element counter
             *
             * @return void
             */
            void increment_edge_counter()
            {
                ++mNumberOfConnectedEdges;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the element counter to zero
             *
             * @return void
             */
            void reset_edge_counter()
            {
                mNumberOfConnectedEdges = 0;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the value of the element counter
             */
            auto get_edge_counter() const
            -> decltype ( mNumberOfConnectedEdges )
            {
                return mNumberOfConnectedEdges;
            }

            //------------------------------------------------------------------------------

            /**
             * set ID of global domain
             *
             * param[in]  aIndex    new index of basis
             *
             * @return void
             */
            void set_domain_id( luint aID )
            {
                mDomainID = aID;
            }

            //------------------------------------------------------------------------------

            /**
             * get index of global domain
             *
             * @return luint global index of basis
             */
            auto get_hmr_id() const
            -> decltype( mDomainID  )
            {
                return mDomainID;
            }

            //------------------------------------------------------------------------------

            /**
             * set index of  global domain
             *
             * param[in]  aIndex    new index of basis
             *
             * @return void
             */
            void set_domain_index( luint aIndex )
            {
                mDomainIndex = aIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the local index that is needed for MTK
             */
            void set_local_index( luint aIndex )
            {
                mLocalIndex = aIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * get index of global domain
             *
             * @return luint global index of basis
             */
            auto get_hmr_index()
            -> decltype( mDomainIndex  )
            {
                return mDomainIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the used flag of this basis to true
             *
             * @return void
             */
            void use()
            {
                mUsedFlag = true;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the used flag of this basis to false
             *
             * @return void
             */
            void unuse()
            {
                mUsedFlag = false;
            }

            //------------------------------------------------------------------------------

            /**
             * tells if this basis is used by curreny proc
             *
             * @return bool
             */
            auto is_used() const
            -> decltype( mUsedFlag )
            {
                return mUsedFlag;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the owned and shared used flag of this basis to true
             *
             * @return void
             */
            void use_owned_and_shared()
            {
                mUsedOwnedAndSharedFlag = true;
            }

            //------------------------------------------------------------------------------

            /**
             * sets the owned and shared used flag of this basis to false
             *
             * @return void
             */
            void unuse_owned_and_shared()
            {
                mUsedOwnedAndSharedFlag = false;
            }

            //------------------------------------------------------------------------------

            /**
             * tells if this basis, owned and shared,is used by curreny proc
             *
             * @return bool
             */
            auto is_use_owned_and_shared() const
            -> decltype( mUsedOwnedAndSharedFlag )
            {
                return mUsedOwnedAndSharedFlag;
            }

            //------------------------------------------------------------------------------

            /**
             * Returns an array of size [N] telling the proc local ijk-position
             * of the basis on the current level.
             *
             * @return luint pointer to array containing ijk-position
             *               careful: element must not go out of scope.
             */
            virtual const luint * get_ijk( ) const = 0;

            // ----------------------------------------------------------------------------

            /**
             * set XYZ coordinates
             *
             * @param[in] aXYZ    array containing coordinates
             *
             * @return void
             */
            virtual void set_xyz( const real * aXYZ ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * get XYZ coordinates
             *
             * @return real*
             */
            virtual const real * get_xyz() const = 0;

            //------------------------------------------------------------------------------

            /**
             * sets the value of the memory index
             *
             * @param[in] aMemoryIndex
             *
             * @return void
             */
            void set_memory_index( luint aMemoryIndex )
            {
                mMemoryIndex = aMemoryIndex;
            }
            //------------------------------------------------------------------------------

            /**
             * returns the value of the memory index
             */
            auto get_memory_index() const
            -> decltype ( mMemoryIndex )
            {
                return mMemoryIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the element container
             * and resets the memory counter
             *
             * @return void
             */
            void init_element_container()
            {
                // assign memory to container
                if ( mNumberOfConnectedElements != 0 )
                {
                    mElements.resize( mNumberOfConnectedElements, nullptr );
                }

                // reset counter
                mNumberOfConnectedElements = 0;
            }

            //------------------------------------------------------------------------------

            /**
             * tell the node that it is connected to the element
             *
             * copies the pointer into mElements and increments counter
             *
             * @return void
             */
            void insert_element( Element* aElement )
            {
                mElements( mNumberOfConnectedElements++ ) = aElement;
            }

            //------------------------------------------------------------------------------

            uint get_num_elements()
            {
                return mElements.size();
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            Element * get_element( uint aIndex )
            {
                return mElements( aIndex );
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element ( const version )
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            const Element * get_element( uint aIndex ) const
            {
                return mElements( aIndex );
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the element container
             * and resets the memory counter
             *
             * @return void
             */
            void init_facet_container()
            {
                if ( mNumberOfConnectedFacets != 0 )
                {
                    // assign memory to container
                    mFacets = new Facet* [ mNumberOfConnectedFacets ];

                    // reset counter
                    mNumberOfConnectedFacets = 0;
                }
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the element container
             * and resets the memory counter
             *
             * @return void
             */
            void delete_facet_container()
            {
                if( mNumberOfConnectedFacets != 0 )
                {
                    delete [] mFacets;

                    // reset counter
                    mNumberOfConnectedFacets = 0;
                }
            }

            //------------------------------------------------------------------------------

            /**
             * tell the node that it is connected to the element
             *
             * copies the pointer into mElements and increments counter
             *
             * @return void
             */
            void insert_facet( Facet* aFacet )
            {
                mFacets[ mNumberOfConnectedFacets++ ] = aFacet;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            Facet * get_facet( uint aIndex )
            {
                return mFacets[ aIndex ];
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element ( const version )
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            const Facet * get_facet( uint aIndex ) const
            {
                return mFacets[ aIndex ];
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the element container
             * and resets the memory counter
             *
             * @return void
             */
            void init_edge_container()
            {
                if ( mNumberOfConnectedEdges != 0 )
                {
                    // assign memory to container
                    mEdges = new Edge* [ mNumberOfConnectedEdges ];

                    // reset counter
                    mNumberOfConnectedEdges = 0;
                }
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the element container
             * and resets the memory counter
             *
             * @return void
             */
            void delete_edge_container()
            {
                if( mNumberOfConnectedEdges != 0 )
                {
                    delete [] mEdges;

                    // reset counter
                    mNumberOfConnectedEdges = 0;
                }
            }

            //------------------------------------------------------------------------------

            /**
             * tell the node that it is connected to the element
             *
             * copies the pointer into mElements and increments counter
             *
             * @return void
             */
            void insert_edge( Edge* aEdge )
            {
                mEdges[ mNumberOfConnectedEdges++ ] = aEdge;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            Edge * get_edge( uint aIndex )
            {
                return mEdges[ aIndex ];
            }

            //------------------------------------------------------------------------------

            /**
             * returns a pointer to the linked element ( const version )
             *
             * @param[in]  aIndex   number of element that is requested
             *
             * @return     Element_Base*    pointer to connected element
             */
            const Edge * get_edge( uint aIndex ) const
            {
                return mEdges[ aIndex ];
            }

            //------------------------------------------------------------------------------

            virtual void set_active_flag()
            {
                MORIS_ERROR( false, "set_active_flag() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void set_refined_flag()
            {
                MORIS_ERROR( false, "set_refined_flag() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void set_deactive_flag()
            {
                MORIS_ERROR( false, "set_deactive_flag() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual bool is_active()
            {
                MORIS_ERROR( false, "is_active() not available for selected basis type." );
                return false;
            }

            //------------------------------------------------------------------------------

            virtual bool is_refined()
            {
                MORIS_ERROR( false, "is_refinded() not available for selected basis type." );
                return false;
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the neighbor container
             *
             * @return void
             */
            virtual void init_neighbor_container( )
            {
                MORIS_ERROR( false, "init_neighbor_container() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the neighbor container
             *
             * @return void
             */
            virtual void delete_neighbor_container( )
            {
                MORIS_ERROR( false, "delete_neighbor_container() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            /**
             * reserves the memory for the neighbor container
             *
             * @return void
             */
            virtual void init_children_container( )
            {
                MORIS_ERROR( false, "init_children_container() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void increment_parent_counter()
            {
                MORIS_ERROR( false, "increment_parent_counter() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void insert_parent( Basis * aParent )
            {
                MORIS_ERROR( false, "insert_parent() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual Basis * get_parent( uint aParentNumber )
            {
                MORIS_ERROR( false, "get_parent() not available for selected basis type." );
                return nullptr;
            }
            //------------------------------------------------------------------------------

            virtual uint get_number_of_parents()
            {
                MORIS_ERROR( false, "get_number_of_parents() not available for selected basis type." );
                return 0;
            }

            //------------------------------------------------------------------------------

            virtual void insert_neighbor( uint aaNeighborNumber,
                    Basis * aNeighbor )
            {
                MORIS_ERROR( false, "insert_neighbor() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual Basis * get_neighbor( uint aNeighborNumber )
            {
                MORIS_ERROR( false, "get_neighbor() not available for selected basis type." );
                return nullptr;
            }

            //------------------------------------------------------------------------------

            virtual void insert_child(
                    uint aChildNumbner,
                    Basis      * aChild )
            {
                MORIS_ERROR( false, "insert_child() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual Basis * get_child( uint aChildNumber )
            {
                MORIS_ERROR( false, "get_child() not available for selected basis type." );
                return nullptr;
            }

            //------------------------------------------------------------------------------

            virtual bool has_children()
            {
                MORIS_ERROR( false, "has_children() not available for selected basis type." );
                return false;
            }

            //------------------------------------------------------------------------------

            virtual void flag_descendants( )
            {
                MORIS_ERROR( false, "flag_descendants() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void unflag_descendants( )
            {
                MORIS_ERROR( false, "unflag_descendants() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void collect_descendants(
                    Cell< Basis* > & aBasisList,
                    luint          & aBasisCount )
            {
                MORIS_ERROR( false, "collect_descendants() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void count_descendants( luint & aBasisCount )
            {
                MORIS_ERROR( false, " count_descendants() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual void set_active_index( luint aIndex )
            {
                MORIS_ERROR( false, "set_active_index() not available for selected basis type." );
            }

            //------------------------------------------------------------------------------

            virtual luint get_active_index()
            {
                MORIS_ERROR( false, "get_active_index() not available for selected basis type." );
                return gNoEntityID;
            }

            //------------------------------------------------------------------------------

            virtual mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
            {
                MORIS_ERROR( false, "get_interpolation() not available for for selected basis type.");
                return nullptr;
            }

            //------------------------------------------------------------------------------

            virtual bool has_interpolation( const uint aBSplineMeshIndex )
            {
                MORIS_ERROR( false, "has_interpolation() not available for for selected basis type.");
                return false;
            }

            //------------------------------------------------------------------------------

            virtual const mtk::Vertex_Interpolation * get_interpolation(  const uint aBSplineMeshIndex ) const
            {
                MORIS_ERROR( false, "get_interpolation() const not available for for selected basis type.");
                return nullptr;
            }

            //------------------------------------------------------------------------------

            /**
             * set the DOFs
             */
            virtual void set_coefficients( const uint                   aBSplineMeshIndex,
                    Cell< mtk::Vertex* > & aDOFs )
            {
                MORIS_ERROR( false, "set_coefficients() not available for for selected basis type.");
            }

            // ----------------------------------------------------------------------------

            /**
             * set the T-Matrix coefficients
             */
            virtual void set_weights( const uint               aBSplineMeshIndex,
                    const Matrix< DDRMat > & aTMatrix )
            {
                MORIS_ERROR( false, "set_weights() not available for for selected basis type.");
            }

            //------------------------------------------------------------------------------
            // virtual Mat< moris_id >
            // get_adof_ids() const
            // {
            //    MORIS_ERROR( false, "get_adof_ids() not available for for selected basis type.");
            //    return Mat< moris_id >(0,0);
            // }

            //------------------------------------------------------------------------------

            // virtual Mat< moris_index >
            // get_adof_indices() const
            //{
            //    MORIS_ERROR( false, "get_adof_indices() not available for for selected basis type.");
            //    return Mat< moris_index >(0,0);
            //}

            //------------------------------------------------------------------------------

            //virtual Matrix< DDUMat >
            //get_adof_owners() const
            // {
            //     MORIS_ERROR( false, "get_adof_owners() not available for for selected basis type.");
            //     return Matrix< DDUMat >(0,0);
            // }

            //------------------------------------------------------------------------------

            //virtual Cell< mtk::Vertex* > &
            //get_adof_pointers()
            // {
            //     MORIS_ERROR( false, "get_adof_pointers() not available for for selected basis type.");
            //     return gEmptyVertexCell;
            // }

            //------------------------------------------------------------------------------

            //virtual const Cell< mtk::Vertex* > &
            //get_adof_pointers() const
            //{
            //    MORIS_ERROR( false, "get_adof_pointers() const not available for for selected basis type.");
            //     return gEmptyVertexCell;
            // }

            //------------------------------------------------------------------------------

            /*virtual const Matrix< DDRMat > *
         get_weights( const uint aOrder ) const
         {
             MORIS_ERROR( false, "get_weights() not available for for selected basis type.");
             return nullptr;
         }*/

            //------------------------------------------------------------------------------

            virtual Matrix< DDRMat > get_coords() const
            {
                MORIS_ERROR( false, "get_coords() not available for for selected basis type.");
                return Matrix< DDRMat >(0,0);
            }

            //------------------------------------------------------------------------------

            virtual void get_basis_local_child_inds( Matrix< DDSMat > & aChildren )
            {
                MORIS_ERROR( false, "get_basis_local_child_inds() not available for for selected basis type.");
            }

            //------------------------------------------------------------------------------

            virtual void init_interpolation( uint aOrder )
            {
                MORIS_ERROR( false, "init_interpolation() not available for for selected basis type.");
            }

    };
    //------------------------------------------------------------------------------

} /* namespace moris */

#endif

