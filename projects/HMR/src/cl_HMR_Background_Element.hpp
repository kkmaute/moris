/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Element.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BACKGROUND_ELEMENT_HPP_
#define SRC_HMR_CL_HMR_BACKGROUND_ELEMENT_HPP_

#include "cl_HMR_Background_Edge.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Background_Facet.hpp"
#include "typedefs.hpp" //COR/src
#include "cl_Cell.hpp" //CNT/src
#include "cl_Bitset.hpp" //CNT/src
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp" //LINALG/src

namespace moris::hmr
{

    /**
     * Background element templated against number of dimensions
     *
     * @tparam N Number of dimensions
     */
    template< uint N >
    class Background_Element : public Background_Element_Base
    {
    private:

        //! Number of children
        static constexpr uint C = ( N == 1 ? 2 : ( N == 2 ? 4 : 8 ) );

        //! Number of neighbors
        static constexpr uint B = ( N == 1 ? 2 : ( N == 2 ? 8 : 26 ) );

        //! Number of faces
        static constexpr uint F = 2 * N;

        //! Number of edges
        static constexpr uint E = ( N == 3 ? 12 : 0 );

        //! Fixed size array containing children.
        //! If the element has no children, it contains null pointers
        Background_Element_Base* mChildren[ C ] = { nullptr };

        //! bitset defining the child index // FIXME why is this an bitset and not an uint.
        Bitset< N >              mChildBitset;

        //! Fixed size array containing neighbors on same level
        //! Active or refined.
        Background_Element_Base* mNeighbors[ B ] = { nullptr };

        //! local ijk position on proc
        luint                    mIJK[ N ];

        //! faces for this element
        Background_Facet *       mFacets[ F ] = { nullptr };

        //! owner bitset
        Bitset< F >              mFacetOwnFlags;

        //! edges for this element
        Background_Edge **       mEdges;

        //! edges bitset
        Bitset< E >              mEdgeOwnFlags;

//--------------------------------------------------------------------------------
    public:
//--------------------------------------------------------------------------------

         /**
          * Default element constructor
          *
          * @param[in]  aParent       Pointer to parent of element,
          *                           null pointer if element is on level zero.
          *
          * @param[in]  aIJK          Pointer to IJK position of element.
          * @param[in]  aID           Global ID of element.
          * @param[in]  aSubDomainID  Proc local subdomain ID of element.
          * @param[in]  aLevel        Level of element.
          * @param[in]  aChildIndex   Tells which child number this element is.
          *                           Pass UINT_MAX if element is on level zero
          * @param[in]  aOwner        For zero level elements, pass UINT_MAX,
          *                           owner is determined later. For all higher
          *                           levels, pass owner of parent.
          */
        Background_Element(       Background_Element_Base * aParent,
                            uint aActivePattern,
                            const luint                   * aIJK,
                            luint aID,
                            uint aLevel,
                            uint aChildIndex,
                            uint aOwner ) : Background_Element_Base( aParent,
                                                                                                aActivePattern,
                                                                                                aID,
                                                                                                aLevel,
                                                                                                aOwner )
        {
            // write child index into bitset
            this->set_child_index( aChildIndex );

            // save ijk position in memory.
            // Needed to calculate ijk of children.
            for( uint k=0; k<N; ++k )
            {
                mIJK[ k ] = aIJK[ k ];
            }

            // reset owner flags
            mFacetOwnFlags.reset();

            // Creates an size 12 Background_edge pointer array in 3D. Does nothing in 2D
            this->init_edge_container();

            // this only applies for elements of level 1 and above
            if( aLevel > 0 )
            {
                // check if element needs to be refined again
                if( aParent->get_min_refimenent_level() > aLevel )
                {
                    // set min refinement level
                    this->set_min_refimenent_level( aParent->get_min_refimenent_level() );

                    // flag this element for refinement
                    this->put_on_refinement_queue();
                }
            }
        }

//--------------------------------------------------------------------------------

        /**
         * Default element destructor.
         * Deletes children if children flag is set.
         */
        ~Background_Element()
        {
            if ( mChildrenFlag )
            {
                for ( auto p:  mChildren )
                {
                    delete p;
                }
            }
            // delete faces
           this->delete_facets();

           // delete edges
           this->delete_edge_container();
        }
//--------------------------------------------------------------------------------

        /**
         * Returns an array of size [N] telling the proc local ijk-position
         * of the element on the current level.
         *
         * @return luint pointer to array containing ijk-position
         *               careful: element must not go out of scope.
         */
        const luint *
        get_ijk( ) const
        {
            return mIJK;
        }

//--------------------------------------------------------------------------------

        /**
         * This function is called by the background mesh during refinement.
         * The child pointer is put into the member array mChildren.
         *
         * @param[in] aChild   pointer of child to be added
         *
         * @return void
         */
        void insert_child(  Background_Element_Base* aChild )
        {
            mChildrenFlag = true;
            mChildren[ aChild->get_child_index() ] = aChild;
        }

//--------------------------------------------------------------------------------

        /**
         * This function is needed by the background mesh during refinement.
         * The ijk-position of a child is needed to calculate the local and
         * global ID of the child element.
         *
         * param[ out ]  Matrix< DDUMat > of size <number of dimensions>
         *                                *< number of children >
         * @return void
         */
        void get_ijk_of_children( Matrix< DDLUMat > & aIJK ) const;

//--------------------------------------------------------------------------------

        /**
         * Recursive function that counts all active descendants of an element.
         * If the element is not refined, the function returns one, not zero.
         * This function is needed by the background mesh in order
         * to determine the size of mActiveElements.
         *
         * @param[inout] aCount   Counter to be incremented
         *
         * @return void
         */
        void
        get_number_of_active_descendants( uint aPattern, luint & aCount ) const;

//--------------------------------------------------------------------------------
        /**
         * Recursive function that counts all descendants of an element plus
         * the element itself.
         *
         * param[inout] aCount   Counter to be incremented
         *
         * @return void
         */
        void
        get_number_of_descendants( luint & aCount ) const;

//--------------------------------------------------------------------------------

        /**
         * To be called after the cell aElementList has been allocated
         * to the size given by  get_number_of_descendants().
         * Returns an array that consists all related elements, including
         * the element itself.
         *
         * @param[inout] aElementList  cell to which the pointers are added
         * @param[inout] aCount   Counter to be incremented
         *
         * @return void
         */
        void
        collect_descendants( Cell< Background_Element_Base* > & aElementList,
                             luint                            & aElementCount );

//--------------------------------------------------------------------------------

        /**
         * To be called after the cell aElementList has been allocated
         * to the size given by  get_number_of_active_descendants().
         * Needed by the background mesh to update mActiveElements.
         *
         * @param[inout] aElementList  cell to which the pointers are added
         * @param[inout] aCount   Counter to be incremented
         *
         * @return void
         *
         */
        void collect_active_descendants( uint aPattern,
                                               Cell< Background_Element_Base *> & aElementList,
                                               luint                            & aElementCount );

        void collect_active_descendants( uint aPattern,
                                               Cell< const Background_Element_Base *> & aElementList,
                                               luint                                  & aElementCount ) const ;
//--------------------------------------------------------------------------------

        /**
         * To be called after the cell aElementList has been allocated
         * to the size given by  get_number_of_active_descendants().
         * Needed by the background mesh to update mActiveElements.
         *
         * @param[in]    aPattern      activation scheme this operation is performed on
         * @param[inout] aElementList  Matrix with memory indices of elements
         * @param[inout] aCount        Counter to be incremented
         *
         * @return void
         *
         */
        void collect_active_descendants_by_memory_index(
                uint aPattern,
                Matrix< DDLUMat >                & aElementList,
                luint                            & aElementCount,
                int                         aNeighborIndex=-1) const;

//--------------------------------------------------------------------------------

        /**
         * Provided my background cell side ordinal, return the neighbor's background cell
         * side ordinal
         *
         * @return int neighbor side ordinal
         *
         */
        int get_neighbor_side_ordinal( int aNeighborIndex) const;

//--------------------------------------------------------------------------------

        /**
         * Provided my background cell side ordinal, return the child cell ordinals
         * on side
         *
         * @return int neighbor child cell ordinal
         *
         */
        void get_child_cell_ordinals_on_side( int                aSideOrdinal,
                                              Matrix<IndexMat> & aChildCellOrdinals) const;

//--------------------------------------------------------------------------------

        uint get_num_children() const
        {
            return C;
        }

//--------------------------------------------------------------------------------

        /**
         * Returns a pointer to a child of an element. If the element
         * has no children, a null pointer will be returned.
         *
         * @param[in] aIndex      Index of requested child (2D: 0-3, 3D: 0-7)
         * @return Background_Element_Base*  pointer to selected child
         */
        Background_Element_Base * get_child( uint aIndex )
        {
            return mChildren[ aIndex ];
        }

//--------------------------------------------------------------------------------

        /**
         * Returns a pointer to a child of an element. If the element
         * has no children, a null pointer will be returned.
         *
         * @param[in] aIndex      Index of requested child (2D: 0-3, 3D: 0-7)
         * @return Background_Element_Base*  pointer to selected child
         */
        Background_Element_Base * get_child( uint aIndex ) const
        {
            return mChildren[ aIndex ];
        }

//--------------------------------------------------------------------------------

        /**
         * This function is called by the background mesh.
         * The pointer of the neighbor is put into the member array mNeighbors.
         *
         * @param[in] aIndex      number of neighbor
         * @param[in] aNeighbor   pointer of neighbor to be added
         *
         * @return void
         */
        void insert_neighbor( uint aIndex,
                              Background_Element_Base * aNeighbor )
        {
            mNeighbors[ aIndex ] = aNeighbor;
        }

//--------------------------------------------------------------------------------

        /**
         * Returns a pointer to a neighbor of an element.
         *
         * @param[ in ] aNeighborIndex       index of requested neighbor
         * @return Background_Element_Base*  pointer to requested neighbor
         */
        Background_Element_Base * get_neighbor( uint aIndex )
        {
            return mNeighbors[ aIndex ];
        }

//--------------------------------------------------------------------------------

        /**
         * Returns a pointer to a neighbor of an element ( const version )
         *
         * @param[ in ] aNeighborIndex       index of requested neighbor
         * @return Background_Element_Base*  pointer to requested neighbor
         */
        const Background_Element_Base* get_neighbor( uint aIndex ) const
        {
            return mNeighbors[ aIndex ];
        }

//--------------------------------------------------------------------------------

        /**
         * Recursive function that loops up to a specified level and counts
         * active and refined elements on that level.
         *
         * @param[in]     aLevel    level to be considered
         * @param[inout]  aElementCount    counter for elements
         */
        void count_elements_on_level( uint aLevel,
                                      luint & aElementCount );
//--------------------------------------------------------------------------------

        /**
         * Recursive function that loops up to a specified level and collects
         * active and refined elements on that level.
         *
         * @param[in]     aLevel          level to be considered
         * @param[inout]  aElementList    cell to which the pointers are added
         * @param[inout]  aElementCount   counter for elements
         */
        void collect_elements_on_level( uint                               aLevel,
                                        Cell< Background_Element_Base* > & aElementList,
                                        luint                            & aElementCount );
//--------------------------------------------------------------------------------

        void collect_neighbors( uint aPattern );

//--------------------------------------------------------------------------------

        uint get_number_of_facets() const
        {
            return F;
        }

//--------------------------------------------------------------------------------

        uint get_number_of_edges() const
        {
            return E;
        }

//--------------------------------------------------------------------------------

        /**
         * Function for debugging that prints all neighbors of the element
         * to the screen.
         *
         * @return void
         *
         */
        void print_neighbors( uint aPattern )
        {
            // print header
            std::fprintf( stdout, "\n Neighbors of Element %lu ( ID: %lu, mem:  %lu, child: %u ): \n\n",
                          ( long unsigned int ) mDomainIndex(aPattern),
                          ( long unsigned int ) mDomainID,
                          ( long unsigned int ) mMemoryIndex,
                          ( unsigned int )      this->get_child_index() );
            for( uint k=0; k<B; ++k )
            {
                // test if neighbor exists
                if( mNeighbors[ k ] )
                {
                    // get id of neighbor
                    long unsigned int tID = mNeighbors[ k ]->get_hmr_id();

                    // get active flag
                    int tActive = mNeighbors[ k ]->is_active( aPattern );

                    // print index and id of neighbor
                    std::fprintf( stdout, "    %2u   id: %5lu     a: %d\n",
                            (unsigned int) k, (long unsigned int) tID, tActive);
                }
                else
                {
                    std::fprintf( stdout, "    %2u  null \n",  (unsigned int) k );
                }
            }

            // print blank line
            std::fprintf( stdout, "\n" );
        }

//--------------------------------------------------------------------------------

        /**
         * Tells which child number an element has.
         * This information is for example needed for the T-Matrix.
         *
         * @return uint    index between 0 and 3 (2D) or 0 and 7 (3D)
         */
        uint get_child_index() const
        {
            return mChildBitset.to_ulong();
        }

//--------------------------------------------------------------------------------

        /**
         * Test a bit in the child index bitset
         *
         * @param[in]    aBit    index of bit to test
         *
         * @return bool
         */
        bool test_child_index_bit( uint aBit ) const
        {
            return mChildBitset.test( aBit );
        }

//--------------------------------------------------------------------------------

        /**
         * creates a bitset that describes the pedigree path
         *
         */
        void encode_pedigree_path( luint            & aAncestorID,
                                    Matrix< DDUMat > & aPedigreeList,
                                    luint            & aCounter )
        {
            // create bitset
            Bitset< gBitsetSize > tBitset;

            // clear bitset
            tBitset.reset();

            // get pointer to element
            Background_Element_Base* tParent = this;

            // get level of this element
            uint tLevel = this->get_level();

            // copy level into list
            aPedigreeList( aCounter++ ) = tLevel;

            if ( tLevel > 0 )
            {
                uint tMax = tLevel * N -1;

                // Position in bitset
                uint tPivot = tMax;

                // reset output bitset
                tBitset.reset();

                // loop over all levels and calculate index
                for( uint k=tLevel; k>0; --k )
                {
                    for( int i=N-1; i>=0; --i ) // do not use uint here
                    {
                        if ( tParent->test_child_index_bit( i ) )
                        {
                            // set this bit
                            tBitset.set( tPivot );
                        }

                        // decrement pivot
                        --tPivot;
                    }
                    // get next element
                    tParent = tParent->get_parent();
                }

                // reset pivot
                tPivot = 0;

                while ( tPivot <= tMax )
                {
                    Bitset< 32 > tChar;
                    for( uint k=0; k<32; ++k )
                    {
                        if( tBitset.test( tPivot++ ) )
                        {
                            tChar.set( k );
                        }
                        if ( tPivot > tMax )
                        {
                            break;
                        }
                    }

                    // copy character into output list
                    aPedigreeList( aCounter++ ) = ( uint ) tChar.to_ulong();
                }
            }

            // copy domain ID from ancestor
            aAncestorID = tParent->get_hmr_id();
        }

//--------------------------------------------------------------------------------

        luint get_length_of_pedigree_path()
        {
            return ceil( ( ( real ) ( this->get_level()*N ) ) / 32 ) + 1;
        }

//--------------------------------------------------------------------------------

        /**
         * Returns a cell with pointers to elements on the same level,
         * if they exist.
         *
         * @param[ in  ] aOrder       degree of neighborship
         * @param[ out ] aNeighbors   cell containing found neighbors
         */
        void get_neighbors_from_same_level( uint aOrder,
                                                  Cell< Background_Element_Base * > & aNeighbors );

//--------------------------------------------------------------------------------
        void delete_facets()
        {
            // loop over all faces
            for( uint f=0; f<F; ++f )
            {
                if( mFacetOwnFlags.test( f ) )
                {
                    if(mFacets[ f ] != nullptr)
                    {
                        delete mFacets[ f ];
                        mFacets[ f ] = nullptr;
                        mFacetOwnFlags.reset( f );
                    }
                }
            }
        }
//--------------------------------------------------------------------------------

        void delete_edges()
        {
            // loop over all edges
            for( uint e=0; e<E; ++e )
            {
                if( mEdgeOwnFlags.test( e ) )
                {
                    if(mEdges[ e ] != nullptr)
                    {
                        delete mEdges[ e ];
                        mEdges[ e ] = nullptr;
                        mEdgeOwnFlags.reset( e );
                    }
                }
            }
        }

//--------------------------------------------------------------------------------
        void reset_neigbors()
        {
            // loop over all faces
            for( uint f=0; f<B; ++f )
            {
                if(mNeighbors[ f ] != nullptr)
                {
                    mNeighbors[ f ] = nullptr;
                }
            }
        }
//--------------------------------------------------------------------------------
        /**
         * create the faces of this element
         */
        void create_facets()
        {
            //                 this->delete_facets();

            // Index on neighbor element
            uint tIndexOnOther[ 6 ] = { 2, 3, 0, 1, 5, 4 } ;

            // loop over all faces
            for( uint f = 0; f < F; ++f )
            {
                uint tOther =  tIndexOnOther[ f ];

                // grab pointer of facet from neighbor, if neighbor exists and is on same level.
                if( mNeighbors[ f ] != nullptr )
                {
                    // test if neighbor lives on same level
                    if(  mNeighbors[ f ]->get_level() == mLevel )
                    {
                        // copy pointer of facet. May be null
                        mFacets[ f ] = mNeighbors[ f ]->get_facet( tOther );
                    }
                }
                // test if facet has not been created yet
                if ( mFacets[ f ] == nullptr )
                {
                    // set flag that this element is responsible for
                    // deleting this face
                    mFacetOwnFlags.set( f );

                    // test if this element is on the edge of the domain
                    if( mNeighbors[ f ] == nullptr )
                    {
                        // this facet has now proc owner, I am leader
                        // create face
                        mFacets[ f ] = new Background_Facet( this, f );
                    }
                    else if( mNeighbors[ f ]->get_level() != mLevel )
                    {
                        // this element belongs to the creator
                        mFacets[ f ] = new Background_Facet( this, f );
                    }
                    else
                    {
                        // element picks owner with lower domain_index
                        mFacets[ f ] = new Background_Facet( this, mNeighbors[ f ], f  );

                        // insert face into neighbor
                        mNeighbors[ f ]->insert_facet(  mFacets[ f ], tOther );
                    }
                }
                else
                {
                    mFacetOwnFlags.reset( f );
                }
            } // end loop over all faces
        }

//--------------------------------------------------------------------------------

        /**
         * creates the edges ( 3D only )
         */
        void create_edges();

//--------------------------------------------------------------------------------

        /**
         * reset the flags of the faces
         */
        void reset_flags_of_facets();

//--------------------------------------------------------------------------------

        /**
         * returns a face of the background element
         */
        Background_Facet * get_facet( uint aIndex )
        {
            return mFacets[ aIndex ];
        }

//--------------------------------------------------------------------------------

        /**
         * inserts a face into the background element
         */
        void insert_facet( Background_Facet * aFace, uint aIndex )
        {
            MORIS_ASSERT( mFacets[ aIndex ] == nullptr, "tried to overwrite existing facet" );
            // copy face to slot
            mFacets[ aIndex ] = aFace;
        }

//--------------------------------------------------------------------------------

        /**
         * returns an edge of the background element ( 3D only )
         */
        Background_Edge * get_edge( uint aIndex );

//--------------------------------------------------------------------------------

        void insert_edge(  Background_Edge * aEdge, uint aIndex );

//-------------------------------------------------------------------------------

        void reset_flags_of_edges();

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_1(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_2(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_3(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_4(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_5(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void get_number_of_active_descendants_on_side_6(
                uint aPattern,
                      luint & aCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_1(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_2(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_3(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_4(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_5(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//--------------------------------------------------------------------------------

        void collect_active_descendants_on_side_6(
                uint aPattern,
                Cell< Background_Element_Base* > & aElementList,
                luint                            & aElementCount );

//-------------------------------------------------------------------------------

        void init_edge_container();

//--------------------------------------------------------------------------------
    private:
//--------------------------------------------------------------------------------

        void set_child_index( uint aIndex );

//-------------------------------------------------------------------------------

        void delete_edge_container();
//--------------------------------------------------------------------------------
    };
//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::init_edge_container()
    {
        // do nothing
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::init_edge_container()
    {
        // assign memory for edge container and fill with null pointers
        mEdges = new Background_Edge * [ 12 ]{};

        // reset owner flags
        mEdgeOwnFlags.reset();
    }
//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::delete_edge_container()
    {
        // do nothing
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::delete_edge_container()
    {
        // delete edges
        this->delete_edges();

        // delete edge container
        delete [] mEdges;
    }

//--------------------------------------------------------------------------------
    template< uint N >
    inline
    void Background_Element< N >::get_neighbors_from_same_level(
            uint aOrder,
            Cell< Background_Element_Base * > & aNeighbors )
    {
        MORIS_ERROR( false, "Don't know how search neighbors on same level.");
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::set_child_index( uint aIndex )
    {
        MORIS_ERROR( false, "Don't know how to set child index.");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 2 >::set_child_index( uint aIndex )
    {
        switch( aIndex )
        {
            case( 0 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.reset( 1 );
                break;
            }
            case( 1 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.reset( 1 );
                break;
            }

            case( 2 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.set( 1 );
                break;
            }

            case( 3 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.set( 1 );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid child index.");
                break;
            }
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::set_child_index( uint aIndex )
    {
        switch( aIndex )
        {
            case( 0 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.reset( 1 );
                mChildBitset.reset( 2 );
                break;
            }
            case( 1 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.reset( 1 );
                mChildBitset.reset( 2 );
                break;
            }
            case( 2 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.set( 1 );
                mChildBitset.reset( 2 );
                break;
            }
            case( 3 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.set( 1 );
                mChildBitset.reset( 2 );
                break;
            }
            case( 4 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.reset( 1 );
                mChildBitset.set( 2 );
                break;
            }
            case( 5 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.reset( 1 );
                mChildBitset.set( 2 );
                break;
            }
            case( 6 ) :
            {
                mChildBitset.reset( 0 );
                mChildBitset.set( 1 );
                mChildBitset.set( 2 );
                break;
            }
            case( 7 ) :
            {
                mChildBitset.set( 0 );
                mChildBitset.set( 1 );
                mChildBitset.set( 2 );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid child index.");
                break;
            }
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::get_ijk_of_children( Matrix< DDLUMat > & aIJK ) const
    {
        MORIS_ERROR( false, "Don't know how to calculate ijk of children.");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template <>
    inline
    void Background_Element< 1 >::get_ijk_of_children(
            Matrix< DDLUMat > & aIJK ) const
    {
        // set size of IJK output
        aIJK.set_size( 1, 2 );

        // child 0
        aIJK( 0, 0 ) = 2*mIJK[0];

        // child 1
        aIJK( 0, 1 ) = aIJK( 0, 0 ) + 1;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 2 >::get_ijk_of_children( Matrix< DDLUMat > & aIJK ) const
    {
        // set size of IJK output
        aIJK.set_size( 2, 4 );

        // child 0
        aIJK( 0, 0 ) = 2*mIJK[0];
        aIJK( 1, 0 ) = 2*mIJK[1];

        // child 1
        aIJK( 0, 1 ) = aIJK( 0, 0 ) + 1;
        aIJK( 1, 1 ) = aIJK( 1, 0 );

        // child 2
        aIJK( 0, 2 ) = aIJK( 0, 0 );
        aIJK( 1, 2 ) = aIJK( 1, 0 ) + 1;

        // child 3
        aIJK( 0, 3 ) = aIJK( 0, 1 );
        aIJK( 1, 3 ) = aIJK( 1, 2 );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::get_ijk_of_children( Matrix< DDLUMat > & aIJK ) const
    {
        // set size of IJK output
        aIJK.set_size( 3, 8 );

        // child 0
        aIJK( 0, 0 ) = 2*mIJK[0];
        aIJK( 1, 0 ) = 2*mIJK[1];
        aIJK( 2, 0 ) = 2*mIJK[2];

        // child 1
        aIJK( 0, 1 ) = aIJK( 0, 0 ) + 1;
        aIJK( 1, 1 ) = aIJK( 1, 0 );
        aIJK( 2, 1 ) = aIJK( 2, 0 );

        // child 2
        aIJK( 0, 2 ) = aIJK( 0, 0 );
        aIJK( 1, 2 ) = aIJK( 1, 0 ) + 1;
        aIJK( 2, 2 ) = aIJK( 2, 0 );

        // child 3
        aIJK( 0, 3 ) = aIJK( 0, 1 );
        aIJK( 1, 3 ) = aIJK( 1, 2 );
        aIJK( 2, 3 ) = aIJK( 2, 0 );

        // child 4
        aIJK( 0, 4 ) = aIJK( 0, 0 );
        aIJK( 1, 4 ) = aIJK( 1, 0 );
        aIJK( 2, 4 ) = aIJK( 2, 0 ) + 1;

        // child 5
        aIJK( 0, 5 ) = aIJK( 0, 0 ) + 1;
        aIJK( 1, 5 ) = aIJK( 1, 0 );
        aIJK( 2, 5 ) = aIJK( 2, 4 );

        // child 6
        aIJK( 0, 6 ) = aIJK( 0, 0 );
        aIJK( 1, 6 ) = aIJK( 1, 0 ) + 1;
        aIJK( 2, 6 ) = aIJK( 2, 4 );

        // child 7
        aIJK( 0, 7 ) = aIJK( 0, 1 );
        aIJK( 1, 7 ) = aIJK( 1, 2 );
        aIJK( 2, 7 ) = aIJK( 2, 4 );
    }
//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants(
            uint aPattern,
                  luint & aCount ) const
    {
        MORIS_ERROR( false, "Don't know how to count active descendants.");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 1 >::get_number_of_active_descendants(
            uint aPattern,
                  luint & aCount ) const
    {
        // test if element is active
        if ( mActiveFlags.test( aPattern ) )
        {
            // increment counter
            ++aCount ;
        }
        else if( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 1 ]->get_number_of_active_descendants( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 2 >::get_number_of_active_descendants(
            uint aPattern,
                  luint & aCount ) const
    {
        // test if element is active
        if ( mActiveFlags.test( aPattern ) )
        {
            // increment counter
            ++aCount ;
        }
        else if( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 1 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 2 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants(
            uint aPattern,
                  luint & aCount ) const
    {
        // test if element is active
        if ( mActiveFlags.test( aPattern ) )
        {
            // increment counter
            ++aCount ;
        }
        else if( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 1 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 2 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 4 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 5 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 6 ]->get_number_of_active_descendants( aPattern, aCount );
            mChildren[ 7 ]->get_number_of_active_descendants( aPattern, aCount );
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::get_number_of_descendants(
            luint & aCount ) const
    {
        MORIS_ERROR( false, "Don't know how to count descendants.");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 1 >::get_number_of_descendants(
            luint & aCount ) const
    {
        // add self to counter
        ++aCount;

        // test if children exist
        if ( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_descendants( aCount );
            mChildren[ 1 ]->get_number_of_descendants( aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 2 >::get_number_of_descendants(
            luint & aCount ) const
    {
        // add self to counter
        ++aCount;

        // test if children exist
        if ( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_descendants( aCount );
            mChildren[ 1 ]->get_number_of_descendants( aCount );
            mChildren[ 2 ]->get_number_of_descendants( aCount );
            mChildren[ 3 ]->get_number_of_descendants( aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::get_number_of_descendants(
            luint & aCount ) const
    {
        // add self to counter
        ++aCount ;

        // test if children exist
        if ( mChildrenFlag )
        {
            // ask children to increment counter
            mChildren[ 0 ]->get_number_of_descendants( aCount );
            mChildren[ 1 ]->get_number_of_descendants( aCount );
            mChildren[ 2 ]->get_number_of_descendants( aCount );
            mChildren[ 3 ]->get_number_of_descendants( aCount );
            mChildren[ 4 ]->get_number_of_descendants( aCount );
            mChildren[ 5 ]->get_number_of_descendants( aCount );
            mChildren[ 6 ]->get_number_of_descendants( aCount );
            mChildren[ 7 ]->get_number_of_descendants( aCount );
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::collect_descendants(
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "Don't know how to collect descendants.");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 1 >::collect_descendants(
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        // add self to list
        aElementList( aElementCount++ ) = this;

        // test if children exist
        if ( mChildrenFlag )
        {
            // add children to list
            mChildren[ 0 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 1 ]->collect_descendants( aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 2 >::collect_descendants(
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        // add self to list
        aElementList( aElementCount++ ) = this;

        // test if children exist
        if ( mChildrenFlag )
        {
            // add children to list
            mChildren[ 0 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 1 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 2 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 3 ]->collect_descendants( aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::collect_descendants(
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        // add self to list
        aElementList( aElementCount++ ) = this;

        // test if children exist
        if ( mChildrenFlag )
        {
            // add children to list
            mChildren[ 0 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 1 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 2 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 3 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 4 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 5 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 6 ]->collect_descendants( aElementList, aElementCount );
            mChildren[ 7 ]->collect_descendants( aElementList, aElementCount );
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::collect_active_descendants(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "Don't know how to collect active descendants.");
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::collect_active_descendants(
            uint aPattern,
            Cell< const Background_Element_Base* > & aElementList,
            luint                            & aElementCount ) const
    {
        MORIS_ERROR( false, "Don't know how to collect active descendants ( const ).");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 1 >::collect_active_descendants(
           uint aPattern,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if ( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 1 >::collect_active_descendants(
           uint aPattern,
           Cell< const Background_Element_Base* > & aElementList,
           luint                            & aElementCount ) const
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if ( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::collect_active_descendants(
           uint aPattern,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if ( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 2 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 3 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::collect_active_descendants(
           uint aPattern,
           Cell< const Background_Element_Base* > & aElementList,
           luint                                  & aElementCount ) const
  {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if ( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 2 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 3 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 3 >::collect_active_descendants(
           uint aPattern,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 2 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 3 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 4 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 5 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 6 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 7 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 3 >::collect_active_descendants(
           uint aPattern,
           Cell< const Background_Element_Base* > & aElementList,
           luint                            & aElementCount ) const
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = this;
       }
       else if( mChildrenFlag )
       {
           // add active children to list
           mChildren[ 0 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 1 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 2 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 3 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 4 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 5 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 6 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
           mChildren[ 7 ]->collect_active_descendants( aPattern, aElementList, aElementCount );
       }
  }

//--------------------------------------------------------------------------------

   template< uint N >
   void Background_Element< N >::collect_active_descendants_by_memory_index(
           uint aPattern,
           Matrix< DDLUMat >                & aElementList,
           luint                            & aElementCount,
           const  int                         aNeighborIndex ) const
   {
       MORIS_ERROR( false, "Don't know how to collect active descendants.");
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::collect_active_descendants_by_memory_index(
           uint aPattern,
           Matrix< DDLUMat >                & aElementList,
           luint                            & aElementCount,
           const  int                        aNeighborIndex ) const
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = mMemoryIndex;
       }
       else if ( mChildrenFlag )
       {
           switch ( aNeighborIndex )
           {
               case( 0 ) :
               {
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 1 ) :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 2 ) :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 3 ) :
               {
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               default :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
           }
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 3 >::collect_active_descendants_by_memory_index(
           uint aPattern,
           Matrix< DDLUMat >                & aElementList,
           luint                            & aElementCount,
           const  int                        aNeighborIndex ) const
   {
       // test if self is active
       if ( mActiveFlags.test( aPattern ) )
       {
           // add self to list
           aElementList( aElementCount++ ) = mMemoryIndex;
       }
       else if( mChildrenFlag )
       {
           switch ( aNeighborIndex )
           {
               case( 0 ) :
               {
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 6 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 7 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 1 ) :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 4 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 6 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 2 ) :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 4 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 5 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 3 ) :
               {
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 5 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 7 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 4 ) :
               {
                   mChildren[ 4 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 5 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 6 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 7 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               case( 5 ) :
               {
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex  );
                   break;
               }
               default :
               {
                   // add active children to list
                   mChildren[ 0 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 1 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 2 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 3 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 4 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 5 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 6 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   mChildren[ 7 ]->collect_active_descendants_by_memory_index( aPattern, aElementList, aElementCount, aNeighborIndex );
                   break;
               }

           }

       }
    }

//--------------------------------------------------------------------------------

   template< uint N >
   int Background_Element< N >::get_neighbor_side_ordinal( const  int aNeighborIndex ) const
   {
       MORIS_ERROR( false, "Don't know how to get neighbor side ordinal.");
       return -1;
   }
//--------------------------------------------------------------------------------
   template <>
   inline
   int Background_Element< 2 >::get_neighbor_side_ordinal( const  int aNeighborIndex ) const
   {
           switch ( aNeighborIndex )
           {
               case( 0 ) :
               {
                   return 2;
                   break;
               }
               case( 1 ) :
               {
                   return 3;
                   break;
               }
               case( 2 ) :
               {
                   return 0;
                   break;
               }
               case( 3 ) :
               {
                   return 1;
                   break;
               }
               default :
               {
                   MORIS_ERROR(0,"Invalid side ordinal specified for this element type (0-3) for a 2D quad type cell");
                   return -1;
                   break;
               }
       }
   }

//--------------------------------------------------------------------------------

   template <>
   inline
   int Background_Element< 3 >::get_neighbor_side_ordinal( const  int aNeighborIndex ) const
   {

       switch ( aNeighborIndex )
       {
           case( 0 ) :
                   {
               return 2;
               break;
                   }
           case( 1 ) :
                   {
               return 3;
               break;
                   }
           case( 2 ) :
                   {
               return 0;
               break;
                   }
           case( 3 ) :
                   {
               return 1;
               break;
                   }
           case( 4 ) :
                   {
               return 5;
               break;
                   }
           case( 5 ) :
                   {
               return 4;
               break;
                   }
           default :
           {
               MORIS_ERROR(0,"Invalid side ordinal specified for this element type (0-5) for a 3D hex type cell ");
               return -1;
               break;
           }
       }
   }
//--------------------------------------------------------------------------------
   template< uint N >
   void Background_Element< N >::get_child_cell_ordinals_on_side( const  int        aSideOrdinal,
                                                                             Matrix<IndexMat> & aChildCellOrdinals) const
   {
       MORIS_ERROR( false, "Don't know how to get child cell ordinals on side.");
   }
//--------------------------------------------------------------------------------
   template <>
   inline
   void Background_Element< 2 >::get_child_cell_ordinals_on_side( const  int         aSideOrdinal,
                                                                              Matrix<IndexMat> & aChildCellOrdinals ) const
   {
           switch ( aSideOrdinal )
           {
               case( 0 ) :
               {
                   aChildCellOrdinals = {{0,1}};
                   break;
               }
               case( 1 ) :
               {
                   aChildCellOrdinals = {{1,3}};
                   break;
               }
               case( 2 ) :
               {
                   aChildCellOrdinals = {{3,2}};
                   break;
               }
               case( 3 ) :
               {
                   aChildCellOrdinals = {{2,0}};
                   break;
               }
               default :
               {
                   MORIS_ERROR(0,"Invalid side ordinal specified for this element type (0-3) for a 2D quad type cell");
                   break;
               }
       }
   }

//--------------------------------------------------------------------------------

   template <>
   inline
   void Background_Element< 3 >::get_child_cell_ordinals_on_side( const  int         aSideOrdinal,
                                                                                Matrix<IndexMat> & aChildCellOrdinals ) const
   {

       switch ( aSideOrdinal )
       {
           case( 0 ) :
                   {
               aChildCellOrdinals = {{0,1,5,4}};
               break;
                   }
           case( 1 ) :
                   {
               aChildCellOrdinals = {{1,3,7,5}};
               break;
                   }
           case( 2 ) :
                   {
               aChildCellOrdinals = {{3,2,6,7}};
               break;
                   }
           case( 3 ) :
                   {
               aChildCellOrdinals = {{2,0,4,6}};
               break;
                   }
           case( 4 ) :
                   {
               aChildCellOrdinals = {{2,3,1,0}};
               break;
                   }
           case( 5 ) :
                   {
               aChildCellOrdinals = {{4,5,7,6}};
               break;
                   }
           default :
           {
               MORIS_ERROR(0,"Invalid side ordinal specified for this element type (0-5) for a 3D hex type cell ");
               break;
           }
       }
   }

//--------------------------------------------------------------------------------

   template< uint N >
   void Background_Element< N >::count_elements_on_level( uint aLevel,
                                                                            luint & aElementCount )
   {
       MORIS_ERROR( false, "Don't know how to count elements on level");
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 1 >::count_elements_on_level( uint aLevel,
                                                                            luint & aElementCount )
   {
       // test if element is on specified level
       if ( mLevel ==  aLevel )
       {
           // increment counter
           ++aElementCount;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // count children
           mChildren[ 0 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 1 ]->count_elements_on_level( aLevel, aElementCount );
       }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::count_elements_on_level( uint aLevel,
                                                                            luint & aElementCount )
   {
       // test if element is on specified level
       if ( mLevel ==  aLevel )
       {
           // increment counter
           ++aElementCount;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // count children
           mChildren[ 0 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 1 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 2 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 3 ]->count_elements_on_level( aLevel, aElementCount );
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 3 >::count_elements_on_level( uint aLevel,
                                                                              luint & aElementCount )
   {
       // test if element is on specified level
       if ( mLevel ==  aLevel )
       {
           // increment counter
           ++aElementCount;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // count children
           mChildren[ 0 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 1 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 2 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 3 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 4 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 5 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 6 ]->count_elements_on_level( aLevel, aElementCount );
           mChildren[ 7 ]->count_elements_on_level( aLevel, aElementCount );
       }
   }

//--------------------------------------------------------------------------------

   template< uint N >
   void Background_Element< N >::collect_elements_on_level(
           uint aLevel,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       MORIS_ERROR( false, "Don't know how to collect elements on level");
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 1 >::collect_elements_on_level(
           uint aLevel,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       // test if element is on specified level
       if ( mLevel ==  aLevel )
       {
           // add this element to list and increment counter
           aElementList( aElementCount++ ) = this;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // collect children
           mChildren[ 0 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 1 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::collect_elements_on_level(
           uint aLevel,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
  {
       // add this element to list and increment counter
       if ( mLevel ==  aLevel )
       {
           // increment counter
           aElementList( aElementCount++ ) = this;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // collect children
           mChildren[ 0 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 1 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 2 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 3 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );
       }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 3 >::collect_elements_on_level(
           uint aLevel,
           Cell< Background_Element_Base* > & aElementList,
           luint                            & aElementCount )
   {
       // add this element to list and increment counter
       if ( mLevel ==  aLevel )
       {
           // increment counter
           aElementList( aElementCount++ ) = this;
       }
       else if ( ( mLevel < aLevel) && mChildrenFlag )
       {
           // collect children
           mChildren[ 0 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 1 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 2 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 3 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 4 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 5 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 6 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );

           mChildren[ 7 ]->collect_elements_on_level(
                   aLevel,
                   aElementList,
                   aElementCount );
       }
   }

//--------------------------------------------------------------------------------

   // fixme: neighbors do not account refinement pattern number
   template< uint N >
   void Background_Element< N >::collect_neighbors( uint aPattern )
   {
       MORIS_ERROR( false, "Don't know how to collect neighbors");
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   template <>
   inline
   void Background_Element< 2 >::collect_neighbors( uint aPattern )
   {
       switch( this->get_child_index() )
       {
           case( 0 ):
           {
               // neighbor 0
               Background_Element_Base* tNeighbor = mParent->get_neighbor( 4 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 4 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       mNeighbors[ 4 ] = tNeighbor;
                   }
               }

               // neighbors 1 and 2
               tNeighbor = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 0 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 5 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       mNeighbors[ 0 ] = tNeighbor;
                       mNeighbors[ 5 ] = tNeighbor;
                   }
               }

               // neighbors 3 to 5
               tNeighbor = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 3 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 1 ] =   mParent->get_child( 1 );
                       mNeighbors[ 7 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       mNeighbors[ 3 ] = tNeighbor;
                       mNeighbors[ 1 ] = mParent->get_child( 1 );
                       mNeighbors[ 7 ] = tNeighbor;
                   }
               }
               else
               {
                   mNeighbors[ 1 ] = mParent->get_child( 1 );
               }

               // neighbor 6
               mNeighbors[ 2 ] = mParent->get_child( 2 );

               // neighbor 7
               mNeighbors[ 6 ] = mParent->get_child( 3 );

               break;
           }
           case( 1 ):
           {
               // neighbors 0 and 1
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 4 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 0 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       mNeighbors[ 4 ] = tNeighbor;
                       mNeighbors[ 0 ] = tNeighbor;
                   }
               }

               // neighbor 2
               tNeighbor = mParent->get_neighbor( 5 );
               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 5 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       mNeighbors[ 5 ] = tNeighbor;
                   }
               }

               // neighbor 3
               mNeighbors[ 3 ] = mParent->get_child( 0 );

               // neighbor 4 to 7
               tNeighbor = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 1 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 7 ] = mParent->get_child( 2 );
                       mNeighbors[ 2 ] = mParent->get_child( 3 );
                       mNeighbors[ 6 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       mNeighbors[ 1 ] = tNeighbor;
                       mNeighbors[ 7 ] = mParent->get_child( 2 );
                       mNeighbors[ 2 ] = mParent->get_child( 3 );
                       mNeighbors[ 6 ] = tNeighbor;
                   }
               }
               else
               {
                   mNeighbors[ 7 ] =   mParent->get_child( 2 );
                   mNeighbors[ 2 ] =   mParent->get_child( 3 );
               }
               break;
           }
           case( 2 ):
           {
               // neighbors 0 to 3
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 4 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 0 ] = mParent->get_child( 0 );
                       mNeighbors[ 5 ] = mParent->get_child( 1 );
                       mNeighbors[ 3 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       mNeighbors[ 4 ] = tNeighbor;
                       mNeighbors[ 0 ] = mParent->get_child( 0 );
                       mNeighbors[ 5 ] = mParent->get_child( 1 );
                       mNeighbors[ 3 ] = tNeighbor;
                   }
               }
               else
               {
                   mNeighbors[ 0 ] =  mParent->get_child( 0 );
                   mNeighbors[ 5 ] =  mParent->get_child( 1 );
               }

               // neighbor 4
               mNeighbors[ 1 ] =  mParent->get_child( 3 );

               // neighbor 5
               tNeighbor = mParent->get_neighbor( 7 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 7 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       mNeighbors[ 7 ] = tNeighbor;
                   }
               }

               // neighbors 6 and 7
               tNeighbor = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 2 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 6 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       mNeighbors[ 2 ] = tNeighbor;
                       mNeighbors[ 6 ] = tNeighbor;
                   }
               }
               break;
           }
           case( 3 ):
           {
               // neighbor 0
               mNeighbors[ 4 ] = mParent->get_child( 0 );

               // neighbor 1
               mNeighbors[ 0 ] = mParent->get_child( 1 );

               // neighbors 2 to 4
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 5 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 3 ] =   mParent->get_child( 2 );
                       mNeighbors[ 1 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       mNeighbors[ 5 ] = tNeighbor;
                       mNeighbors[ 3 ] = mParent->get_child( 2 );
                       mNeighbors[ 1 ] = tNeighbor;
                   }
               }
               else
               {
                   mNeighbors[ 3 ] =  mParent->get_child( 2 );
               }

               // neighbors 5 and 6
               tNeighbor = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 7 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 2 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       mNeighbors[ 7 ] = tNeighbor;
                       mNeighbors[ 2 ] = tNeighbor;
                   }
               }

               // neighbor 7
               tNeighbor = mParent->get_neighbor( 6 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       mNeighbors[ 6 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       mNeighbors[ 6 ] = tNeighbor;
                   }
               }
               break;
           }
       }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline
    void Background_Element< 3 >::collect_neighbors( uint aPattern )
    {

       switch( this->get_child_index() )
       {
       case( 0 ) :
           {
               // link to siblings
               mNeighbors[  1 ] = mParent->get_child(  1 );
               mNeighbors[  2 ] = mParent->get_child(  2 );
               mNeighbors[  5 ] = mParent->get_child(  4 );
               mNeighbors[ 12 ] = mParent->get_child(  3 );
               mNeighbors[ 15 ] = mParent->get_child(  5 );
               mNeighbors[ 16 ] = mParent->get_child(  6 );
               mNeighbors[ 24 ] = mParent->get_child(  7 );

               // get neighbor 0 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 0 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 11 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 14 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor;
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 3 of parent
               tNeighbor = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 3 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 13 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 17 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor;
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 4 of parent
               tNeighbor = mParent->get_neighbor( 4 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 4 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor->get_child( 4 );
                       mNeighbors[  7 ] = tNeighbor->get_child( 5 );
                       mNeighbors[  8 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor;
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 6 of parent
               tNeighbor = mParent->get_neighbor( 6 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 6 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 6 of parent
                       mNeighbors[  6 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 6 of parent
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 9 of parent
               tNeighbor = mParent->get_neighbor( 9 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 9 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 9 of parent
                       mNeighbors[  9 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 9 of parent
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 10 of parent
               tNeighbor = mParent->get_neighbor( 10 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 10 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 10 of parent
                       mNeighbors[ 10 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 10 of parent
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 18 of parent
               tNeighbor = mParent->get_neighbor( 18 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 18 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 18 of parent
                       mNeighbors[ 18 ] = tNeighbor->get_child( 7 );
                   }
                   else
                   {
                       // link to neighbor 18 of parent
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               break;
           }
           case( 1 ) :
           {
               // link to siblings
               mNeighbors[  2 ] = mParent->get_child(  3 );
               mNeighbors[  3 ] = mParent->get_child(  0 );
               mNeighbors[  5 ] = mParent->get_child(  5 );
               mNeighbors[ 13 ] = mParent->get_child(  2 );
               mNeighbors[ 16 ] = mParent->get_child(  7 );
               mNeighbors[ 17 ] = mParent->get_child(  4 );
               mNeighbors[ 25 ] = mParent->get_child(  6 );

               // get neighbor 0 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 0 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 10 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 14 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor;
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 1 of parent
               tNeighbor = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 1 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 12 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 15 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor;
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 4 of parent
               tNeighbor = mParent->get_neighbor( 4 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 4 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor->get_child( 5 );
                       mNeighbors[  8 ] = tNeighbor->get_child( 7 );
                       mNeighbors[  9 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor;
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 6 of parent
               tNeighbor = mParent->get_neighbor( 6 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 6 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 6 of parent
                       mNeighbors[  6 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 6 of parent
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 7 of parent
               tNeighbor = mParent->get_neighbor( 7 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 7 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 7 of parent
                       mNeighbors[  7 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 7 of parent
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 11 of parent
               tNeighbor = mParent->get_neighbor( 11 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 11 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 11 of parent
                       mNeighbors[ 11 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 11 of parent
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 19 of parent
               tNeighbor = mParent->get_neighbor( 19 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 19 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 19 of parent
                       mNeighbors[ 19 ] = tNeighbor->get_child( 6 );
                   }
                   else
                   {
                       // link to neighbor 19 of parent
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               break;
           }
           case( 2 ) :
           {

               // link to siblings
               mNeighbors[  0 ] = mParent->get_child(  0 );
               mNeighbors[  1 ] = mParent->get_child(  3 );
               mNeighbors[  5 ] = mParent->get_child(  6 );
               mNeighbors[ 11 ] = mParent->get_child(  1 );
               mNeighbors[ 14 ] = mParent->get_child(  4 );
               mNeighbors[ 15 ] = mParent->get_child(  7 );
               mNeighbors[ 23 ] = mParent->get_child(  5 );

               // get neighbor 2 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 2 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 12 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 16 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor;
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 3 of parent
               tNeighbor = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 3 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 10 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 17 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor;
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 4 of parent
               tNeighbor = mParent->get_neighbor( 4 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 4 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor->get_child( 6 );
                       mNeighbors[  6 ] = tNeighbor->get_child( 4 );
                       mNeighbors[  7 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor;
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 8 of parent
               tNeighbor = mParent->get_neighbor( 8 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 8 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 8 of parent
                       mNeighbors[  8 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 8 of parent
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 9 of parent
               tNeighbor = mParent->get_neighbor( 9 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 9 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 9 of parent
                       mNeighbors[  9 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 9 of parent
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 13 of parent
               tNeighbor = mParent->get_neighbor( 13 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 13 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 13 of parent
                       mNeighbors[ 13 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 13 of parent
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 21 of parent
               tNeighbor = mParent->get_neighbor( 21 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 21 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 21 of parent
                       mNeighbors[ 21 ] = tNeighbor->get_child( 5 );
                   }
                   else
                   {
                       // link to neighbor 21 of parent
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }
               break;
           }
           case( 3 ) :
           {
               // link to siblings
               mNeighbors[  0 ] = mParent->get_child(  1 );
               mNeighbors[  3 ] = mParent->get_child(  2 );
               mNeighbors[  5 ] = mParent->get_child(  7 );
               mNeighbors[ 10 ] = mParent->get_child(  0 );
               mNeighbors[ 14 ] = mParent->get_child(  5 );
               mNeighbors[ 17 ] = mParent->get_child(  6 );
               mNeighbors[ 22 ] = mParent->get_child(  4 );

               // get neighbor 1 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 1 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 11 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 15 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor;
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 2 of parent
               tNeighbor = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 2 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 13 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 16 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor;
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 4 of parent
               tNeighbor = mParent->get_neighbor( 4 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 4 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor->get_child( 7 );
                       mNeighbors[  6 ] = tNeighbor->get_child( 5 );
                       mNeighbors[  9 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 4 of parent
                       mNeighbors[  4 ] = tNeighbor;
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 7 of parent
               tNeighbor = mParent->get_neighbor( 7 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 7 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 7 of parent
                       mNeighbors[  7 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 7 of parent
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 8 of parent
               tNeighbor = mParent->get_neighbor( 8 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 8 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 8 of parent
                       mNeighbors[  8 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 8 of parent
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 12 of parent
               tNeighbor = mParent->get_neighbor( 12 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 12 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 12 of parent
                       mNeighbors[ 12 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 12 of parent
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 20 of parent
               tNeighbor = mParent->get_neighbor( 20 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 20 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 20 of parent
                       mNeighbors[ 20 ] = tNeighbor->get_child( 4 );
                   }
                   else
                   {
                       // link to neighbor 20 of parent
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               break;
           }
           case( 4 ) :
           {
               // link to siblings
               mNeighbors[  1 ] = mParent->get_child(  5 );
               mNeighbors[  2 ] = mParent->get_child(  6 );
               mNeighbors[  4 ] = mParent->get_child(  0 );
               mNeighbors[  7 ] = mParent->get_child(  1 );
               mNeighbors[  8 ] = mParent->get_child(  2 );
               mNeighbors[ 12 ] = mParent->get_child(  7 );
               mNeighbors[ 20 ] = mParent->get_child(  3 );

               // get neighbor 0 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 0 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor->get_child( 6 );
                       mNeighbors[  6 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 11 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor;
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 3 of parent
               tNeighbor = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 3 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor->get_child( 5 );
                       mNeighbors[  9 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 13 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor;
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 5 of parent
               tNeighbor = mParent->get_neighbor( 5 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 5 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 15 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 16 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor;
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 10 of parent
               tNeighbor = mParent->get_neighbor( 10 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 10 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 10 of parent
                       mNeighbors[ 10 ] = tNeighbor->get_child( 7 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 10 of parent
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 14 of parent
               tNeighbor = mParent->get_neighbor( 14 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 14 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 14 of parent
                       mNeighbors[ 14 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 14 of parent
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 17 of parent
               tNeighbor = mParent->get_neighbor( 17 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 17 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 17 of parent
                       mNeighbors[ 17 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 17 of parent
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 22 of parent
               tNeighbor = mParent->get_neighbor( 22 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 22 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 22 of parent
                       mNeighbors[ 22 ] = tNeighbor->get_child( 3 );
                   }
                   else
                   {
                       // link to neighbor 22 of parent
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               break;
           }
           case( 5 ) :
           {

               // link to siblings
               mNeighbors[  2 ] = mParent->get_child(  7 );
               mNeighbors[  3 ] = mParent->get_child(  4 );
               mNeighbors[  4 ] = mParent->get_child(  1 );
               mNeighbors[  8 ] = mParent->get_child(  3 );
               mNeighbors[  9 ] = mParent->get_child(  0 );
               mNeighbors[ 13 ] = mParent->get_child(  6 );
               mNeighbors[ 21 ] = mParent->get_child(  2 );

               // get neighbor 0 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 0 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 0 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor->get_child( 7 );
                       mNeighbors[  6 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 10 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 0 of parent
                       mNeighbors[  0 ] = tNeighbor;
                       mNeighbors[  6 ] = tNeighbor;
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 1 of parent
               tNeighbor = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 1 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor->get_child( 4 );
                       mNeighbors[  7 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 12 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor;
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 5 of parent
               tNeighbor = mParent->get_neighbor( 5 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 5 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 16 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 17 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor;
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 11 of parent
               tNeighbor = mParent->get_neighbor( 11 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 11 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 11 of parent
                       mNeighbors[ 11 ] = tNeighbor->get_child( 6 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 11 of parent
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 14 of parent
               tNeighbor = mParent->get_neighbor( 14 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 14 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 14 of parent
                       mNeighbors[ 14 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 14 of parent
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 15 of parent
               tNeighbor = mParent->get_neighbor( 15 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 15 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 15 of parent
                       mNeighbors[ 15 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 15 of parent
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 23 of parent
               tNeighbor = mParent->get_neighbor( 23 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 23 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 23 of parent
                       mNeighbors[ 23 ] = tNeighbor->get_child( 2 );
                   }
                   else
                   {
                       // link to neighbor 23 of parent
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }
               break;
           }
           case( 6 ) :
           {
               // link to siblings
               mNeighbors[  0 ] = mParent->get_child(  4 );
               mNeighbors[  1 ] = mParent->get_child(  7 );
               mNeighbors[  4 ] = mParent->get_child(  2 );
               mNeighbors[  6 ] = mParent->get_child(  0 );
               mNeighbors[  7 ] = mParent->get_child(  3 );
               mNeighbors[ 11 ] = mParent->get_child(  5 );
               mNeighbors[ 19 ] = mParent->get_child(  1 );

               // get neighbor 2 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 2 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor->get_child( 4 );
                       mNeighbors[  8 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 12 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor;
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 3 of parent
               tNeighbor = mParent->get_neighbor( 3 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 3 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor->get_child( 7 );
                       mNeighbors[  9 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 10 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 18 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 3 of parent
                       mNeighbors[  3 ] = tNeighbor;
                       mNeighbors[  9 ] = tNeighbor;
                       mNeighbors[ 10 ] = tNeighbor;
                       mNeighbors[ 18 ] = tNeighbor;
                   }
               }

               // get neighbor 5 of parent
               tNeighbor = mParent->get_neighbor( 5 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 5 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 14 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 15 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor;
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 13 of parent
               tNeighbor = mParent->get_neighbor( 13 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 13 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 13 of parent
                       mNeighbors[ 13 ] = tNeighbor->get_child( 5 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 13 of parent
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 16 of parent
               tNeighbor = mParent->get_neighbor( 16 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 16 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 16 of parent
                       mNeighbors[ 16 ] = tNeighbor->get_child( 0 );
                       mNeighbors[ 24 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 16 of parent
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }

               // get neighbor 17 of parent
               tNeighbor = mParent->get_neighbor( 17 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 17 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 17 of parent
                       mNeighbors[ 17 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 17 of parent
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 25 of parent
               tNeighbor = mParent->get_neighbor( 25 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 25 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 25 of parent
                       mNeighbors[ 25 ] = tNeighbor->get_child( 1 );
                   }
                   else
                   {
                       // link to neighbor 25 of parent
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }
               break;
           }
           case( 7 ) :
           {
               // link to siblings
               mNeighbors[  0 ] = mParent->get_child(  5 );
               mNeighbors[  3 ] = mParent->get_child(  6 );
               mNeighbors[  4 ] = mParent->get_child(  3 );
               mNeighbors[  6 ] = mParent->get_child(  1 );
               mNeighbors[  9 ] = mParent->get_child(  2 );
               mNeighbors[ 10 ] = mParent->get_child(  4 );
               mNeighbors[ 18 ] = mParent->get_child(  0 );

               // get neighbor 1 of parent
               Background_Element_Base* tNeighbor
                   = mParent->get_neighbor( 1 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 1 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor->get_child( 6 );
                       mNeighbors[  7 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 11 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 19 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 1 of parent
                       mNeighbors[  1 ] = tNeighbor;
                       mNeighbors[  7 ] = tNeighbor;
                       mNeighbors[ 11 ] = tNeighbor;
                       mNeighbors[ 19 ] = tNeighbor;
                   }
               }

               // get neighbor 2 of parent
               tNeighbor = mParent->get_neighbor( 2 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 2 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor->get_child( 5 );
                       mNeighbors[  8 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 13 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 21 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 2 of parent
                       mNeighbors[  2 ] = tNeighbor;
                       mNeighbors[  8 ] = tNeighbor;
                       mNeighbors[ 13 ] = tNeighbor;
                       mNeighbors[ 21 ] = tNeighbor;
                   }
               }

               // get neighbor 5 of parent
               tNeighbor = mParent->get_neighbor( 5 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 5 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor->get_child( 3 );
                       mNeighbors[ 14 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 17 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 22 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 5 of parent
                       mNeighbors[  5 ] = tNeighbor;
                       mNeighbors[ 14 ] = tNeighbor;
                       mNeighbors[ 17 ] = tNeighbor;
                       mNeighbors[ 22 ] = tNeighbor;
                   }
               }

               // get neighbor 12 of parent
               tNeighbor = mParent->get_neighbor( 12 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 12 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 12 of parent
                       mNeighbors[ 12 ] = tNeighbor->get_child( 4 );
                       mNeighbors[ 20 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 12 of parent
                       mNeighbors[ 12 ] = tNeighbor;
                       mNeighbors[ 20 ] = tNeighbor;
                   }
               }

               // get neighbor 15 of parent
               tNeighbor = mParent->get_neighbor( 15 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 15 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 15 of parent
                       mNeighbors[ 15 ] = tNeighbor->get_child( 2 );
                       mNeighbors[ 23 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 15 of parent
                       mNeighbors[ 15 ] = tNeighbor;
                       mNeighbors[ 23 ] = tNeighbor;
                   }
               }

               // get neighbor 16 of parent
               tNeighbor = mParent->get_neighbor( 16 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 16 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 16 of parent
                       mNeighbors[ 16 ] = tNeighbor->get_child( 1 );
                       mNeighbors[ 25 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 16 of parent
                       mNeighbors[ 16 ] = tNeighbor;
                       mNeighbors[ 25 ] = tNeighbor;
                   }
               }

               // get neighbor 24 of parent
               tNeighbor = mParent->get_neighbor( 24 );

               // test if neighbor exists
               if ( tNeighbor != nullptr )
               {
                   // test if neighbor 24 has children
                   if ( tNeighbor->has_children( aPattern ) )
                   {
                       // link to children of neighbor 24 of parent
                       mNeighbors[ 24 ] = tNeighbor->get_child( 0 );
                   }
                   else
                   {
                       // link to neighbor 24 of parent
                       mNeighbors[ 24 ] = tNeighbor;
                   }
               }
               break;
           }
       }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::reset_flags_of_facets()
    {
        MORIS_ERROR( false, "Don't know how to reset face flags");
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::reset_flags_of_facets()
    {
        mFacets[ 0 ]->unflag();
        mFacets[ 1 ]->unflag();
        mFacets[ 2 ]->unflag();
        mFacets[ 3 ]->unflag();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::reset_flags_of_facets()
    {
        mFacets[ 0 ]->unflag();
        mFacets[ 1 ]->unflag();
        mFacets[ 2 ]->unflag();
        mFacets[ 3 ]->unflag();
        mFacets[ 4 ]->unflag();
        mFacets[ 5 ]->unflag();
    }

//--------------------------------------------------------------------------------

    /**
     * returns an edge of the background element ( 3D only )
     */
    template< uint N >
    Background_Edge * Background_Element< N >::get_edge( uint aIndex )
    {
        MORIS_ERROR( false, "get_edge() is only available for the 3D element" );
        return nullptr;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    Background_Edge * Background_Element< 3 >::get_edge( uint aIndex )
    {
        return mEdges[ aIndex ];
    }

//--------------------------------------------------------------------------------

    /**
     * inserts an edge into the element ( 3D only )
     */
    template< uint N >
    void Background_Element< N >::insert_edge( Background_Edge * aEdge, uint aIndex )
    {
        MORIS_ERROR( false, "insert_edge() is only available for the 3D element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::insert_edge( Background_Edge * aEdge, uint aIndex )
    {
        MORIS_ASSERT( mEdges[ aIndex ] == nullptr, "tried to overwrite edge" );

        if( aEdge != nullptr )
        {
            mEdges[ aIndex ] = aEdge;
            aEdge->insert_element( this, aIndex );
        }
    }

//--------------------------------------------------------------------------------

    /**
     * inserts an edge into the element ( 3D only )
     */
    template< uint N >
    void Background_Element< N >::create_edges()
    {
        MORIS_ERROR( false, "create_edges() is only available for the 3D element" );
    }

//--------------------------------------------------------------------------------

    /**
     * resets the edge flags
     */
    template< uint N >
    void Background_Element< N >::reset_flags_of_edges()
    {
        MORIS_ERROR( false, "reset_flags_of_edges() is only available for the 3D element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::reset_flags_of_edges()
    {
        for( uint k=0; k<12; ++k )
        {
            mEdges [ k ]->unflag();
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_1(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_1() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_2(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_2() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_3(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_3() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_4(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_4() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_5(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_5() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::get_number_of_active_descendants_on_side_6(
            uint aPattern,
                  luint & aCount )
    {
        MORIS_ERROR( false, "get_number_of_active_descendants_on_side_6() not available for this element" );
    }

//--------------------------------------------------------------------------------

    template<>
    inline
    void Background_Element< 2 >::get_number_of_active_descendants_on_side_1(
            uint aPattern,
                  luint & aCount )
    {
       if( this->is_active( aPattern ) )
       {
           // add self to list
           ++aCount;
       }
       else
       {
           // add children 0 and 1 to list
           mChildren[ 0 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
           mChildren[ 1 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
       }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::get_number_of_active_descendants_on_side_2(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 1 and 3 to list
            mChildren[ 1 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::get_number_of_active_descendants_on_side_3(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 1 and 3 to list
            mChildren[ 2 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
        }
     }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::get_number_of_active_descendants_on_side_4(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 1 and 3 to list
            mChildren[ 0 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
            mChildren[ 2 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
        }
   }

//--------------------------------------------------------------------------------

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_1(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 0 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
            mChildren[ 1 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
            mChildren[ 4 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
            mChildren[ 5 ]->get_number_of_active_descendants_on_side_1( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_2(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 1 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
            mChildren[ 5 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
            mChildren[ 7 ]->get_number_of_active_descendants_on_side_2( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_3(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 2 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
            mChildren[ 6 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
            mChildren[ 7 ]->get_number_of_active_descendants_on_side_3( aPattern, aCount );
        }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_4(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 0 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
            mChildren[ 2 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
            mChildren[ 4 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
            mChildren[ 6 ]->get_number_of_active_descendants_on_side_4( aPattern, aCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_5(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 0 ]->get_number_of_active_descendants_on_side_5( aPattern, aCount );
            mChildren[ 1 ]->get_number_of_active_descendants_on_side_5( aPattern, aCount );
            mChildren[ 2 ]->get_number_of_active_descendants_on_side_5( aPattern, aCount );
            mChildren[ 3 ]->get_number_of_active_descendants_on_side_5( aPattern, aCount );
        }
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::get_number_of_active_descendants_on_side_6(
            uint aPattern,
                  luint & aCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            ++aCount;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 4 ]->get_number_of_active_descendants_on_side_6( aPattern, aCount );
            mChildren[ 5 ]->get_number_of_active_descendants_on_side_6( aPattern, aCount );
            mChildren[ 6 ]->get_number_of_active_descendants_on_side_6( aPattern, aCount );
            mChildren[ 7 ]->get_number_of_active_descendants_on_side_6( aPattern, aCount );
        }
    }

//--------------------------------------------------------------------------------

    template< uint N >
    void Background_Element< N >::collect_active_descendants_on_side_1(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_1() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::collect_active_descendants_on_side_2(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_2() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::collect_active_descendants_on_side_3(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_3() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::collect_active_descendants_on_side_4(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_4() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    void Background_Element< N >::collect_active_descendants_on_side_5(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_5() not available for this element" );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template< uint N >
    inline
    void Background_Element< N >::collect_active_descendants_on_side_6(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        MORIS_ERROR( false, "collect_active_descendants_on_side_6() not available for this element" );
    }

//--------------------------------------------------------------------------------

    template<>
    inline
    void Background_Element< 2 >::collect_active_descendants_on_side_1(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 0 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
            mChildren[ 1 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::collect_active_descendants_on_side_2(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 1 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
            mChildren[ 3 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::collect_active_descendants_on_side_3(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 3 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
            mChildren[ 2 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 2 >::collect_active_descendants_on_side_4(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 2 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
            mChildren[ 0 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
        }
    }

//--------------------------------------------------------------------------------

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_1(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 0 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
            mChildren[ 1 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
            mChildren[ 5 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
            mChildren[ 4 ]->collect_active_descendants_on_side_1( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_2(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 1 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
            mChildren[ 3 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
            mChildren[ 7 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
            mChildren[ 5 ]->collect_active_descendants_on_side_2( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_3(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 3 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
            mChildren[ 2 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
            mChildren[ 6 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
            mChildren[ 7 ]->collect_active_descendants_on_side_3( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_4(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 2 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
            mChildren[ 0 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
            mChildren[ 4 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
            mChildren[ 6 ]->collect_active_descendants_on_side_4( aPattern, aElementList, aElementCount );
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_5(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 2 ]->collect_active_descendants_on_side_5( aPattern, aElementList, aElementCount );
            mChildren[ 3 ]->collect_active_descendants_on_side_5( aPattern, aElementList, aElementCount );
            mChildren[ 1 ]->collect_active_descendants_on_side_5( aPattern, aElementList, aElementCount );
            mChildren[ 0 ]->collect_active_descendants_on_side_5( aPattern, aElementList, aElementCount );
        }
     }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline
    void Background_Element< 3 >::collect_active_descendants_on_side_6(
            uint aPattern,
            Cell< Background_Element_Base* > & aElementList,
            luint                            & aElementCount )
    {
        if( this->is_active( aPattern ) )
        {
            // add self to list
            aElementList( aElementCount++ ) = this;
        }
        else
        {
            // add children 0 and 1 to list
            mChildren[ 4 ]->collect_active_descendants_on_side_6( aPattern, aElementList, aElementCount );
            mChildren[ 5 ]->collect_active_descendants_on_side_6( aPattern, aElementList, aElementCount );
            mChildren[ 7 ]->collect_active_descendants_on_side_6( aPattern, aElementList, aElementCount );
            mChildren[ 6 ]->collect_active_descendants_on_side_6( aPattern, aElementList, aElementCount );
        }
    }

//--------------------------------------------------------------------------------
} /* namespace moris */

#include "fn_HMR_Background_Element_Neighbors_2D.hpp"
#include "fn_HMR_Background_Element_Neighbors_3D.hpp"
#include "fn_HMR_Background_Element_Edges_3D.hpp"

#endif /* SRC_HMR_CL_HMR_BACKGROUND_ELEMENT_HPP_ */

