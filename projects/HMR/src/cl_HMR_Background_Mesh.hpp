/*
 * cl_HMR_Background_Mesh.hpp
 *
 *  Created on: May 1, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_BACKGROUND_MESH_HPP_
#define SRC_HMR_CL_HMR_BACKGROUND_MESH_HPP_

#include "cl_HMR_Background_Edge.hpp" //HMR/src
#include "cl_HMR_Background_Element.hpp" //HMR/src
#include "cl_HMR_Background_Element_Base.hpp" //HMR/src
#include "cl_HMR_Background_Facet.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Domain.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "HMR_Globals.hpp" //HMR/src
#include "HMR_Tools.hpp" //HMR/src
#include "assert.hpp"
#include "cl_Communication_Tools.hpp" //COM/src
#include "cl_Communication_Manager.hpp" //COM/src

#include "typedefs.hpp" //COR/src
#include "cl_Cell.hpp" //CON/src

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_Matrix.hpp" //LINALG/src


namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

      /**
        * \brief Background Mesh Class which is templated against dimension.
        *        To be created by the factory.
        */
        template < uint N >
        class Background_Mesh : public Background_Mesh_Base
        {
            //! contains ijk lookup tables for the whole mesh
            const Domain< N >  mDomain;

            //! contains ijk lookup tables for the proc local subdomain
            Domain< N >  mMySubDomain;

            //! contains the length of an element
            real mElementLength[ gMaxNumberOfLevels ][ N ];

            //! contains the coordinates of first node on aura
            real mDomainOffset[ N ];

//--------------------------------------------------------------------------------
        public:
//--------------------------------------------------------------------------------
            /**
             * constructor for templated mesh class
             *
             * param[in]  aParameters   pointer to user defined settings container
             */
            Background_Mesh( const Parameters * aParameters ) : Background_Mesh_Base( aParameters ),
                                                                mDomain( aParameters->get_domain_ijk(),  mPaddingSize )
            {
                // create mesh decomposition
                this->decompose_mesh();

                // calculate element with, length and height for each level
                this->calculate_element_length();

                // calculate coordinate of first node on proc ( is in aura )
                this->calculate_domain_offset();

                // initialize first layer of elements
                this->initialize_coarsest_elements();

                // indices for elements on coarsest level within proc domain
                this->create_coarsest_frame();

                // set element properties on coarsest level // create aura and inverse aura
                // aura size = padding size
                this->finalize_coarsest_elements();

                // synchronize with other procs
                this->synchronize_coarsest_aura();

                // populates padding element container
                this->collect_coarsest_padding_elements();

                // create list of active elements
                this->collect_active_elements();

                // create list of active elements including aura
                this->collect_active_elements_including_aura();

                // calculate neighborhood
                this->collect_neighbors();

                // calculate indices for elements
                this->update_element_indices();
            }

//--------------------------------------------------------------------------------

            /**
             * Mesh destructor. Needs to destroy all element pointers on coarsest
             * level. Higher level elements are destroyed implicitly by element
             * destructor.
             */
            ~Background_Mesh()
            {
                // delete pointers in element cell
                for ( auto p:  mCoarsestElementsIncludingAura )
                {
                    delete p;
                }
            }

//--------------------------------------------------------------------------------

            /**
             * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
             *                                       * < max number of levels >
             *
             * @return         Matrix< DDLUMat > number of elements per direction on
             *                              proc, including aura
             */
            Matrix< DDLUMat > get_number_of_elements_per_direction_on_proc() const
            {
                Matrix< DDLUMat > aMat( N, gMaxNumberOfLevels );

                for( uint l=0; l<gMaxNumberOfLevels; ++l )
                {
                    for ( uint k=0; k<N; ++k )
                    {
                        aMat( k, l ) = mMySubDomain.mNumberOfElementsPerDimension[ l ][ k ];
                    }
                }
                return aMat;
            }

//--------------------------------------------------------------------------------

            /**
             * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
             *                                       * < max number of levels >
             *
             * @return         Matrix< DDLUMat > number of elements per direction
             *                              within whole mesh, including aura
             */
            Matrix< DDLUMat > get_number_of_elements_per_direction() const
            {
                Matrix< DDLUMat > aMat( N, gMaxNumberOfLevels );

                for( uint l=0; l<gMaxNumberOfLevels; ++l )
                {
                    for ( uint k=0; k<N; ++k )
                    {
                        aMat( k, l ) = mDomain.mNumberOfElementsPerDimension[ l ][ k ];
                    }
                }
                return aMat;
            }

//--------------------------------------------------------------------------------

            /**
             * Returns a Matrix< DDLUMat > containing the ijk positions of the calculation
             *                        domain on the proc
             *
             * @return Matrix< DDLUMat >
             */
            Matrix< DDLUMat > get_subdomain_ijk() const
            {
                Matrix< DDLUMat > aMat( 2, N );

                for ( uint k=0; k<N; ++k )
                {
                    aMat( 0, k ) = mMySubDomain.mFrameIJK[ 0 ][ k ][ 0 ];
                    aMat( 1, k ) = mMySubDomain.mFrameIJK[ 0 ][ k ][ 1 ];
                }

                return aMat;
            }

//--------------------------------------------------------------------------------

            /**
             * Returns the ijk-offset of domain of current proc.
             * This value is needed to transform global IDs to local ones and
             * vice versa
             *
             * @return Matrix< DDLUMat > of dimension  < number of dimensions >
             *                                  * <max number of levels>
             */
            Matrix< DDLUMat > get_subdomain_offset_of_proc()
            {
                uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

                Matrix< DDLUMat > aIJK( tNumberOfDimensions, gMaxNumberOfLevels );

                for( uint l=0; l<gMaxNumberOfLevels; ++l )
                {
                    for( uint k=0; k<tNumberOfDimensions; ++k )
                    {
                        aIJK( k, l ) =  mMySubDomain.mAuraIJK[ l ][ k ][ 0 ];
                    }
                }
                return aIJK;
            }
//--------------------------------------------------------------------------------

            /**
             * calculates the node coordinates of an element
             *
             * @param[in]   aElement    Element to be processed
             * @param[out]  aNodeCoords Matrix containing the node coordinates
             *                          ( 3 x 2^n )
             */
            void calc_corner_nodes_of_element( const Background_Element_Base   * aElement,
                                                     Matrix< DDRMat >          & aNodeCoords );

//--------------------------------------------------------------------------------

            /**
             * calculates the coordinates of the center of the element
             *
             * @param[in]   aElement    Element to be processed
             * @param[out]  aNodeCoords Matrix containing the node coordinates
             *                          ( 3 x 1 )
             */
            void calc_center_of_element( const Background_Element_Base   * aElement,
                                               Matrix< DDRMat >          & aNodeCoords );

//--------------------------------------------------------------------------------

            /**
             * returns a pointer to an element on level Zero
             *
             * @param[in] aI proc local i-position of element
             * @param[in] aJ proc local j-position of element
             *
             * @return pointer to element on top level
             *
             */
            Background_Element_Base * get_coarsest_element_by_ij( const luint & aI,
                                                                  const luint & aJ )
            {
                return mCoarsestElementsIncludingAura( this->calc_subdomain_id_of_element( 0, aI, aJ ) );
            }

//--------------------------------------------------------------------------------

            /**
             * returns a pointer to an element on level Zero
             *
             * @param[in] aI proc local i-position of element
             * @param[in] aJ proc local j-position of element
             *
             * @return pointer to element on top level
             *
             */
            Background_Element_Base * get_coarsest_element_by_ijk( const luint & aI,
                                                                   const luint & aJ,
                                                                   const luint & aK )
            {
                return mCoarsestElementsIncludingAura( this->calc_subdomain_id_of_element( 0, aI, aJ, aK ) );
            }

//--------------------------------------------------------------------------------

            /**
             * returns the offset of the current proc
             *
             * @return Matrix< DDRMat >
             */
            Matrix< DDRMat > get_domain_offset()
            {
                Matrix< DDRMat > aMat( N, 1 );
                for( uint k=0; k<N; ++k )
                {
                    aMat( k ) = mDomainOffset[ k ];
                }
                return aMat;
            }

//--------------------------------------------------------------------------------
        protected:
//--------------------------------------------------------------------------------

            /**
             * Calculates a local ID from a given level
             * and global ID
             *
             * @param[in]    aLevel   level of element to be investigated
             * @param[in]    aID      global ID of element to be investigated
             *
             * @return       luint    local ID on submesh of proc
             */
            luint calc_subdomain_id_from_global_id( const uint  & aLevel,
                                                    const luint & aID) const;

//--------------------------------------------------------------------------------

            /**
             * Unets element active flag, sets element flag,
             * initializes children, and removes refinement queue flag
             *
             * param[in]    aElement  element to be refined
             * param[in]    aKeepState  Indicates
             *
             * @return      vpid
             */
            void refine_element( Background_Element_Base* aElement, const bool aKeepState );

//--------------------------------------------------------------------------------

            /**
             *
             * Protected function that checks neighbors and inserts them into
             * aElement. This function is to be used for elements on level zero
             * only.
             *
             * @param[inout]  aElement    element to be processed
             *
             * @return void
             *
             */
            void collect_neighbors_on_level_zero( );

//--------------------------------------------------------------------------------

            /**
             * calculate global ID from level and from i-position ( 1D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             *
             * @return luint global ID of element
             *
             */
            luint calc_domain_id_of_element( const uint  & aLevel,
                                             const luint & aI ) const ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculate global ID from level and from ij-position ( 2D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             * @param[in] aJ     j-index of element
             *
             * @return luint global ID of element
             *
             */
            luint calc_domain_id_of_element( const uint  & aLevel,
                                             const luint & aI,
                                             const luint & aJ ) const ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculate global ID from level and from ijk-position ( 3D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             * @param[in] aJ     j-index of element
             * @param[in] aK     k-index of element
             *
             * @return luint global ID of element
             *
             */
            luint calc_domain_id_of_element(
                    const uint  & aLevel,
                    const luint & aI,
                    const luint & aJ,
                    const luint & aK ) const ;
//--------------------------------------------------------------------------------

            /**
             * subroutine for collect_side_set that collects elements on coarsest
             * level for a side
             */
            void collect_coarsest_elements_on_side(
                    const uint                       & aSideOrdinal,
                    Cell< Background_Element_Base* > & aCoarsestElementsOnSide );

//--------------------------------------------------------------------------------
        private:
//--------------------------------------------------------------------------------

            /**
             * calculates the element length lookup table
             *
             * @return void
             */
            void calculate_element_length()
            {
                // get domain dimensions from settings
                Matrix< DDRMat > tDomainDimensions = mParameters->get_domain_dimensions();

                // get number of elements on coarsest level from settings
                Matrix< DDLUMat > tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

                // calculate width for first level
                for( uint k=0; k<N; ++k )
                {
                    mElementLength[ 0 ][ k ] = tDomainDimensions( k ) / ( ( real ) ( tNumberOfElements( k ) ) );
                }

                // loop over all higher levels
                for( uint l=1; l<gMaxNumberOfLevels; ++l )
                {
                    // loop over all dimensions
                    for( uint k=0; k<N; ++k )
                    {
                        // calculate length of element
                        mElementLength[ l ][ k ] = 0.5*mElementLength[ l-1 ][ k ];
                    }
                }
            }

//--------------------------------------------------------------------------------

            /**
             * calculates the coordinate of first node on proc
             * ( including aura )
             *
             * @return void
             */
            void calculate_domain_offset()
            {
                // get domain offset
                Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

                // get padding size
                real tPaddingSize = ( real ) mParameters->get_padding_size();

                // subtract padding size from offset
                for( uint k=0; k<N; ++k )
                {
                    mDomainOffset[ k ] = tParametersOffset( k ) - tPaddingSize * mElementLength[ 0 ][ k ];
                }
            }

//--------------------------------------------------------------------------------
            /**
             * Calculates ijk lookup table for current proc.
             *
             * @return  void
             */
            void decompose_mesh()
            {
                // print output info
                if ( mParameters->is_verbose() && par_rank() == 0 )
                {
                    std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                    if ( par_size() == 1 )
                    {
                        std::fprintf( stdout, "  decomposing mesh over %u proc\n", ( unsigned int ) par_size() ) ;
                    }
                    else
                    {
                        std::fprintf( stdout, "  decomposing mesh over %u procs\n", ( unsigned int ) par_size() ) ;
                    }
                    std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                    std::fprintf( stdout, "\n" );
                }

                if ( mParameters->is_verbose() )
                {
                    // wait until all procs are here ( because of output )
                    barrier();
                }

                // calculate neighbors
                create_proc_cart( N ,
                                  mProcDims,
                                  mMyProcCoords,
                                  mMyProcNeighbors );

                // get number of elements per dimension from settings
                auto tNumberOfElementsPerDimension = mParameters->get_number_of_elements_per_dimension();

                // calculate number of elements per dimension
                Matrix< DDLUMat > tNumberOfElementsPerDimensionOnProc( N, 1 );
                for ( uint k=0; k<N; ++k )
                {
                    tNumberOfElementsPerDimensionOnProc( k ) = tNumberOfElementsPerDimension ( k ) /  mProcDims( k ) ;

                    if( par_rank() == 0 )
                    {
                        // make sure that cart size is OK
                        if ( tNumberOfElementsPerDimension( k ) % mProcDims( k ) != 0 )
                        {
                            MORIS_ASSERT( 0, "proc size incompatible to defined elements per dimension" );
                        }
                    }
                }

                // calculate decomposition domain
                // set owned and shared limits
                Matrix< DDLUMat > tDomainIJK( 2, N );

                for ( uint k=0; k<N; ++k )
                {
                    tDomainIJK( 0, k ) =   mPaddingSize + tNumberOfElementsPerDimensionOnProc ( k ) * mMyProcCoords ( k );

                    tDomainIJK( 1, k ) =   tDomainIJK( 0, k ) + tNumberOfElementsPerDimensionOnProc ( k ) - 1;
                }

                // create proc domain
                Domain< N > tMySubDomain( tDomainIJK, mPaddingSize );

                // move domain object to member
                mMySubDomain = std::move( tMySubDomain );

                // print proc area
                if ( mParameters->is_verbose() )
                {
                    moris_id tMyRank = par_rank() ;

                    uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

                    std::string tString = "  proc " + std::to_string( tMyRank );

                    // add dots to the string for pretty output
                    if ( tMyRank < 10 )
                    {
                        tString +=" ... :" ;
                    }
                    else if ( tMyRank < 100 )
                    {
                        tString +=" .. :" ;
                    }
                    else if ( tMyRank < 1000 )
                    {
                        tString +=" . :" ;
                    }
                    else if ( tMyRank < 10000 )
                    {
                        tString +="  :" ;
                    }
                    else
                    {
                        tString +=" :" ;
                    }

                    if ( tNumberOfDimensions == 1 )
                    {
                        std::fprintf( stdout, "%s owns i domain ", tString.c_str() ) ;
                    }
                    else if ( tNumberOfDimensions == 2 )
                    {
                        std::fprintf( stdout, "%s owns i-j domain ", tString.c_str()  ) ;
                    }
                    else if( tNumberOfDimensions == 3 )
                    {
                        std::fprintf( stdout, "%s owns i-j-k domain ", tString.c_str() ) ;
                    }

                    // print ijk domain
                    std::fprintf( stdout, "%lu-%lu",
                            ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 0 ],
                            ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ 0 ][ 1 ] );
                    for ( uint k=1; k<N; ++k )
                    {
                        std::fprintf( stdout, ", %lu-%lu ",
                                ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ k ][ 0 ],
                                ( long unsigned int ) mMySubDomain.mDomainIJK[ 0 ][ k ][ 1 ] );
                    }
                    std::fprintf( stdout, "\n\n" );
                }

                // test if settings are OK
                if( par_rank() == 0 &&  par_size() != 1 )
                {
                    // get number of dimensions from settings
                    uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

                    Matrix< DDLUMat > tProcSplit( tNumberOfDimensions, 1 );

                    bool tError = false;

                    // loop over all dimensions and copy elements
                    for( uint k=0; k<tNumberOfDimensions; ++k )
                    {
                        tProcSplit( k ) = mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ k ] - 2*mPaddingSize;

                        if ( tProcSplit( k ) < mPaddingSize )
                        {
                            tError = true;
                        }
                    }

                    // test if error occurred
                    if ( tError )
                    {
                        std::fprintf( stdout, "     ERROR: Mesh too coarse for selected split and padding size.\n" );
                        for( uint k=0; k<tNumberOfDimensions; ++k )
                        {
                            if (  tProcSplit( k ) < mPaddingSize )
                            {
                                // print in correct grammar
                                if ( tProcSplit( k ) == 1 )
                                {
                                    std::fprintf( stdout,
                                            "            In dimension %u, each proc domain gets 1 element, need at least %lu.\n",
                                            ( unsigned int ) k,
                                            ( long unsigned int ) mPaddingSize );
                                }
                                else
                                {
                                    std::fprintf( stdout,
                                            "            In dimension %u, each proc domain gets %lu elements, need at least %lu.\n",
                                            ( unsigned int ) k,
                                            ( long unsigned int ) tProcSplit( k ),
                                            ( long unsigned int ) mPaddingSize );
                                }
                            }
                        }
                        std::fprintf( stdout, "\n" );
                        exit( -1 );
                    }
                }
            }

//--------------------------------------------------------------------------------
            /**
             * In this subroutine, the coarsest layer of elements on the proc,
             * including the aura, is generated. The element pointers are created
             * and stored in a fixed size cell, mCoarsestElementsIncludingAura.
             * New elements are created with the assumption that they are active,
             * the ID of the proc owning that element is not known yet, and set to
             * MORIS_UINT_MAX.
             *
             * @return  void
             */
            void initialize_coarsest_elements();

//--------------------------------------------------------------------------------
            /**
             * This function loops over all elements on the coarsest level,
             * identifies padding elements and calculates element ownership. It
             * also identifies the elements on the coarsest level that belong to
             * theaura. In fact, there are two auras. The normal aura around the
             * proc domain, which contains the elements that belong to the neighbor
             * proc and are shared with the current proc, and the inverse aura,
             * which contains the elements that belong to the current proc,
             * and are shared with the neighbor. In most cases, padding elements are
             * not considered for the aura.
             *
             * @return  void
             */
            void finalize_coarsest_elements();

//--------------------------------------------------------------------------------
            /**
             * Sometimes, it is more helpful to only access the elements that are
             * within the calculation domain of the proc, excluding the aura.
             * This function creates a list of local indices of these elements.
             *
             * @return  void
             */
            void create_coarsest_frame();

//--------------------------------------------------------------------------------

            /**
             * Calculates number of elements on first level.
             * This function is only needed once.
             * Therefore, template specialization is not required.
             *
             * @return   Matrix< DDLUMat > of dumension < number of dimensions>
             *                       containing number of elements per direction
             *                       on coarsest proc, including aura
             */
            Matrix< DDLUMat > get_number_of_subdomain_elements_per_direction_on_level_zero()
            {
                Matrix< DDLUMat > aNumberOfElements( N, 1 );
                for( uint k=0; k<N; ++k )
                {
                    aNumberOfElements( k ) =  mMySubDomain.mNumberOfElementsPerDimension[ 0 ][ k ];
                }
                return aNumberOfElements;
            }

//--------------------------------------------------------------------------------

            /**
             * Internal function called from initialize_coarsest_elements()/
             * Takes the pointer of the element and puts it into the correct
             * position in mCoarsestElementsIncludingAura
             *
             *
             */
            void insert_zero_level_element( const luint                   & aPosition,
                                                  Background_Element_Base * aElement )
            {
                mCoarsestElementsIncludingAura( aPosition ) = aElement ;
            }

//-------------------------------------------------------------------------------

            /**
             * calculate subdomain ID from level and from i-position ( 1D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             *
             * @return luint subdomain ID of element
             *
             */
            luint calc_subdomain_id_of_element( const uint  & aLevel,
                                                const luint & aI ) const ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculate subdomain ID from level and from ij-position ( 2D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             * @param[in] aJ     j-index of element
             *
             * @return luint subdomain ID of element
             *
             */
            luint calc_subdomain_id_of_element( const uint  & aLevel,
                                                const luint & aI,
                                                const luint & aJ ) const ;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            /**
             * calculate subdomain ID from level and from ijk-position ( 3D case )
             *
             * @param[in] aLevel level of element
             * @param[in] aI     i-index of element
             * @param[in] aJ     j-index of element
             * @param[in] aK     k-index of element
             *
             * @return luint subdomain ID of element
             *
             */
            luint calc_subdomain_id_of_element( const uint  & aLevel,
                                                const luint & aI,
                                                const luint & aJ,
                                                const luint & aK ) const;

//--------------------------------------------------------------------------------
            /**
             * Internal function called during element refinement.
             * Calculates IDs and Subdomain IDs for given ijk-positions
             *
             * @param[in]  aLevel         level of child elements to be considered
             * @param[in]  aIJK           matrix of dimension <number of dimensions>
             *                           * <number of children>
             * @param[out] aIDs           global IDs for each child element
             * @param[out] aSubdomainIDs  proc local IDs for each child element
             *
             * @return void
             */
            void calc_element_ids( const uint              & aLevel,
                                   const Matrix< DDLUMat > & aIJK,
                                         Matrix< DDLUMat > & aIDs ) const;

//--------------------------------------------------------------------------------

            void check_queued_element_for_padding( Background_Element_Base * aElement  )
            {
                // only do something if this element belongs to me
                if ( aElement->get_owner() == mMyRank )
                {
                    // get local ijk position of element
                    const luint * tIJK = aElement->get_ijk();

                    // get level of element
                    uint tLevel = aElement->get_level();

                    // perform aura check
                    bool tIsPaddingCandidate = false;

                    // loop over all dimensions
                    for( uint k=0; k<N; ++k )
                    {
                        // claculate global coordinate of element
                        luint tI = tIJK[ k ] + mMySubDomain.mAuraIJK[ tLevel ][ k ][ 0 ];

                        // test if element is candidate for padding test
                        tIsPaddingCandidate =
                                tIsPaddingCandidate ||
                               ( ( tI < mDomain.mDomainIJK[ tLevel ][ k ][ 0 ] + mPaddingRefinement )
                            ||   ( tI > mDomain.mDomainIJK[ tLevel ][ k ][ 1 ] - mPaddingRefinement ) );
                    }

                    if ( tIsPaddingCandidate )
                    {
                        // get neighbors from samel level
                        Cell< Background_Element_Base* > tNeighbors;
                        aElement->get_neighbors_from_same_level( mPaddingRefinement,
                                tNeighbors );

                        // loop over all neighbors
                        for ( auto tNeighbor : tNeighbors )
                        {
                            // test if neighbor os padding
                            if ( tNeighbor->is_padding() && ! tNeighbor->is_queued_for_refinement() )
                            {
                                // flag padding element for refinement
                                tNeighbor->put_on_refinement_queue();
                            }
                        }
                    }
                }
            }

//--------------------------------------------------------------------------------
        }; /* Background_Mesh */

//--------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::initialize_coarsest_elements()
        {
            MORIS_ERROR( false, "Do not know how initialize elements\n");
        }

//--------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::finalize_coarsest_elements()
        {
            MORIS_ERROR( false, "Don't know how to finalize coarsest level.");
        }

//--------------------------------------------------------------------------------

        template < uint N >
        luint Background_Mesh< N >::calc_domain_id_of_element(
                const uint  & aLevel,
                const luint & aI ) const
        {
            MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called.");
            return 0;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < uint N >
        luint Background_Mesh< N >::calc_domain_id_of_element(
                const uint  & aLevel,
                const luint & aI,
                const luint & aJ ) const
        {
            MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called.");
            return 0;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        template < uint N >
        luint Background_Mesh< N >::calc_domain_id_of_element(
                const uint  & aLevel,
                const luint & aI,
                const luint & aJ,
                const luint & aK ) const
        {
            MORIS_ERROR( false, "wrong function calc_domain_id_of_element() called.");
            return 0;
        }

//--------------------------------------------------------------------------------

        template < uint N >
        luint Background_Mesh< N >::calc_subdomain_id_of_element( const uint  & aLevel,
                                                                  const luint & aI ) const
        {
            MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called.");
            return 0;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < uint N >
        luint Background_Mesh< N >::calc_subdomain_id_of_element( const uint  & aLevel,
                                                                  const luint & aI,
                                                                  const luint & aJ ) const
        {
            MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called.");
            return 0;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template < uint N >
        luint Background_Mesh< N >::calc_subdomain_id_of_element( const uint  & aLevel,
                                                                  const luint & aI,
                                                                  const luint & aJ,
                                                                  const luint & aK ) const
        {
            MORIS_ERROR( false, "wrong function calc_subdomain_id_of_element() called.");
            return 0;
        }

//--------------------------------------------------------------------------------

         template < uint N >
         void Background_Mesh< N >::calc_element_ids(
                 const uint         & aLevel,
                 const Matrix< DDLUMat > & aIJK,
                 Matrix< DDLUMat >       & aIDs ) const
         {
             MORIS_ERROR( false, "Don't know how to calculate ids yet.");
         }

//--------------------------------------------------------------------------------

         template < uint N >
         luint Background_Mesh< N >::calc_subdomain_id_from_global_id(
                 const uint         & aLevel,
                 const luint        & aID) const
         {
             MORIS_ERROR( false, "Don't know how to calculate IDs yet.");
             return 0;
         }

//--------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::refine_element( Background_Element_Base * aElement, const bool aKeepState )
        {
            MORIS_ERROR( false, "Don't know how to refine element." );
        }

//--------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::create_coarsest_frame()
        {
            MORIS_ERROR( false, "Don't know how to create coarsest frame.");
        }

//-------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::collect_neighbors_on_level_zero()
        {
            MORIS_ERROR( false, "Don't know how to collect_neighbors_on_level_zero.");
        }


//-------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::calc_corner_nodes_of_element(
                const Background_Element_Base   * aElement,
                Matrix< DDRMat >                       & aNodeCoords )
        {
            MORIS_ERROR( false,  "Do not know how calculate corner nodes\n" );
        }

//-------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::calc_center_of_element(
                const Background_Element_Base  * aElement,
                Matrix< DDRMat >                      & aNodeCoords )
        {
            MORIS_ERROR( false,  "Do not know how calculate center of element\n");
        }

//-------------------------------------------------------------------------------

        template < uint N >
        void Background_Mesh< N >::collect_coarsest_elements_on_side(
                const uint                       & aSideOrdinal,
                Cell< Background_Element_Base* > & aCoarsestElementsOnSide )
        {
            MORIS_ERROR( false,  "Do not know how to collect coarsest elements on side \n");
        }

//--------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#include "cl_HMR_Background_Mesh_2D.hpp"
#include "cl_HMR_Background_Mesh_3D.hpp"

#endif /* SRC_HMR_CL_HMR_BACKGROUND_MESH_HPP_ */
