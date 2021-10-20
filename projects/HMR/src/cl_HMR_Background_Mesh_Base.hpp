/*
 * cl_HMR_Base.hpp
 *
 *  Created on: May 7, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_BACKGROUND_MESH_BASE_HPP_
#define SRC_HMR_CL_HMR_BACKGROUND_MESH_BASE_HPP_

#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "HMR_Globals.hpp" //HMR/src
#include "assert.hpp"
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp" //LINALG/src



namespace moris
{
    namespace hmr
    {
        //--------------------------------------------------------------------------------
        /**
         * \brief Base class for templated BackgroundMesh
         *
         * The background mesh contains a hierarchical data structure, a quadtree
         * in 2D, or an octree in 3D, that manages active and refined elements of
         * the mesh. The class also contains all information on parallel mesh
         * decomposition, and manages parallel communication. class is templated
         * against the dimension, 2D or 3D. The background mesh does not contain
         * any node or basis information.
         */
        class Background_Mesh_Base
        {
            protected:
                //! ref to user defined settings
                const Parameters * mParameters;

                //! dimensions as set in settings
                const uint         mNumberOfDimensions = MORIS_UINT_MAX;

                //! max polynomial as set in settings
                const uint         mMaxPolynomial = MORIS_UINT_MAX;

                //! padding refinement radius for B-Splines
                const uint         mPaddingRefinement = MORIS_UINT_MAX;

                //! padding size defines the thickness of the aura and the invisible
                //! elements at the edge of the domain. It is the maximum of the
                //! buffer size, the maximum polynomial and filter width
                const luint        mPaddingSize = MORIS_LUINT_MAX;

                //! size of refinement buffer as set by user
                luint              mBufferSize = MORIS_LUINT_MAX;

                //! defined by the dimension. 2 in 1D, 4 in 2D, 8 in 3D
                const uint         mNumberOfChildrenPerElement = MORIS_UINT_MAX;

                //! rank of this proc
                const moris_id     mMyRank = gNoID;

                //! neighbors on same level (active or not)
                uint               mNumberOfNeighborsPerElement = MORIS_UINT_MAX;

                //! subdomain IDs of elements shared with neighbor, owned by neighbor
                Cell< Matrix< DDLUMat > > mCoarsestAura;

                //! subdomain IDs of elements shared with neighbor, owned by myself
                Cell< Matrix< DDLUMat > > mCoarsestInverseAura;

                //! defines how many procs are used in i, j and k-direction
                Matrix< DDUMat >          mProcDims;

                //! i-j-k coordinate of current proc
                Matrix< DDUMat >          mMyProcCoords;

                //! Mat containing IDs of neighbors and the proc itself
                //! if there is no neighbor, value is MORIS_UINT_MAX
                Matrix< IdMat >           mMyProcNeighbors;

                // all elements on level zero on proc including aura
                Cell< Background_Element_Base * > mCoarsestElementsIncludingAura;

                //! active elements within proc domain ( excluding aura )
                Cell< Background_Element_Base * > mActiveElements;

                //! active elements within proc domain ( including aura )
                Cell< Background_Element_Base * > mActiveElementsIncludingAura;

                //! local indices with active and refined elements in proc domain (without aura)
                Cell< Background_Element_Base * > mCoarsestElements;

                //! cell containing coarsest padding elements
                Cell< Background_Element_Base * > mCoarsestPaddingElements;

                //! list of elements to be refined
                Cell< Background_Element_Base * > mRefinementQueue;

                //! lookup table containing number of elements per level
                //! updated through count_elements
                luint  mNumberOfElementsPerLevel[ gMaxNumberOfLevels ];

                //! maximum number of levels
                //! updated through count_elements
                uint mMaxLevel = 0;

                //! pattern this mesh operates on
                uint mActivePattern = 0;

                luint mMaxElementDomainIndex = 0;

                //--------------------------------------------------------------------------------
            public:
                //--------------------------------------------------------------------------------

                /**
                 * Base class constructor. Reserves size of
                 * mCoarsestElementsIncludingAura
                 * mCoarsestAura and mCoarsestInverseAura
                 *
                 * @param[in] aParameters   ref to container of user defined settings
                 */
                Background_Mesh_Base( const Parameters * aParameters );

                //--------------------------------------------------------------------------------

                /**
                 * virtual destructor of base mesh class. Does nothing.
                 */
                virtual ~Background_Mesh_Base()
                {

                }

                //--------------------------------------------------------------------------------

                /**
                 * return a pointer to the parameter object
                 */
                const Parameters * get_parameters() const
                {
                    return mParameters;
                }

                //--------------------------------------------------------------------------------

                /**
                 * Prints element information of coarsest level.
                 * Useful for debugging.
                 *
                 * @return void
                 */
                void print_level_zero()
                {
                    for (luint e = 0; e<mCoarsestElementsIncludingAura.size(); ++e )
                    {
                        std::fprintf( stdout, "  el: %lu id: %lu a: %d  r: %d  o: %u\n",
                                e,
                                mCoarsestElementsIncludingAura( e )->get_hmr_id(),
                                mCoarsestElementsIncludingAura( e )->is_active( mActivePattern ),
                                mCoarsestElementsIncludingAura( e )->is_refined( mActivePattern ),
                                mCoarsestElementsIncludingAura( e )->get_owner() );
                    }
                }
                //--------------------------------------------------------------------------------

                /**
                 * Returns a pointer to an element in the list of active elements.
                 * Note that the index of an element potentially changes after
                 * refinement.
                 *
                 * @param[in] aIndex   index on element list
                 *                     ( might change after refinement )
                 *
                 * @return pointer to element in mActiveElements
                 *
                 */
                Background_Element_Base * get_element( const luint & aIndex )
                {
                    return mActiveElements( aIndex );
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns a pointer to an element in the list of active elements.
                 * Note that the index of an element potentially changes after
                 * refinement. ( const version )
                 *
                 * @param[in] aIndex   index on element list
                 *                     ( might change after refinement )
                 *
                 * @return pointer to element in mActiveElements
                 *
                 */
                const Background_Element_Base * get_element( const luint & aIndex ) const
                {
                    return mActiveElements( aIndex );
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns a pointer to an element in the list of active elements.
                 * Note that the index of an element potentially changes after
                 * refinement.
                 *
                 * @param[in] aIndex   index on element list
                 *                     ( might change after refinement )
                 *
                 * @return pointer to element in mActiveElements
                 *
                 */
                Background_Element_Base * get_element_including_aura( const luint & aIndex )
                {
                    return mActiveElementsIncludingAura( aIndex );
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns a pointer to an element in the list of active elements.
                 * Note that the index of an element potentially changes after
                 * refinement. ( const version )
                 *
                 * @param[in] aIndex   index on element list
                 *                     ( might change after refinement )
                 *
                 * @return pointer to element in mActiveElements
                 *
                 */
                const Background_Element_Base * get_element_including_aura( const luint & aIndex ) const
                {
                    return mActiveElementsIncludingAura( aIndex );
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns a pointer to an element in the list of active elements.
                 * Note that the index of an element potentially changes after
                 * refinement.
                 *
                 * @param[in] aIndex   index on element list
                 *                     ( might change after refinement )
                 *
                 * @return pointer to element in mActiveElementsIncludingAura
                 *
                 */
                Background_Element_Base * get_element_from_proc_domain_including_aura( const luint & aIndex )
                {
                    return mActiveElementsIncludingAura( aIndex );
                }

                //--------------------------------------------------------------------------------

                /**
                 * returns a pointer to an element on level Zero
                 *
                 * @param[in] aSubdomainID  local ID ( = local index )
                 *                          of top level element
                 *
                 * @return pointer to element on top level
                 *
                 */
                Background_Element_Base * get_coarsest_element_by_subdomain_id( const luint & aSubdomainID )
                {
                    return mCoarsestElementsIncludingAura( aSubdomainID );
                }
                //--------------------------------------------------------------------------------

                /**
                 * returns a pointer to an element on level Zero
                 *
                 * @param[in] aSubdomainIndex  Index in array on domain
                 *
                 * @return pointer to element on top level
                 *
                 */
                Background_Element_Base* get_coarsest_element_by_subdomain_index( const luint & aSubdomainIndex )
                {
                    return mCoarsestElements( aSubdomainIndex );
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
                virtual Background_Element_Base * get_coarsest_element_by_ij( const luint & aI,
                        const luint & aJ ) = 0;

                //--------------------------------------------------------------------------------

                /**
                 * returns a pointer to an element on level Zero
                 *
                 * @param[in] aI proc local i-position of element
                 * @param[in] aJ proc local j-position of element
                 * @param[in] aK proc local k-position of element
                 *
                 * @return pointer to element on top level
                 *
                 */
                virtual Background_Element_Base * get_coarsest_element_by_ijk( const luint & aI,
                        const luint & aJ,
                        const luint & aK ) = 0;

                //--------------------------------------------------------------------------------
                /**
                 * Refines the element, adds children, sets refined switch and so on
                 *
                 * @param[in] aElement      pointer to  element to be processed
                 * @param[in] aKeepState    Indicates
                 *
                 * @return void
                 *
                 */
                virtual void refine_element( Background_Element_Base * aElement, const bool aKeepState ) = 0;

                //--------------------------------------------------------------------------------

                /**
                 * collects all active elements (excluding aura) and saves their
                 * pointers in mActiveElements
                 *
                 * @return void
                 */
                void collect_active_elements();

                //--------------------------------------------------------------------------------

                /**
                 * collects all elements (including aura) on proc
                 *
                 * @return void
                 */
                void collect_all_elements( Cell< Background_Element_Base* >  & aElementList );

                //--------------------------------------------------------------------------------
                /**
                 * collects all active elements (including aura) and saves their
                 * pointers in mActiveElementsIncludingAura
                 */
                void collect_active_elements_including_aura();

                //--------------------------------------------------------------------------------

                /**
                 * output for debugging
                 *
                 * @return void
                 */
                void print_active_elements();

                //--------------------------------------------------------------------------------
                /**
                 * output for debugging
                 *
                 * @return void
                 */
                void print_active_elements_including_aura();

                //--------------------------------------------------------------------------------

                //void
                //set_buffer_size( const uint & aBuffer )
                //{
                //    mBufferSize = aBuffer;
                //}

                //--------------------------------------------------------------------------------

                /**
                 * returns a moris::mat with the local IDs (= local indices)
                 * of top level aura elements belonging to the the neighbor proc.
                 *
                 * @param[in]     aNeighborProc   neighbor number of proc
                 *                                ( not to be confused with proc rank )
                 *
                 * @return        Matrix< DDLUMat > containing local IDs (= local indices)
                 *                             of top level elements
                 *
                 */
                Matrix< DDLUMat > get_subdomain_ids_of_coarsest_aura( const uint & aNeighborProc ) const
                        {
                    return mCoarsestAura( aNeighborProc );
                        }

                //--------------------------------------------------------------------------------
                /**
                 * Returns a moris::mat with the local IDs of coarsest elements
                 * that are to be shared with the neighbor proc, owned by myself.
                 *
                 * @param[in]     aNeighborProc   neighbor number of proc
                 *                                ( not to be confused with proc rank )
                 *
                 * @return        Matrix< DDLUMat > containing local IDs
                 *
                 */
                Matrix< DDLUMat > get_subdomain_ids_of_coarsest_inverse_aura( const uint & aNeighborProc ) const
                        {
                    return mCoarsestInverseAura( aNeighborProc );
                        }
                //--------------------------------------------------------------------------------

                /**
                 * Returns a mat containing the IDs of active elements on proc.
                 * Calls collect_active_elements before, to use actual pattern.
                 * Useful for debugging.
                 *
                 * @param[out]   aElementIDs  : Matrix< DDLUMat > containing global IDs of
                 *                              active elements on proc (without aura)
                 *
                 * @return       void
                 */
                void get_active_elements_on_proc( Matrix< DDLUMat > & aElementIDs );

                //--------------------------------------------------------------------------------

                /**
                 * Returns a mat containing the IDs of active elements on proc.
                 * Calls collect_active_elements_including_aura before, to use actual pattern.
                 * Useful for debugging.
                 *
                 * @param[out]   aElementIDs  : Matrix< DDLUMat > containing global IDs of
                 *                              active elements on proc, including aura
                 *
                 * @return         void
                 *
                 */
                void get_active_elements_on_proc_including_aura( Matrix< DDLUMat > & aElementIDs );

                //--------------------------------------------------------------------------------
                /**
                 * tells how many elements are on the coarsest level (including aura)
                 *
                 * return @luint
                 */
                luint get_number_of_coarsest_elements_on_proc_including_aura()
                {
                    return mCoarsestElementsIncludingAura.size();
                }

                //--------------------------------------------------------------------------------

                /**
                 * tells how many elements are on the coarsest level (without aura)
                 *
                 * return @luint
                 */
                luint
                get_number_of_coarsest_elements_on_proc()
                {
                    return mCoarsestElements.size();
                }

                //--------------------------------------------------------------------------------

                /**
                 * returns the size of mRefinementQueue
                 *
                 * @return luint
                 */
                luint
                get_number_of_queued_elements()
                {
                    return mRefinementQueue.size();
                }

                //--------------------------------------------------------------------------------

                /**
                 * applies this refinment queue to another pattern.
                 * Pattern can change if second pattern is not refined as first pattern
                 *
                 * @param[ in ]     uint active pattern
                 */
                void apply_refinement_queue_to_pattern( const uint aPattern );

                //--------------------------------------------------------------------------------

                /**
                 * synchronizes and processes the refinement queue
                 * and updates active element table
                 *
                 *@param[ in ]     uint active pattern
                 *
                 * @return         bool telling if at least one element has been refined
                 */
                bool perform_refinement( const uint aPattern );

                //--------------------------------------------------------------------------------

                /**
                 * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
                 *                                       * < max number of levels >
                 *
                 * @return         Matrix< DDLUMat > number of elements per direction on
                 *                              proc, including aura
                 */
                virtual Matrix< DDLUMat > get_number_of_elements_per_direction_on_proc() const = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Returns a Matrix< DDLUMat > of the dimension < number of dimensions >
                 *                                       * < max number of levels >
                 *
                 * @return         Matrix< DDLUMat > number of elements per direction
                 *                              within whole mesh, including aura
                 */
                virtual Matrix< DDLUMat > get_number_of_elements_per_direction() const = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Returns a Matrix< DDLUMat > containing the ijk positions of the calculation
                 *                        domain on the proc
                 *
                 * @return Matrix< DDLUMat >
                 */
                virtual Matrix< DDLUMat > get_subdomain_ijk() const = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Returns the number of active elements in the cell mActiveElements.
                 * Needs to be called after collect_active_elements()
                 *
                 * @return  luint  size of cell mActiveElements
                 *
                 */
                luint get_number_of_active_elements_on_proc()
                {
                    MORIS_ASSERT( mActiveElements.size() > 0, "No active elements found on mesh" );
                    return mActiveElements.size();
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns the number of active elements in the cell
                 * mActiveElementsIncludingAura. Needs to be called after
                 * collect_active_elements_including_aura()
                 *
                 * @return  luint  size of cell mActiveElementsIncludingAura
                 */
                luint get_number_of_active_elements_on_proc_including_aura()
                {
                    return mActiveElementsIncludingAura.size();
                }

                //--------------------------------------------------------------------------------

                /**
                 * Returns the number of padding elements.
                 * Needs to be called after collect_active_elements()
                 *
                 * @return  luint  size of cell mActiveElements
                 *
                 */
                luint get_number_of_padding_elements_on_proc()
                {
                    MORIS_ASSERT( mCoarsestPaddingElements.size() > 0, "No padding elements found on mesh" );

                    luint tNumberOfCoarsestPaddingElements = mCoarsestPaddingElements.size();

                    // initialize padding counter
                    luint tPaddingCount = 0;

                    // loop over all coarsest padding elements
                    for( luint k = 0; k < tNumberOfCoarsestPaddingElements; ++k )
                    {
                        // count descendants
                        mCoarsestPaddingElements( k )->get_number_of_descendants( tPaddingCount );
                    }

                    return tPaddingCount;
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
                virtual Matrix< DDLUMat > get_subdomain_offset_of_proc() = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Returns cartesian ijk-coordinates of current proc.
                 *
                 * @return Matrix< DDUMat >
                 */
                const Matrix< DDUMat > & get_proc_coords() const
                        {

                    return mMyProcCoords;
                        }
                //--------------------------------------------------------------------------------

                /**
                 * Returns the number of procs per direction
                 *
                 * @return Matrix< DDUMat >
                 */
                const Matrix< DDUMat > & get_proc_dims() const
                        {
                    return mProcDims;
                        }

                //--------------------------------------------------------------------------------

                /**
                 * Protected funcition that adds neighbors to elements.
                 * Internally calls collect_active_elements_including_aura()
                 *
                 * @return void
                 */
                void collect_neighbors();

                //--------------------------------------------------------------------------------

                /**
                 * calculates a consecutive index for all elements on domain.
                 * Must be called after collect_active_elements.
                 */
                void update_element_indices();

                //--------------------------------------------------------------------------------

                /**
                 * Get maximum level of mesh on proc.
                 * Must be called after count_elements()
                 *
                 * @return uint   max level of mesh on local proc
                 */
                auto get_max_level() const -> decltype( mMaxLevel )
                {
                    return mMaxLevel;
                }

                //--------------------------------------------------------------------------------

                /**
                 * collects elements on level aLevel and writes them into a cell
                 * including elements on aura
                 *
                 * @param[in]    aLevel        level to be considered
                 * @param[out]   aElementList  cell in which pointers are written
                 *
                 * @return       void
                 */
                void collect_elements_on_level_including_aura(
                        const uint                       & aLevel,
                        Cell< Background_Element_Base* > & aElementList );

                //--------------------------------------------------------------------------------

                /**
                 * collects elements on level aLevel and writes them into a cell
                 * aura is not regarded
                 *
                 * @param[in]    aLevel        level to be considered
                 * @param[out]   aElementList  cell in which pointers are written
                 *
                 * @return       void
                 */
                void collect_elements_on_level_within_proc_domain(
                        const uint                       & aLevel,
                        Cell< Background_Element_Base* > & aElementList );

                //--------------------------------------------------------------------------------

                /**
                 * Creates a list of active elements shared with a neighbor.
                 * Needed for node ownership
                 *
                 *  @param[in]   aProcNeighbor
                 *  @param[in]   aMode      0:  collect aura
                 *                          1:  collect inverse aura
                 *                          2:  collect aura and inverse aura
                 *
                 *  @param[out]  aElementList     cell with element pointers on aura
                 */
                void collect_active_elements_from_aura(
                        const uint                       & aProcNeighbor,
                        const uint                       & aMode,
                        Cell< Background_Element_Base* > & aElementList );

                //--------------------------------------------------------------------------------

                /**
                 * returns number of proc neighbors ( 1D: 3, 2D: 9, 3D: 27 )
                 *
                 * return uint
                 */
                uint get_number_of_proc_neighbors() const
                {
                    return mMyProcNeighbors.length();
                }

                //--------------------------------------------------------------------------------

                /**
                 * returns the matrix containing the proc neighbor IDs
                 *
                 * @return Matrix< DDUMat >
                 */
                const Matrix< IdMat > & get_proc_neigbors() const
                        {
                    return mMyProcNeighbors;
                        }

                //--------------------------------------------------------------------------------

                /**
                 * calculates the node coordinates of an element
                 *
                 * @param[in]   aElement    Element to be processed
                 * @param[out]  aNodeCoords Matrix containing the node coordinates
                 *                          ( 3 x 2^n )
                 */
                virtual void calc_corner_nodes_of_element(
                        const Background_Element_Base * aElement,
                        Matrix< DDRMat >              & aNodeCoords ) = 0;

                //--------------------------------------------------------------------------------

                /**
                 * calculates the coordinates of the center of the element
                 *
                 * @param[in]   aElement    Element to be processed
                 * @param[out]  aNodeCoords Matrix containing the node coordinates
                 *                          ( 3 x 1 )
                 */
                virtual void calc_center_of_element(
                        const Background_Element_Base * aElement,
                        Matrix< DDRMat >              & aNodeCoords ) = 0;

                //--------------------------------------------------------------------------------

                /**
                 * returns the number of neihgbors per element
                 *
                 * @return uint
                 */
                auto get_number_of_neighbors_per_element() const
                -> decltype( mNumberOfNeighborsPerElement )
                {
                    return mNumberOfNeighborsPerElement;
                }

                //--------------------------------------------------------------------------------

                /* Element_Base *
            decode_pedigree_path(
                    const luint       & aAncestorID
                    const Matrix< DDUMat > & aPedigreeList,
                    luint             & aCounter ); */
                Background_Element_Base * decode_pedigree_path(
                        const luint            & aAncestorID,
                        const Matrix< DDUMat > & aPedigreeList,
                        luint                  & aCounter );

                //--------------------------------------------------------------------------------

                /**
                 * creates a VTKfile of the mesh on the proc for debugging
                 */
                void save_to_vtk( const std::string & aFilePath );

                //--------------------------------------------------------------------------------

                /**
                 * returns the offset of the current proc
                 *
                 * @return Matrix< DDRMat >
                 */
                virtual Matrix< DDRMat > get_domain_offset() = 0;


                //--------------------------------------------------------------------------------

                /**
                 * tell the other procs that we need an element, and therefore
                 * its DOFs
                 */
                void synchronize_t_matrix_flags();

                // ----------------------------------------------------------------------------

                /**
                 * sets the active pattern to another value
                 */
                void set_activation_pattern( const uint & aPattern )
                {
                    if( mActivePattern != aPattern )
                    {
                        MORIS_ERROR( aPattern < gNumberOfPatterns, "Invalid Pattern index.");

                        MORIS_LOG_INFO( "%s Select activation pattern %u.",
                                proc_string().c_str(),
                                ( unsigned int ) aPattern );

                        mActivePattern = aPattern;

                        this->update_database();
                    }
                }

                // -----------------------------------------------------------------------------

                /**
                 * returns the active refinement pattern of the background mesh
                 */
                uint get_activation_pattern() const
                {
                    return mActivePattern;
                }

                // -----------------------------------------------------------------------------

                /**
                 * clones one pattern into another. Unlike copy pattern,
                 * this pattern can extend the buffer if desired.
                 */
                void clone_pattern( const uint & aSource, const uint & aTarget );

                // -----------------------------------------------------------------------------

                /**
                 * creates a union of two patterns
                 */
                void unite_patterns( const uint & aSourceA,
                        const uint & aSourceB,
                        const uint & aTarget );

                // -----------------------------------------------------------------------------

                /**
                 * creates a union of two patterns
                 */
                void unite_patterns( const moris::Cell< uint > & aSourcePattern,
                        const uint                  aTarget );

                // -----------------------------------------------------------------------------

                /**
                 * copies a pattern to another slot
                 */
                void copy_pattern( const uint & aSource,
                        const uint & aTarget );

                // -----------------------------------------------------------------------------

                /**
                 * resets the activation pattern
                 */
                void reset_pattern( const uint & aPattern );

                // -----------------------------------------------------------------------------

                /**
                 * DELETE ME
                 */
                void print_all_elements()
                {
                    Cell< Background_Element_Base* > tAllElements;
                    this->collect_all_elements( tAllElements );
                    luint tCount0 = 0;
                    std::cout << "Active Pattern: " << mActivePattern << std::endl;
                    for( auto tElement: tAllElements )
                    {
                        if ( ! tElement->is_padding() )
                        {
                            std::cout << "#Element " << tCount0++ << " " << tElement->get_hmr_id() <<
                                    " " << tElement->is_active( mActivePattern ) << std::endl;
                        }
                    }
                    std::cout << std::endl;
                }

                //--------------------------------------------------------------------------------

                /**
                 * counts elements on level aLevel and writes them into a cell
                 *
                 * @param[in]    aLevel        level to be considered
                 * @param[out]   aElementList  cell in which pointers are written
                 *
                 * @return       void
                 */
                void collect_elements_on_level( const uint                             & aLevel,
                        Cell< Background_Element_Base* > & aElementList );

                //--------------------------------------------------------------------------------

                /**
                 * Loops over all elements in frame and counts number of elements on
                 * level
                 *
                 * @param[in]    aLevel        level to be considered
                 *
                 * @return       luint         number of elements on this level
                 */
                luint count_elements_on_level( const uint& aLevel );

                // -----------------------------------------------------------------------------

                /**
                 * updates the database according to selected pattern
                 */
                void update_database();

                //------------------------------------------------------------------------------

                /**
                 * creates the faces of the background elements ( for 2D )
                 */
                void create_facets();

                //------------------------------------------------------------------------------

                /**
                 * creates the faces and edges of the background elements ( for 3D )
                 */
                void create_faces_and_edges();

                // -----------------------------------------------------------------------------

                /**
                 * delete facets
                 */
                void delete_faces();

                // -----------------------------------------------------------------------------

                /**
                 * delete edges
                 */
                void delete_edges();
                // -----------------------------------------------------------------------------

                /**
                 * reset neigbors
                 */
                void reset_neigbors();

                //------------------------------------------------------------------------------

                /**
                 * Returns the maximum number of elements on the whole mesh
                 * Padding elements do not count.
                 */
                moris_id get_max_element_domain_index() const
                {
                    return mMaxElementDomainIndex;
                }

                //--------------------------------------------------------------------------------

                /**
                 * creates a cell of all elements which have to be refined.
                 * Means the refinement queue flag of this elements is set to true.
                 *
                 * @return true if elements are flagged
                 *
                 * FIXME: This routine performs two if-checks for each element.
                 *        One check can be saved if the refinement step is performed
                 *        within this function.
                 */
                bool collect_refinement_queue();

                //------------------------------------------------------------------------------

                /**
                 * exposes the refinement queue
                 */
                Cell< Background_Element_Base *  > & get_refinement_queue()
                        {
                    return mRefinementQueue;
                        }

                //------------------------------------------------------------------------------

                /**
                 * Collect background elements on side set.
                 * Side set numbers see Exodus II : A Finite Element Data Model, p. 13
                 */
                void collect_side_set_elements(  const uint                               & aPattern,
                        const uint                               & aSideOrdinal,
                        Cell< Background_Element_Base *  > & aElements );
                //------------------------------------------------------------------------------

                /**
                 * returns the number of children per element
                 */
                uint get_number_of_children_per_element() const
                {
                    return mNumberOfChildrenPerElement;
                }

                //------------------------------------------------------------------------------

                /**
                 * reset the minumum refinement levels of all elements
                 */
                void reset_min_refinement_levels();

                //------------------------------------------------------------------------------

                virtual void get_element_in_bounding_box_memory_index( const uint                     & aPattern,
                        const moris::Matrix< DDRMat >  & aPoint,
                        const moris::Matrix< DDRMat >  & aBoundingBoxSize,
                        moris::Matrix< DDLUMat > & aElementMemoryIndex ) = 0;

                //------------------------------------------------------------------------------
            protected:
                //------------------------------------------------------------------------------

                /**
                 * Calculates a local ID from a given level
                 * and global ID
                 *
                 * @param[in]    aLevel   level of element to be investigated
                 * @param[in]    aID      global ID of element to be investigated
                 *
                 * @return       luint    local ID on submesh of proc
                 */
                virtual luint calc_subdomain_id_from_global_id( const uint  & aLevel,
                        const luint & aID) const = 0;

                //--------------------------------------------------------------------------------

                /**
                 *  After the coarsest mesh has been generated, this function tells
                 * the neighbor procs which elements are to be tagged as padding.
                 *
                 * @return void
                 *
                 */
                void synchronize_coarsest_aura();

                //--------------------------------------------------------------------------------

                /**
                 * Internal function needed to determine the required size of the cell
                 * mActiveElements before it is filled
                 *
                 * @return luint number of active elements on proc (without aura)
                 *
                 */
                luint count_active_elements() const ;

                //--------------------------------------------------------------------------------

                /**
                 * Internal function needed to determine the required size of the cell
                 * mCoarsestElementsIncludingAuras before it is filled
                 *
                 * @return luint number of active elements on proc, including aura
                 *
                 */
                luint count_active_elements_including_aura() const ;

                //--------------------------------------------------------------------------------

                /**
                 * Internal function needed to count all elements on proc
                 * ( active, refined, including aura)
                 *
                 * @return luint number of all elements on proc, including aura
                 *
                 */
                luint count_all_elements_including_aura() const ;

                //--------------------------------------------------------------------------------

                /**
                 * Internal function that calculates the level of a given pedigree ID
                 *
                 * @param[in]    aPedigreeID  pedigree ID of an arbitrary element
                 *
                 * @return       uint          level of this ID
                 *
                 */
                uint calc_level_of_pedigree_id( const luint & aPedigreeID );

                //--------------------------------------------------------------------------------

                /**
                 * Internal function that calculates the pedigree tree from an ID
                 * and stores it into a matrix
                 *
                 * @param[in]    aLevel         level of element
                 * @param[in]    aPedigreeID    pedigree ID of element
                 * @param[out]   aPedigreeTree  Matrix< DDUMat > containing the child path
                 *
                 * @return void
                 */
                /* void
            calc_pedigree_tree_from_pedigree_id(
                    const luint & aPedigreeID,
                    Matrix< DDUMat >   & aPedigreeTree); */

                //--------------------------------------------------------------------------------
                /**
                 *
                 * calculate child index from i-position ( 1D case )
                 *
                 * @param[in] aI   i-index of element
                 *
                 * @return uint child index (0-3) of element
                 *
                 */
                uint calc_child_index( const luint & aI );

                //--------------------------------------------------------------------------------

                /**
                 *
                 * calculate child index from i-position ( 2D case )
                 *
                 * @param[in] aI   i-index of element
                 * @param[in] aJ   j-index of element
                 *
                 * @return uint child index (0-3) of element
                 *
                 */
                uint calc_child_index( const luint & aI,
                        const luint & aJ );

                //--------------------------------------------------------------------------------

                /**
                 *
                 * calculate child index from i-position ( 3D case )
                 *
                 * @param[in] aI   i-index of element
                 * @param[in] aJ   j-index of element
                 * @param[in] aK   k-index of element
                 *
                 * @return uint child index (0-3) of element
                 *
                 */
                uint calc_child_index( const luint & aI,
                        const luint & aJ,
                        const luint & aK);

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
                virtual void collect_neighbors_on_level_zero() = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Loops over all elements on level including aura and counts
                 * level
                 *
                 * @param[in]    aLevel        level to be considered
                 *
                 * @return       luint         number of elements on this level
                 */
                luint count_elements_on_level_including_aura( const uint& aLevel );

                //--------------------------------------------------------------------------------

                /**
                 * counts the number of elements on level and updates
                 * mNumberOfElementsPerLevel and mMaxLevel
                 *
                 * @return  void
                 *
                 */
                void count_elements();

                //--------------------------------------------------------------------------------

                /**
                 * Populates mCoarsestPaddingElements.
                 * Called durung mesh initialization.
                 *
                 * @return void
                 */
                void collect_coarsest_padding_elements();

                //--------------------------------------------------------------------------------

                /**
                 * calculate global ID from level and from ij-position ( 1D case )
                 *
                 * @param[in] aLevel level of element
                 * @param[in] aI     i-index of element
                 *
                 * @return luint global ID of element
                 *
                 */
                virtual luint calc_domain_id_of_element( const uint  & aLevel,
                        const luint & aI ) const = 0;

                //--------------------------------------------------------------------------------

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
                virtual luint calc_domain_id_of_element( const uint  & aLevel,
                        const luint & aI,
                        const luint & aJ ) const = 0;

                //--------------------------------------------------------------------------------

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
                virtual luint calc_domain_id_of_element( const uint  & aLevel,
                        const luint & aI,
                        const luint & aJ,
                        const luint & aK) const = 0;

                //--------------------------------------------------------------------------------

                /**
                 * subroutine for collect_side_set that collects elements on coarsest
                 * level for a side
                 */
                virtual void collect_coarsest_elements_on_side( const uint                             & aSideOrdinal,
                        Cell< Background_Element_Base* > & aCoarsestElementsOnSide ) = 0;

                //--------------------------------------------------------------------------------
            private:
                //--------------------------------------------------------------------------------

                /**
                 * Private function that checks if any active element in the
                 * aura and the inverse aura is flagged for refinement.
                 * Flagged elements are communicated with the neighbor so that the
                 * refinement queue will be synchronized when created
                 *
                 * @return void
                 */
                void synchronize_refinement_queue();

                //--------------------------------------------------------------------------------

                /**
                 * Checks if element is a padding element. If element is padding get its neighbors.
                 * If neighbors are padding add them to refinement queue
                 *
                 * @return void
                 */
                virtual void check_queued_element_for_padding( Background_Element_Base * aElement ) = 0;

                //--------------------------------------------------------------------------------

                /**
                 * Private function needed to create consecutive element index
                 * over whole proc domain
                 */
                void synchronize_local_indices();

                //------------------------------------------------------------------------------

                /**
                 * calculates a staircase buffer
                 *
                 */
                void create_staircase_buffer();

                /**
                 * calculates a staircase buffer for element.
                 * Collects the neighbors of the parent element. Checks the IJK distance of the parent element and these neighbors.
                 * If the distance is smaller than the aHalfBuffer size, put this neighbor element on the refinement queue.
                 * I this neighbor element is put on refinement queue, apply this function to the neighbor element.
                 *
                 * @param[in] aElement            A background element
                 * @param[in] aElementCounter     A Counter
                 * @param[in] aHalfBuffer         The buffer size. We input only haf the buffer size, because we only have to check these neighbors
                 *
                 */
                void create_staircase_buffer_for_element(       Background_Element_Base * aElement,
                        luint                   & aElementCounter,
                        const uint                    & aHalfBuffer );
                //------------------------------------------------------------------------------
        };
    } /* namespace hmr */
} /* namespace moris */
#endif /* SRC_HMR_CL_HMR_BACKGROUND_MESH_BASE_HPP_ */
