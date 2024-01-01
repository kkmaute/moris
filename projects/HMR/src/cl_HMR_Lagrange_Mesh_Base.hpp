/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Mesh_Base.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_

#include <string>

#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Basis.hpp"
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Edge.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Facet.hpp"
#include "cl_HMR_Facet_Cluster.hpp"
#include "cl_HMR_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Side_Set.hpp"
#include "cl_HMR_STK.hpp" //HMR/src
#include "moris_typedefs.hpp" //COR/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Side_Sets_Info.hpp"

#include "cl_Matrix.hpp" //LINALG/src

namespace moris::hmr
{
    // forward declaration of B-Spline mesh
    class BSpline_Mesh_Base;

    /**
     * \brief   Base class for Lagrange_Mesh
     *
     */
    class Lagrange_Mesh_Base : public Mesh_Base
    {
        Vector< bool >mBSplineMeshIsTrivialInterpoaltion;

        // @fixme: confirm that this is not identical to mAllNodesOnProc
        //! Cell containing used Nodes
        Vector< Basis * > mNodes;

        //! Cells containing real field data

        Vector< std::string >      mRealScalarFieldLabels;
        Vector< mtk::EntityRank >  mRealScalarFieldRanks;
        Vector< Matrix< DDRMat > > mRealScalarFieldData;
        Vector< Matrix< DDRMat > > mRealScalarFieldBSplineCoeffs;
        Vector< uint >             mRealScalarFieldBSplineOrders;

        luint mNumberOfUsedAndOwnedNodes = 0;
        luint mNumberOfUsedNodes         = 0;

    protected:

        uint mNumBSplineMeshes = 0;
        Vector< BSpline_Mesh_Base * > mBSplineMeshes;
        Vector< Lagrange_Mesh_Base * > mLagrangeMeshForTMatrix;

        //! IDs for MTK
        moris_id mMaxFacetDomainIndex = 0;
        moris_id mMaxEdgeDomainIndex = 0;
        moris_id mMaxNodeDomainIndex = 0;

    public:

        //! Cell containing facets
        Vector< Facet * > mFacets;

    private:

        //! Cell containing edges. Only populated in 3D
        Vector< Edge * >  mEdges;

        //! pointer to sidesets on database object
        Vector< Side_Set > * mSideSets = nullptr;

    public:

        /**
         * Default Mesh constructor
         *
         * @param[in] aParameters         container of user defined settings
         * @param[in] aBackgroundMesh   pointer to background mesh
         * @param[in] aBackgroundMesh   pointer to B-Spline mesh
         * @param[in] aOrder            polynomial degree of mesh
         */
        Lagrange_Mesh_Base ( const Parameters * aParameters,
                Background_Mesh_Base          * aBackgroundMesh,
                Vector< BSpline_Mesh_Base *  >  & aBSplineMeshes,
                uint aOrder,
                uint aActivationPattern );

        // ----------------------------------------------------------------------------

        /**
         * Virtual destructor. Does nothing.
         */
        virtual ~Lagrange_Mesh_Base(){};

        // ----------------------------------------------------------------------------

        /**
         * This function is called by the constructor, but can also be called
         * after the Lagrange mesh is generated, and the background mesh is
         * refined.
         */
        void update_mesh();

        // ----------------------------------------------------------------------------

        /**
         * called by field constructor
         */
        uint create_real_scalar_field_data(
                const std::string& aLabel,
                mtk::EntityRank    aEntityRank = mtk::EntityRank::NODE );

        /**
         * Gets the number of bases per element on this mesh
         *
         * @return Number of bases per element
         */
        uint get_number_of_bases_per_element() const
        {
            return mNumberOfBasesPerElement;
        }

        // ----------------------------------------------------------------------------

        /**
         * Returns a pointer to the field Data Array. Needed for MTK output.
         */
        Vector< Matrix< DDRMat > > & get_real_scalar_field_data()
        {
            return mRealScalarFieldData;
        }

        // ----------------------------------------------------------------------------

        Matrix< DDRMat > & get_real_scalar_field_data( uint aFieldIndex )
        {
            return mRealScalarFieldData( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        const Matrix< DDRMat > & get_real_scalar_field_data( uint aFieldIndex ) const
        {
            return mRealScalarFieldData( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        mtk::EntityRank get_real_scalar_field_rank( uint aFieldIndex ) const
        {
            return mRealScalarFieldRanks( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        Matrix< DDRMat > & get_real_scalar_field_coeffs( uint aFieldIndex )
        {
            return mRealScalarFieldBSplineCoeffs( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        const Matrix< DDRMat > & get_real_scalar_field_coeffs( uint aFieldIndex ) const
        {
            return mRealScalarFieldBSplineCoeffs( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        const std::string & get_real_scalar_field_label( uint aFieldIndex  ) const
        {
            return mRealScalarFieldLabels( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        void set_real_scalar_field_label( uint aFieldIndex,
                const std::string & aLabel )
        {
            mRealScalarFieldLabels( aFieldIndex ) = aLabel;
        }

        // ----------------------------------------------------------------------------

        uint get_real_scalar_field_bspline_order( uint aFieldIndex ) const
        {
            return mRealScalarFieldBSplineOrders( aFieldIndex );
        }

        // ----------------------------------------------------------------------------

        void set_real_scalar_field_bspline_order( uint aFieldIndex,
                uint aOrder )
        {
            mRealScalarFieldBSplineOrders( aFieldIndex ) = aOrder;
        }

        // ----------------------------------------------------------------------------

        void reset_fields();

        // ----------------------------------------------------------------------------

        /**
         * returns the number of fields
         */
        uint get_number_of_real_scalar_fields() const
        {
            return mRealScalarFieldData.size();
        }

        // ----------------------------------------------------------------------------

        /**
         * internal test that makes sure that each node is generated once
         */
        bool test_for_double_nodes();

        // ----------------------------------------------------------------------------

        /**
         * returns the number of nodes owned and shared on current proc
         */
        auto get_number_of_nodes_on_proc()
        const -> decltype ( mNumberOfUsedNodes )
        {
            return mNumberOfUsedNodes;
        }

        // ----------------------------------------------------------------------------

        //            uint get_number_of_nodes_on_proc_without_aura()
        //            {
        //                uint tNumberElementsWithoutAura = this->get_number_of_elements();
        //
        //                uint tNumBasisPerElement = this->get_number_of_bases_per_element();
        //
        //                Matrix< DDSMat > tBasisExistMat( )
        //
        //                for( uint Ik = 0; Ik < tNumberElementsWithoutAura; Ik++ )
        //                {
        //                    Element * tElement = this->get_element( Ik );
        //
        //                    for( uint Ii = 0; Ii < tNumberElementsWithoutAura; Ii++ )
        //                    {
        //                        moris_index tBasisIndex = tElement->get_basis( Ii )->get_index();
        //                    }
        //                }
        //                return mNumberOfUsedNodes;
        //            }

        // ----------------------------------------------------------------------------

        /**
         * returns a node pointer
         */
        Basis * get_node_by_index( uint aIndex )
        {
            MORIS_ASSERT( aIndex < mNodes.size(), "Requested node %-5i does not exist", aIndex );
            return mNodes( aIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * returns a node pointer ( const version )
         */
        const Basis * get_node_by_index( uint aIndex ) const
        {
            MORIS_ASSERT( aIndex < mNodes.size(), "Requested node %-5i does not exist", aIndex );
            return mNodes( aIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * returns a node pointer
         */
        Basis * get_node_by_index_including_aura( uint aIndex )
        {
            MORIS_ASSERT( aIndex < mAllBasisOnProc.size(), "Requested node %-5i does not exist", aIndex );
            return mAllBasisOnProc( aIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * returns a node pointer ( const version )
         */
        const Basis * get_node_by_index_including_aura( uint aIndex ) const
        {
            MORIS_ASSERT( aIndex < mAllBasisOnProc.size(), "Requested node %-5i does not exist", aIndex );
            return mAllBasisOnProc( aIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * returns the refinement pattern index of the B-Spline mesh
         */
        auto get_bspline_pattern( const uint aMeshIndex ) const
        -> decltype ( mBSplineMeshes( aMeshIndex )->get_activation_pattern() )
        {
            return  mBSplineMeshes( aMeshIndex )->get_activation_pattern() ;
        }

        // ----------------------------------------------------------------------------

        /**
         * Dumps the subdomain of the proc to a file.
         * By default, mesh is written as exodus.
         * Can also write vtk and gmsh files.
         * @param[in]  aFilePath  Path to where the mesh is written to.
         */
        void save_to_file( const std::string & aFilePath );

        // ----------------------------------------------------------------------------

        /**
         * creates a VTKfile of the mesh on the proc for debugging
         */
        void save_to_vtk( const std::string & aFilePath );

        // ----------------------------------------------------------------------------

        /**
         * creates a VTKfile of the mesh on the proc for debugging
         */
        void save_to_gmsh( const std::string & aFilePath );

        // ----------------------------------------------------------------------------

        /**
         * Creates the MTK output object
         */
        STK * create_stk_object( const double aTimeStep=0.0 );

        // ----------------------------------------------------------------------------

        /**
         * links the elements on the mesh with their twins on the
         * B-Spline mesh. This is needed for the T-Matrix etc
         */
        void link_twins();

        // ----------------------------------------------------------------------------

        /**
         * calculates system wide unique node indices for MTK
         *
         * @return void
         */
        void calculate_node_indices();

        // ----------------------------------------------------------------------------

        /**
         * when creating a 3D mesh with the number_aura_flag = true,
         * some few node Ids at a junction of 4 procs and with a rare special refinement around this junction
         * some of the node Ids are not communicated to the diagonal processor.
         * This function will communicate these few missed node ids.
         * This node Ids are node Ids if a a hanging node on the edge of the aura which are missing in the actual aura
         * because the aura element misses the refinement.
         * The test [Lagrange_Mesh_4_proc_problem] is producing exactly this scenario and shall serve as an example.
         *
         * @return void
         */
        void communicate_missed_node_indices();

        // ----------------------------------------------------------------------------

        /**!
         * Calculates symmetric node sharing information. Symmetric in this context
         * means that all processors know all the processors a node is shared with.
         *
         * @return void
         */
        void calculate_node_sharing();

        // ----------------------------------------------------------------------------

        /**
         * returns the number of active basis for the linked B-Spline mesh
         */
        luint get_number_of_bsplines_on_proc( uint aMeshIndex ) const                       //FIXME
        {
            // check that the requested BSpline-Mesh exists
            if ( mBSplineMeshes( aMeshIndex ) == nullptr )
            {
                return 0;
            }
            else
            {
                // return mBSplineMeshes( aMeshIndex )->get_number_of_active_basis_on_proc();
                return mBSplineMeshes( aMeshIndex )->get_number_of_indexed_basis();
            }
        }

        // ----------------------------------------------------------------------------

        Basis * get_bspline( const uint aMeshIndex,
                uint aBasisIndex )                      //FIXME
        {
            // check that the requested BSpline-Mesh exists
            MORIS_ASSERT( mBSplineMeshes( aMeshIndex ) != nullptr,
                    "Lagrange_Mesh_Base::get_bspline() - no B-Spline mesh paired with requested mesh index" );

            //                return mBSplineMeshes( aMeshIndex )->get_active_basis( aBasisIndex );
            return mBSplineMeshes( aMeshIndex )->get_basis_by_index( aBasisIndex );
        }

        // ----------------------------------------------------------------------------

        void save_faces_to_vtk( const std::string & aPath );

        // ----------------------------------------------------------------------------

        void save_edges_to_vtk( const std::string & aPath );

        // ----------------------------------------------------------------------------

        uint get_number_of_edges() const
        {
            return mEdges.size();
        }

        // ----------------------------------------------------------------------------

        uint get_number_of_facets() const
        {
            return mFacets.size();
        }

        // ----------------------------------------------------------------------------

        Facet * get_facet( uint aIndex )
        {
            return mFacets( aIndex );
        }

        // ----------------------------------------------------------------------------

        const Facet * get_facet( uint aIndex ) const
        {
            return mFacets( aIndex );
        }

        // ----------------------------------------------------------------------------

        Edge * get_edge( uint aIndex )
        {
            return mEdges( aIndex );
        }

        // ----------------------------------------------------------------------------

        const Edge * get_edge( uint aIndex ) const
        {
            return mEdges( aIndex );
        }

        // ----------------------------------------------------------------------------

        void create_facets();

        // ----------------------------------------------------------------------------

        void create_facet_clusters();

        // ----------------------------------------------------------------------------

        void create_edges();

        // ----------------------------------------------------------------------------

        moris_id get_max_element_id() const
        {
            // plus 1 is not needed, since this is actually
            // the number of entities
            return mBackgroundMesh->get_max_element_domain_index();
        }

        // ----------------------------------------------------------------------------

        moris_id get_max_facet_id() const
        {
            // plus 1 is not needed, since this is actually
            // the number of entities
            return mMaxFacetDomainIndex;
        }

        // ----------------------------------------------------------------------------

        moris_id get_max_edge_id() const
        {
            // plus 1 is not needed, since this is actually
            // the number of entities
            return mMaxEdgeDomainIndex;
        }

        // ----------------------------------------------------------------------------

        moris_id get_max_node_id() const
        {
            // plus 1 is not needed, since this is actually
            // the number of entities
            return mMaxNodeDomainIndex;
        }

        // ----------------------------------------------------------------------------

        /**
         * return the number of B-Spline meshes
         */
        uint get_number_of_bspline_meshes() const
        {
            return mNumBSplineMeshes;
        }

        // ----------------------------------------------------------------------------

        void check_if_bspline_mesh_is_trivial_interpolation()
        {
            mBSplineMeshIsTrivialInterpoaltion.resize( mNumBSplineMeshes, false );

            for( uint Ik = 0; Ik < mNumBSplineMeshes; Ik++ )
            {
                if( mBSplineMeshes( Ik ) == nullptr )
                {
                    mBSplineMeshIsTrivialInterpoaltion( Ik ) = true;
                }
            }
        }

        // ----------------------------------------------------------------------------

        bool get_bspline_mesh_is_trivial_interpolation( const uint aMeshIndex )
        {
            return mBSplineMeshIsTrivialInterpoaltion( aMeshIndex );
        }
        // ----------------------------------------------------------------------------
        /**
         * return the order of the underlying bspline mesh
         */
        uint get_bspline_order( const uint aMeshIndex )
        {
            return mBSplineMeshes( aMeshIndex )->get_min_order();
        }

        // ----------------------------------------------------------------------------
        /**
         * return the underlying bspline mesh
         */
        BSpline_Mesh_Base * get_bspline_mesh( const uint aMeshIndex )
        {
            return mBSplineMeshes( aMeshIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * return the underlying bspline mesh ( const version )
         */
        const BSpline_Mesh_Base * get_bspline_mesh( const uint aMeshIndex ) const
        {
            return mBSplineMeshes( aMeshIndex );
        }

        // ----------------------------------------------------------------------------

        /**
         * return the lagrange elements in the support of this basis
         */
        void get_my_elements_in_basis_support( const uint                    aMeshIndex,
                const uint                    aBasisIndex,
                Matrix< IndexMat > & aElementIndices )
        {
            Basis * tBasis = nullptr;

            if( mBSplineMeshes( aMeshIndex )!= nullptr)
            {
                tBasis = mBSplineMeshes( aMeshIndex )->get_basis_by_index( aBasisIndex );
            }
            else
            {
                tBasis = this->get_node_by_index_including_aura( aBasisIndex );
            }

            uint tNumElements = tBasis->get_num_elements();

            uint tCount = 0;

            for( uint Ik = 0; Ik < tNumElements; Ik ++ )
            {
                Background_Element_Base * tBackgroundElement = tBasis->get_element( Ik )->get_background_element();

                this->get_num_active_elements( tBackgroundElement, tCount );
            }

            Vector< Background_Element_Base * > tBackgroundActiveElements( tCount, nullptr );

            aElementIndices.set_size(  1, tCount, MORIS_SINT_MAX );

            tCount = 0;

            for( uint Ik = 0; Ik < tNumElements; Ik ++ )
            {
                Background_Element_Base * tBackgroundElement = tBasis->get_element( Ik )->get_background_element();

                this->get_active_elements( tBackgroundElement, tBackgroundActiveElements, tCount );
            }

            for( uint Ik = 0; Ik < tBackgroundActiveElements.size(); Ik ++ )
            {
                Element * tElement = this->get_element_by_memory_index( tBackgroundActiveElements( Ik )->get_memory_index() );

                aElementIndices( Ik ) = tElement->get_index();
            }

            MORIS_ASSERT( aElementIndices.max() != MORIS_SINT_MAX, "get_my_elements_in_basis_support(); Some invalid indices in list");

        }

        // ----------------------------------------------------------------------------

        /**
         * calculates nodes in bounding box
         */
        void calculate_nodes_indices_in_bounding_box( const moris::Matrix< DDRMat >   & aPoint,
                const moris::Matrix< DDRMat >   & aBoundingBoxSize,
                moris::Matrix< IndexMat > & aNodeIndices )
        {
            moris::Matrix< DDLUMat > tElementMemoryIndex;

            mBackgroundMesh->get_element_in_bounding_box_memory_index( mActivationPattern,
                    aPoint,
                    aBoundingBoxSize,
                    tElementMemoryIndex );

            aNodeIndices.set_size( tElementMemoryIndex.numel() * mNumberOfBasesPerElement, 1 );

            for( uint Ik = 0; Ik < tElementMemoryIndex.numel(); Ik ++ )
            {
                Element * tElement = this->get_element_by_memory_index( tElementMemoryIndex( Ik ) );

                // loop over all nodes connected to element
                for( uint k = 0; k < mNumberOfBasesPerElement; ++k  )
                {
                    // unflag basis
                    tElement->get_basis( k )->unflag();
                }
            }

            uint tCounter = 0;
            for( uint Ik = 0; Ik < tElementMemoryIndex.numel(); Ik ++ )
            {
                Element * tElement = this->get_element_by_memory_index( tElementMemoryIndex( Ik ) );

                // loop over all nodes connected to element
                for( uint k = 0; k < mNumberOfBasesPerElement; ++k  )
                {
                    // get node pointer
                    Basis * tNode = tElement->get_basis( k );

                    if( !tNode->is_flagged() )
                    {
                        aNodeIndices( tCounter++, 0 ) = tNode->get_index();

                        tNode->flag();
                    }
                }
            }

            aNodeIndices.resize( tCounter, 1 );
        }

        // ----------------------------------------------------------------------------

        /**
         * @brief Get the number of active background elements on b-spline mesh with a given mesh index
         *
         * @param aDiscretizationMeshIndex b-spline background mesh index
         * @return luint number of active elements on that b-spline mesh
         */
        luint get_num_active_bg_elements_on_discretization_mesh_index_including_aura( moris_index const aDiscretizationMeshIndex )
        {
            return mBSplineMeshes( aDiscretizationMeshIndex )->get_background_mesh()->get_number_of_active_elements_on_proc_including_aura();
        }

        // ----------------------------------------------------------------------------

        /**
         * @brief Get a list of indices of active background elements on b-spline mesh with a given mesh index
         *
         * @param aDiscretizationMeshIndex b-spline background mesh index
         * @param aElementIDs Matrix< DDLUMat > to fill with list of active BG element indices
         */
        void get_active_bg_element_indices_on_discretization_mesh_index_including_aura(
                moris_index const aDiscretizationMeshIndex,
                Matrix< DDLUMat > & aElementIDs )
        {
            mBSplineMeshes( aDiscretizationMeshIndex )->get_background_mesh()->get_active_elements_on_proc_including_aura( aElementIDs );
        }

        // ----------------------------------------------------------------------------

        void get_num_active_elements( Background_Element_Base * aBackgroundElement,
                uint                    & aCounter )
        {
            MORIS_ASSERT( aBackgroundElement != nullptr, "get_num_active_elements(); Background element is nullptr");

            bool tIsPadding = aBackgroundElement->is_padding();

            if( !tIsPadding )
            {
                bool tIsActive = aBackgroundElement->is_active( mActivationPattern );

                if( !tIsActive )
                {
                    uint tNumChildren = aBackgroundElement->get_num_children();

                    for( uint Ii = 0; Ii < tNumChildren; Ii ++ )
                    {
                        Background_Element_Base * tBackgroundElement = aBackgroundElement->get_child( Ii );

                        this->get_num_active_elements( tBackgroundElement, aCounter );
                    }
                }
                else
                {
                    aCounter++;
                }
            }
        }

        // ----------------------------------------------------------------------------

        void get_active_elements( Background_Element_Base           * aBackgroundElement,
                Vector< Background_Element_Base * > & aActiveBackgroundElements,
                uint                              & aCounter )
        {
            bool tIsPadding = aBackgroundElement->is_padding();

            if( !tIsPadding )
            {
                bool tIsActive = aBackgroundElement->is_active( mActivationPattern );

                if( !tIsActive )
                {
                    uint tNumChildren = aBackgroundElement->get_num_children();

                    for( uint Ii = 0; Ii < tNumChildren; Ii ++ )
                    {
                        Background_Element_Base * tBackgroundElement = aBackgroundElement->get_child( Ii );

                        this->get_active_elements( tBackgroundElement, aActiveBackgroundElements, aCounter );
                    }
                }
                else
                {
                    aActiveBackgroundElements( aCounter++) = aBackgroundElement;
                }
            }
        }

        // ----------------------------------------------------------------------------
        /**
         * dumps the coefficients into a binary file.
         * Format
         * < number of nodes >
         *
         * for each node
         * < node index >
         * < node id >
         * < number of bsplines >
         * < IDs of bsplines >
         * < interpolation weights >
         *
         */
        void save_coeffs_to_binary_file( const uint aOrder, const std::string & aFilePath );

        // -----------------------------------------------------------------------------

        virtual void calculate_t_matrices( bool aBool = true) = 0;

        // ----------------------------------------------------------------------------

        void set_side_sets(  Vector< Side_Set > & aSideSets )
        {
            mSideSets = & aSideSets;
        }

        // ----------------------------------------------------------------------------

        mtk::MtkSideSetInfo & get_side_set_info( const uint aIndex )
        {
            Vector< Side_Set > & tSets = *mSideSets;

            // set pointer of output object
            tSets( aIndex ).mInfo.mElemIdsAndSideOrds = & tSets( aIndex ).mElemIdsAndSideOrds;

            return tSets( aIndex ).mInfo;
        }

        // ----------------------------------------------------------------------------

        uint get_number_of_side_sets() const
        {
            if( mSideSets != nullptr )
            {
                return mSideSets->size();
            }
            else
            {
                return 0;
            }
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the lagrange elements within one bspline element
         *
         * @param aBspElementIndex index of the b-spline element
         * @param aDiscretizationMeshIndex discretization mesh index
         * @param aCells output: list of Lagrange elements as mtk::cells that sit inside the B-spline element
         */
        void get_elements_in_bspline_element(
                moris_index const aBspElementIndex,
                moris_index const aDiscretizationMeshIndex,
                Vector< mtk::Cell * > & aCells );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the lagrange elements inside the bspline elements for the whole mesh including the aura
         *
         * @param aDiscretizationMeshIndex discretization mesh index
         * @param aCells list of lists of Lagrange elements (mtk::cells) inside each B-spline element
         * @param aCellIndices list of lists of Lagrange elements (indices) inside each B-spline element
         */
        void
        get_lagrange_elements_in_bspline_elements(
                moris_index const                          aDiscretizationMeshIndex,
                Vector< Vector< mtk::Cell* > >&  aCells,
                Vector< Vector< moris_index > >& aCellIndices,
                Vector< moris_index >&                aLagToBspCellIndices,
                Vector< uint >&                       aBspCellRefineLevels,
                Vector< mtk::Cell* >&                 aBsplineCells );

        // -----------------------------------------------------------------------------

        /**
         * @brief evaluate the b-spline basis function at a given lagrange element
         *
         * @param aDiscretizationMeshIndex
         * @param aBSplineCellIndex
         * @param aLagrangeCell
         * @param aBsplineBasis
         * @param aWeights
         */
        virtual void get_extended_t_matrix(
                moris_index                                 aDiscretizationMeshIndex,
                moris_index                                 aBSplineCellIndex,
                Element&                                    aLagrangeCell,
                Vector< Vector< mtk::Vertex* > >& aBsplineBasis,
                Vector< Matrix< DDRMat > >&            aWeights ) = 0;

        // -----------------------------------------------------------------------------

        /**
         * @brief get the information to express extended basis in term of root basis
         *
         * @param[in] aDiscretizationMeshIndex b-spline mesh index
         * @param[in] aRootBSplineCellIndex  the root b-sp
         * @param[in] aExtendedBSplineCellIndex
         * @param[out] aRootBsplineBasis
         * @param[out] aExtendedBsplineBasis
         * @param[out] aWeights
         */
        virtual void get_L2_projection_matrix(
                moris_index                                 aDiscretizationMeshIndex,
                const Element*                              aRootBSplineCell,
                const Element*                              aExtendedBSplineCell,
                Vector< Vector< const mtk::Vertex* > >& aRootBsplineBasis,
                Vector< const mtk::Vertex* >&                aExtendedBsplineBasis,
                Vector< Matrix< DDRMat > >&            aWeights ) = 0;

        // -----------------------------------------------------------------------------

        /**
         * collect Lagrange elements on an BSpline interpolation element
         */
        void
        get_elements_in_interpolation_cluster(
                moris_index const aElementIndex,
                moris_index const aDiscretizationMeshIndex,
                Vector< mtk::Cell * > & aCells);

        // ----------------------------------------------------------------------------

        /**
         * @brief collect Lagrange elements on an BSpline interpolation element and side ordinal
         *
         * @param aBsplineElementIndex
         * @param aDiscretizationMeshIndex
         * @param aSideOrdinal
         * @param aCells
         */
        void
        get_elements_in_bspline_element_and_side_ordinal(
                moris_index const            aBsplineElementIndex,
                moris_index const            aDiscretizationMeshIndex,
                moris_index const            aSideOrdinal,
                Vector< mtk::Cell * > & aCells );

        /**
          * collect Lagrange elements on an BSpline interpolation element and side ordinal
          */
        void
        get_elements_in_interpolation_cluster_and_side_ordinal(
                moris_index const          aElementIndex,
                moris_index const          aDiscretizationMeshIndex,
                moris_index const          aSideOrdinal,
                Vector< mtk::Cell* >& aCells );

        // ----------------------------------------------------------------------------

        // Hack for femdoc. Only tested in serial and linear meshes.
        void nodes_renumbering_hack_for_femdoc();

    protected:

        /**
         * deletes the Lagrange meshes for T-matrices
         */

        void delete_t_matrix_lagrange_mesh();

        /**
         * returns a pointer to the child of an element if it exists
         *
         * @param[in]  aElement     pointer to Lagrange element
         * @param[in]  aChildIndex  index of requested child
         *
         * @return Element *
         */
        Element * get_child(       Element * aElement,
                uint aChildIndex );

        // ----------------------------------------------------------------------------

        /**
         * creates a facet pointer
         */
        virtual Facet * create_facet( Background_Facet * aFacet );

        // ----------------------------------------------------------------------------

        /**
         * creates an edge pointer
         */
        virtual Edge * create_edge( Background_Edge * aEdge );

        // ----------------------------------------------------------------------------

        void delete_facets();

        // ----------------------------------------------------------------------------

        void delete_edges();

        // ----------------------------------------------------------------------------
    private:
        // ----------------------------------------------------------------------------

        /**
         * creates all nodes on coarsest level
         *
         * @return void
         */
        void create_nodes_on_level_zero();

        // ----------------------------------------------------------------------------

        /**
         * creates all nodes on higher levels
         *
         * @return void
         */
        void create_nodes_on_higher_levels();

        // ----------------------------------------------------------------------------

        /**
         * creates Lagrange nodes for this mesh
         *
         * @return void
         */
        void create_nodes();

        // ----------------------------------------------------------------------------

        /**
         * Delete nodes which are only connected to padding elements.
         * We don't need them.
         *
         * @return void
         */
        void delete_nodes_on_padding_elements();

        // ----------------------------------------------------------------------------

        /**
         * loops over all elements and stores nodes in
         * mAllBasisOnProc
         */
        void collect_nodes();

        // ----------------------------------------------------------------------------

        /**
         * calculates system wide node IDs
         *
         * @return void
         */
        void calculate_node_ids();

        // ----------------------------------------------------------------------------

        /**
         * Updates the field mNodes. Called by update
         */
        void update_node_list();

        // ----------------------------------------------------------------------------

        /**
         * calculates XZY coordinates for each node
         *
         * @return void
         */
        virtual void calculate_node_coordinates() = 0;

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique node ID (1D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         *
         * @return uint          domain wide unique ID
         */
        virtual luint calculate_node_id( uint aLevel,
                luint aI ) = 0;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculates domain wide unique node ID (2D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @return uint          domain wide unique ID
         */
        virtual luint calculate_node_id( uint aLevel,
                luint aI,
                luint aJ ) = 0;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * calculates domain wide unique node ID (3D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @param[in]  aK        proc local k-position of node
         * @return uint          domain wide unique ID
         */
        virtual luint calculate_node_id( uint aLevel,
                luint aI,
                luint aJ,
                luint aK ) = 0;

        // ----------------------------------------------------------------------------

        void synchronize_facet_ids( uint aOwnedCount );

        // ----------------------------------------------------------------------------

        void negotiate_edge_ownership();

        // ----------------------------------------------------------------------------

        void synchronize_edge_ids( uint aOwnedCount );
    };
}

#endif /* SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_ */
