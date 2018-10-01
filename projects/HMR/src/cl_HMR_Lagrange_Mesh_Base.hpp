/*
 * cl_HMR_Mesh.hpp
 *
 *  Created on: May 15, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_

#include <string>

#include "typedefs.hpp" //COR/src
#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Basis.hpp"
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_STK.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Facet.hpp"
#include "cl_HMR_Edge.hpp"
namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------

        // forward declaration of B-Spline mesh
        class BSpline_Mesh_Base;
        //class Facet;

// ----------------------------------------------------------------------------

        /**
         * \brief   Base class for Lagange_Mesh
         *
         */
        class Lagrange_Mesh_Base : public Mesh_Base
        {

            //! pointer to B-Spline mesh
            BSpline_Mesh_Base * mBSplineMesh = nullptr;

            // @fixme: confirm that this is not identical to mAllNodesOnProc
            //! Cell containing used Nodes
            Cell< Basis * >     mNodes;

            //! B-Spline pattern this mesh refers to
            //uint  mBSplinePattern = 0;

            //! Cell containing nodal field data
            //! fixme: this has to be changed.
            // The Mesh is not supposed to store data
            Cell< Matrix< DDRMat > > mFieldData;

            //! Cell containing nodal field Labels
            Cell< std::string > mFieldLabels;

            luint mNumberOfUsedAndOwnedNodes = 0;
            luint mNumberOfUsedNodes = 0;

            //! Cell containing facets
            Cell< Facet * > mFacets;

            //! Cell containing edges. Only populated in 3D
            Cell< Edge * >  mEdges;

// ----------------------------------------------------------------------------
        protected:
// ----------------------------------------------------------------------------

            //! IDs for MTK
            moris_id mMaxFacetDomainIndex = 0;
            moris_id mMaxEdgeDomainIndex = 0;
            moris_id mMaxNodeDomainIndex = 0;

// ----------------------------------------------------------------------------
        public:
// ----------------------------------------------------------------------------

            /**
             * Default Mesh constructor
             *
             * @param[in] aParameters         container of user defined settings
             * @param[in] aBackgroundMesh   pointer to background mesh
             * @param[in] aBackgroundMesh   pointer to B-Spline mesh
             * @param[in] aOrder            polynomial degree of mesh
             */
            Lagrange_Mesh_Base (
                 const Parameters     * aParameters,
                 Background_Mesh_Base * aBackgroundMesh,
                 BSpline_Mesh_Base    * aBSplineMesh,
                 const uint           & aOrder,
                 const uint           & aActivationPattern );

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
            void
            update_mesh();

// ----------------------------------------------------------------------------

            /**
             * called by field constructor
             */
            uint
            create_field_data( const std::string & aLabel );

// ----------------------------------------------------------------------------

            /**
             * Returns a pointer to the field Data Array. Needed for MTK output.
             */
            Cell< Matrix< DDRMat > > &
            get_field_data()
            {
                return mFieldData;
            }

// ----------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_field_data( const uint & aFieldIndex )
            {
                return mFieldData( aFieldIndex );
            }

// ----------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_field_data( const uint & aFieldIndex ) const
            {
                return mFieldData( aFieldIndex );
            }

// ----------------------------------------------------------------------------

            std::string &
            get_field_label( const uint & aFieldIndex  )
            {
                return mFieldLabels( aFieldIndex );
            }

// ----------------------------------------------------------------------------

            /**
             * sets a field to given matrix. Needed by MTK output
             */
            /*void
            set_field_data( const uint& aIndex, const Matrix< DDRMat > & aData )
            {
                MORIS_ERROR( aIndex < mFieldData.size(),
                             "Field does not exist" );
                mFieldData( aIndex ) = aData;
            } */

            void
            reset_fields();

            /* void
            add_field( const std::string & aLabel,
                       const Matrix< DDRMat > & aData ); */

// ----------------------------------------------------------------------------

           /**
             * returns the number of fields
             */
            uint
            get_number_of_fields() const
            {
                return mFieldData.size();
            }

// ----------------------------------------------------------------------------

            /**
             * internal test that makes sure that each node is generated once
             */
            bool
            test_for_double_nodes();

// ----------------------------------------------------------------------------

            /**
             * returns the number of nodes owned and shared on current proc
             */
            auto
            get_number_of_nodes_on_proc() const -> decltype ( mNumberOfUsedNodes )
            {
                return mNumberOfUsedNodes;
            }

// ----------------------------------------------------------------------------

            /**
             * returns a node pointer
             */
            Basis*
            get_node_by_index( const uint & aIndex )
            {
                return mNodes( aIndex );
            }

// ----------------------------------------------------------------------------

            /**
             * returns a node pointer ( const version )
             */
            const Basis*
            get_node_by_index( const uint & aIndex ) const
            {
                return mNodes( aIndex );
            }

// ----------------------------------------------------------------------------

            /**
             * returns the refinement pattern index of the B-Spline mesh
             */
            auto
            get_bspline_pattern() const
                -> decltype ( mBSplineMesh->get_activation_pattern() )
            {
                return mBSplineMesh->get_activation_pattern() ;
            }


// ----------------------------------------------------------------------------

            /**
             * Dumps the subdomain of the proc to a file.
             * By default, mesh is written as exodus.
             * Can also write vtk and gmsh files.
             * @param[in]  aFilePath  Path to where the mesh is written to.
             */
            void
            save_to_file( const std::string & aFilePath );

// ----------------------------------------------------------------------------

            /**
             * creates a VTKfile of the mesh on the proc for debugging
             */
            void
            save_to_vtk( const std::string & aFilePath );

// ----------------------------------------------------------------------------

            /**
             * creates a VTKfile of the mesh on the proc for debugging
             */
            void
            save_to_gmsh( const std::string & aFilePath );

// ----------------------------------------------------------------------------

            /**
             * Creates the MTK output object
             */
            STK *
            create_stk_object();

// ----------------------------------------------------------------------------

            /**
             * links the elements on the mesh with their twins on the
             * B-Spline mesh. This is needed for the T-Matrix etc
             */
            void
            link_twins();

// ----------------------------------------------------------------------------

            /**
             * calculates system wide unique node indices for MTK
             *
             * @return void
             */
            void
            calculate_node_indices();

// ----------------------------------------------------------------------------

            /**
             * returns the number of active basis for the linked B-Spline mesh
             */
            luint
            get_number_of_bsplines_on_proc()
            {
                return mBSplineMesh->get_number_of_active_basis_on_proc();
            }

// ----------------------------------------------------------------------------

            Basis *
            get_bspline( const uint & aIndex )
            {
                return mBSplineMesh->get_active_basis( aIndex );
            }

// ----------------------------------------------------------------------------

            void
            save_faces_to_vtk( const std::string & aPath );

// ----------------------------------------------------------------------------

            void
            save_edges_to_vtk( const std::string & aPath );


// ----------------------------------------------------------------------------

            uint
            get_number_of_edges() const
            {
                return mEdges.size();
            }

// ----------------------------------------------------------------------------

            uint
            get_number_of_facets() const
            {
                return mFacets.size();
            }

// ----------------------------------------------------------------------------

            Facet *
            get_facet( const uint & aIndex )
            {
                return mFacets( aIndex );
            }

// ----------------------------------------------------------------------------

            const Facet *
            get_facet( const uint & aIndex ) const
            {
                return mFacets( aIndex );
            }

// ----------------------------------------------------------------------------

            Edge *
            get_edge( const uint & aIndex )
            {
                return mEdges( aIndex );
            }

// ----------------------------------------------------------------------------

            const Edge *
            get_edge( const uint & aIndex ) const
            {
                return mEdges( aIndex );
            }

// ----------------------------------------------------------------------------

            void
            create_facets();

// ----------------------------------------------------------------------------

            void
            create_edges();

// ----------------------------------------------------------------------------

            moris_id
            get_max_element_id() const
            {
                // plus 1 is not needed, since this is actually
                // the number of entities
                return mBackgroundMesh->get_max_element_domain_index();
            }

// ----------------------------------------------------------------------------

            moris_id
            get_max_facet_id() const
            {
                // plus 1 is not needed, since this is actually
                // the number of entities
                return mMaxFacetDomainIndex;
            }

// ----------------------------------------------------------------------------

            moris_id
            get_max_edge_id() const
            {
                // plus 1 is not needed, since this is actually
                // the number of entities
                return mMaxEdgeDomainIndex;
            }

// ----------------------------------------------------------------------------

            moris_id
            get_max_node_id() const
            {
                // plus 1 is not needed, since this is actually
                // the number of entities
                return mMaxNodeDomainIndex;
            }

// ----------------------------------------------------------------------------

            /**
             * return the order of the linked bspline mesh
             */
            uint
            get_bspline_order() const
            {
                return mBSplineMesh->get_order();
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
             * < interpolation weihhts >
             *
             */
            void
            save_coeffs_to_binary_file( const std::string & aFilePath );

// ----------------------------------------------------------------------------
        protected:
// ----------------------------------------------------------------------------


            /**
             * returns a pointer to the child of an element if it exists
             *
             * @param[in]  aElement     pointer to Lagrange element
             * @param[in]  aChildIndex  index of requested child
             *
             * @return Element *
             */
            Element *
            get_child( Element * aElement,
                       const uint            & aChildIndex );



// ----------------------------------------------------------------------------

            /**
             * creates a facet pointer
             */
            virtual
            Facet *
            create_facet( Background_Facet * aFacet );

// ----------------------------------------------------------------------------

            /**
             * creates an edge pointer
             */
            virtual Edge *
            create_edge( Background_Edge * aEdge );

// ----------------------------------------------------------------------------

            void
            delete_facets();

// ----------------------------------------------------------------------------

            void
            delete_edges();

// ----------------------------------------------------------------------------
        private:
// ----------------------------------------------------------------------------

            /**
             * creates all nodes on coarsest level
             *
             * @return void
             */
            void
            create_nodes_on_level_zero();

// ----------------------------------------------------------------------------

            /**
             * creates all nodes on higher levels
             *
             * @return void
             */
            void
            create_nodes_on_higher_levels();

// ----------------------------------------------------------------------------

            /**
             * creates Lagrange nodes for this mesh
             *
             * @return void
             */
            void
            create_nodes();

 // ----------------------------------------------------------------------------

            /**
             * Delete nodes which are only connected to padding elements.
             * We don't need them.
             *
             * @return void
             */
            void
            delete_nodes_on_padding_elements();

// ----------------------------------------------------------------------------

            /**
             * loops over all elements and stores nodes in
             * mAllBasisOnProc
             */
            void
            collect_nodes();

// ----------------------------------------------------------------------------

            /**
             * calculates system wide node IDs
             *
             * @return void
             */
            void
            calculate_node_ids();

// ----------------------------------------------------------------------------

            /**
             * Updates the field mNodes. Called by update
             */
            void
            update_node_list();

// ----------------------------------------------------------------------------

            /**
             * calculates XZY coordinates for each node
             *
             * @return void
             */
            virtual void
            calculate_node_coordinates() = 0;

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
            virtual luint
            calculate_node_id(
                    const uint  & aLevel,
                    const luint & aI ) = 0;

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
            virtual luint
            calculate_node_id(
                    const uint  & aLevel,
                    const luint & aI,
                    const luint & aJ ) = 0;

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
            virtual luint
            calculate_node_id(
                    const uint  & aLevel,
                    const luint & aI,
                    const luint & aJ,
                    const luint & aK ) = 0;

// ----------------------------------------------------------------------------

            void
            synchronize_facet_ids( const uint & aOwnedCount );

// ----------------------------------------------------------------------------

            void
            synchronize_edge_ids( const uint & aOwnedCount );

// ----------------------------------------------------------------------------

            //void
            //link_facet_children_2d();

// ----------------------------------------------------------------------------

            //void
            //link_facet_children_3d();

// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
        };

// ----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_ */
