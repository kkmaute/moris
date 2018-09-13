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
#include "cl_Mat.hpp" //LNA/src
#include "cl_Database.hpp" //MTK/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_MTK.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "HMR_Tools.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------

        // forward declaration of B-Spline mesh
        class BSpline_Mesh_Base;

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
            //! fixme: avoid copying data from field object
            Cell< Mat< real > > mFieldData;

            //! Cell containing nodal field Labels
            Cell< std::string > mFieldLabels;

            luint mNumberOfUsedAndOwnedNodes = 0;
            luint mNumberOfUsedNodes = 0;
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
            Mat< real > &
            create_field_data( const std::string & aLabel );

// ----------------------------------------------------------------------------

            /**
             * Returns a pointer to the field Data Array. Needed for MTK output.
             */
            Cell< Mat< real > > *
            get_field_data()
            {
                return & mFieldData;
            }

// ----------------------------------------------------------------------------

            Mat< real > &
            get_field_data( const uint & aFieldIndex )
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
            set_field_data( const uint& aIndex, const Mat< real > & aData )
            {
                MORIS_ERROR( aIndex < mFieldData.size(),
                             "Field does not exist" );
                mFieldData( aIndex ) = aData;
            } */

            void
            reset_fields();

            void
            add_field( const std::string & aLabel,
                       const Mat< real > & aData );

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
            MTK *
            create_mtk_object();

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

        };

// ----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_MESH_BASE_HPP_ */
