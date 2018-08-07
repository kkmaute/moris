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
// ----------------------------------------------------------------------------
        protected:
// ----------------------------------------------------------------------------

            //! counter for nodes this proc owns
            luint mNumberOfOwnedNodes = 0;

            //! number of nodes used by this proc
            luint mNumberOfNodes = 0;

// ----------------------------------------------------------------------------
        public:
// ----------------------------------------------------------------------------

            /**
             * Default Mesh constructor
             *
             * @param[in] aParameters         container of user defined settings
             * @param[in] aBackgroundMesh   pointer to background mesh
             * @param[in] aOrder            polynomial degree of mesh
             */
            Lagrange_Mesh_Base (
                 const Parameters       * aParameters,
                 Background_Mesh_Base * aBackgroundMesh,
                 const uint           & aOrder );

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
             * internal test that makes sure that each node is generated once
             */
            bool
            test_for_double_nodes();

// ----------------------------------------------------------------------------
            /**
             * returns the number of nodes owned and shared on current proc
             */
            auto
            get_number_of_nodes_on_proc() -> decltype ( mNumberOfNodes )
            {
                return mNumberOfNodes;
            }

// ----------------------------------------------------------------------------

            /**
             * returns the number of nodes
             * ( refers to initialization or last call of update_mesh)
             *
             * @return luint
             */
/*            auto
            get_number_of_nodes_including_aura() const
                ->decltype( mNumberOfAllBasis )
            {
                return mNumberOfAllBasis;
            } */

// ----------------------------------------------------------------------------

           /**
             * tells how many nodes belong to an element
             *
             * @return luint
             */
            /* auto
            get_number_of_nodes_per_element() const ->decltype( mNumberOfBasisPerElement )
            {
                return mNumberOfBasisPerElement ;
            } */

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
            link_twins( BSpline_Mesh_Base* aBSplineMesh );

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
             * calculates system wide unique node indices for MTK
             *
             * @return void
             */
            void
            calculate_node_indices();

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
