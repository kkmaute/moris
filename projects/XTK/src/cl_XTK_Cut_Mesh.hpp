/*
 * cl_XTK_Cut_Mesh.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_CUT_MESH_HPP_
#define SRC_XTK_CL_XTK_CUT_MESH_HPP_

// XTKL: Linear Algebra includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"

// XTKL: Mesh includes
#include "cl_Mesh_Enums.hpp" // For entity rank
#include "cl_MTK_Mesh.hpp"

// XTKL: Xtk includes
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "cl_XTK_Downward_Inheritance.hpp"
#include "cl_XTK_Interface_Element.hpp"

using namespace moris;

namespace xtk
{
    // forward declaration
    class Model;
    class Mesh_Cleanup;

    class Cut_Mesh
    {
    public:
        // ----------------------------------------------------------------------------------

        Cut_Mesh();

        // ----------------------------------------------------------------------------------

        Cut_Mesh(Model *aModel,
                 moris::uint aModelDim);

        // ----------------------------------------------------------------------------------

        Cut_Mesh(Model *aModel,
                 moris::size_t aNumSimpleMesh,
                 moris::size_t aModelDim);

        // ----------------------------------------------------------------------------------

        ~Cut_Mesh();

        // ----------------------------------------------------------------------------------

        /**
         * Allocate space for new children meshes
         * @param aNumSimpleMesh - number of simple meshes to allocate space for
         * @param aModelDim - model dimension (1D, 2D or 3D)
         */
        void
        inititalize_new_child_meshes(
            moris::size_t aNumNewChildMesh,
            moris::size_t aModelDim);

        // ----------------------------------------------------------------------------------

        /**
         * generate_templated_mesh genereates a templated mesh on a given child mesh with the specified edge intersection permutation (if required by a template)
         *
         * @param[in] aChildMeshIndex  - specifies which simple mesh to generate the template on
         * @param[in] aTemplate - Specifies the template used
         * @param[in] aPermutation - If required, specifies the edge intersection pattern
         */
        void
        generate_templated_mesh(
            moris::size_t aChildMeshIndex,
            enum TemplateType aTemplate);

        // ----------------------------------------------------------------------------------

        /**
         * generate_templated_mesh tells all simple meshes to generate a specified templated mesh
         *
         * @param[in] aTemplate  - Specifies the template to use
         */
        void
        generate_templated_mesh(
            enum TemplateType aTemplate);

        // ----------------------------------------------------------------------------------

        /*
         * Generate a templated mesh for a subset of the children meshes.
         */
        void
        generate_templated_mesh(
            Matrix<IndexMat> const &aChildMeshIndices,
            enum TemplateType aTemplate);

        // ----------------------------------------------------------------------------------

        /*
         * Converts existing tet4 child mesh to tet10s
         *
         */
        void
        convert_cut_mesh_to_tet10s();

        // ----------------------------------------------------------------------------------

        /**
         * @ brief Sets up template ancestry with parametric information
         * @param[in] aChildMeshIndex        - simple mesh index
         * @param[in] aTemplate       - specifies the template ancestry to use
         * @param[in] aParentEntities - cell of row vectors of parent entity indices
         */
        void
        initialize_new_mesh_from_parent_element(
            moris::size_t aChildMeshIndex,
            enum TemplateType aTemplate,
            Matrix<IndexMat> &aNodeIndices,
            Cell<Matrix<IndexMat>> &aParentEntities);

        // ----------------------------------------------------------------------------------

        //Modify Template mesh is only useful for unit tests without node vars
        void
        modify_templated_mesh(
            moris::size_t aChildMeshIndex,
            enum TemplateType aTemplate);

        // ----------------------------------------------------------------------------------
        /**
         * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
         *       - 1 means the provided aDPrime1Ind is an XTK index
         *
         * aDPrime2Ind must be XTK local index
         */
        void
        add_entity_to_intersect_connectivity(
            moris::size_t aChildMeshIndex,
            moris::size_t aDPrime1Ind,
            moris::size_t aDPrime2Ind,
            moris::size_t aReturnType);

        // ----------------------------------------------------------------------------------

        /*!
         * add an interface element
         */
        void
        add_interface_element(
            Interface_Element const & aInterfaceElement);

        // ----------------------------------------------------------------------------------

        /*!
         * Get all interface elements
         */
        moris::Cell<Interface_Element> &
        get_interface_elements();

        // ----------------------------------------------------------------------------------

        /*!
         * Get the ids of all interface elements
         */
        moris::Matrix<moris::IdMat>
        get_interface_element_ids();

        // ----------------------------------------------------------------------------------

        /*!
         * Get index of extracted (prism) interface elements
         */
        moris::Matrix<moris::IndexMat>
        get_extracted_interface_elements_loc_inds();

        // ----------------------------------------------------------------------------------

        /*!
         * Set node indices in a child mesh
         * @param[in] aChildMeshIndex - Child mesh index
         * @param[in] aNodeInd - Node indices
         */
        void
        set_node_index(
            moris::size_t const &aChildMeshIndex,
            Matrix<IndexMat> &aNodeInd);

        // ----------------------------------------------------------------------------------

        /*!
         * Set node ids in a child mesh
         * @param[in] aChildMeshIndex - Child mesh index
         * @param[in] aNodeInd - Node ids
         */
        void
        set_node_ids(
            moris::size_t const &aChildMeshIndex,
            moris::Matrix<moris::IdMat> &aNodeIds);

        // ----------------------------------------------------------------------------------

        /*!
         * Set node ids in a child mesh
         * @param[in] aChildMeshIndex - Child mesh index
         * @param[in] aNodeInd - Node ids
         */
        void
        add_node_ids(
            moris::size_t const &aChildMeshIndex,
            Matrix<IdMat> &aNodeIds);

        // ----------------------------------------------------------------------------------

        /*!
         * Set element ids in a child mesh
         * @param[in] aChildMeshIndex - Child mesh index
         * @param[in] aElementIdOffset - Element Id offset
         */
        void
        set_child_element_ids(
            moris::size_t const &aChildMeshIndex,
            moris::moris_id &aElementIdOffset);

        // ----------------------------------------------------------------------------------

        /*!
         * Set element indicess in a child mesh
         * @param[in] aChildMeshIndex - Child mesh index
         * @param[in] aElementIdOffset - Element Ind offset
         */
        void
        set_child_element_inds(
            moris::size_t const &aChildMeshIndex,
            moris::moris_index &aElementIndOffset);

        // ----------------------------------------------------------------------------------

        /*!
         * Get element Ids in a child mesh
         */
        moris::Matrix<moris::IdMat> const &
        get_element_ids(
            moris::size_t const &aChildMeshIndex);

        // ----------------------------------------------------------------------------------

        /*!
         * Get all element ids in cut mesh
         */
        moris::Matrix<moris::IdMat>
        get_all_element_ids();

        // ----------------------------------------------------------------------------------

        /*!
         * Get all element indices in cut mesh
         */
        moris::Matrix<moris::IndexMat>
        get_all_element_inds();

        // ----------------------------------------------------------------------------------

        /*!
         * Get element Ids in the cut mesh of a given id
         */
        Cell<moris::Matrix<moris::IdMat>>
        get_child_elements_by_phase(
            uint aNumPhases,
            moris::mtk::Mesh const &aBackgroundMeshData);

        // ----------------------------------------------------------------------------------

        /*!
         * Get element Inds from a child mesh
         */
        moris::Matrix<moris::IndexMat> const &
        get_element_inds(
            moris::size_t const &aChildMeshIndex) const;

        // ----------------------------------------------------------------------------------

        /*!
         * Get node indices for a given child mesh
         */
        moris::Matrix<moris::IndexMat> const &
        get_node_indices(
            moris::size_t aChildMeshIndex);

        // ----------------------------------------------------------------------------------

        void
        set_child_element_topology(
            enum CellTopology aChildCellTopo);

        // ----------------------------------------------------------------------------------

        enum CellTopology
        get_child_element_topology();

        // ----------------------------------------------------------------------------------

        moris::size_t
        get_num_child_meshes() const;

        // ----------------------------------------------------------------------------------

        moris::size_t
        get_num_entities(
            moris::size_t aChildMeshIndex,
            enum EntityRank aEntityRank) const;

        // ----------------------------------------------------------------------------------

        moris::size_t
        get_num_entities(
            enum EntityRank aEntityRank) const;

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_parent_element_index(
            moris::size_t aChildMeshIndex);

        // ----------------------------------------------------------------------------------

        /*!
        * Returns all the children elements connected to a provided face index in a single child mesh
        *
        * @param[in]  aChildMeshIndex         - Child Mesh Index
        * @param[in]  aParentFaceIndex        - Parent Face Index
        * @param[out] aChildrenElementId      - Children Element Global Id
        * @param[out] aChildrenElementCMInd   - Child element index local to child mesh
        * @param[out] aFaceOrdinal            - Face Ordinal relative to element
        */
        void get_child_elements_connected_to_parent_facet(
            moris::moris_index const &aChildMeshIndex,
            moris::moris_index const &aParentFaceIndex,
            moris::Matrix<moris::IdMat> &aChildrenElementId,
            moris::Matrix<moris::IndexMat> &aChildrenElementCMInd,
            moris::Matrix<moris::IndexMat> &aFaceOrdinal) const;

        // Outputting Mesh information

        // ----------------------------------------------------------------------------------

        void pack_cut_mesh_by_phase(
            moris::size_t const &aMeshIndex,
            moris::size_t const &aNumPhases,
            Cell<moris::Matrix<moris::IdMat>> &aElementIds,
            Cell<moris::Matrix<moris::IdMat>> &aElementCMInds) const;

        // ----------------------------------------------------------------------------------

        /*!
        * pack interface side sets, if phase index not provided all are packaged
        * flag  = 0 - ids
        * flag  = 1- indices
        * flag  = other - cm indices
        */
        moris::Matrix<moris::IdMat>
        pack_interface_sides(
            moris_index aGeometryIndex,
            moris_index aPhaseIndex0,
            moris_index aPhaseIndex1,
            moris_index aIndexFlag = 0) const;

        // ----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        pack_interface_sides_loc_inds() const;

        // ----------------------------------------------------------------------------------

        /*
         * Get full element to node connectivity (ids). Full here means for all children meshes
         */
        moris::Matrix<moris::IdMat>
        get_full_element_to_node_glob_ids();

        // ----------------------------------------------------------------------------------

        /*
        * Get full element to node connectivity (indices). Full here means for all children meshes
        */
        moris::Matrix<moris::IdMat>
        get_full_element_to_node_loc_inds();

        // ----------------------------------------------------------------------------------

        /*
        * Get full element to node glob ids by phase
        */
        moris::Cell<moris::Matrix<moris::IdMat>>
        get_full_element_to_node_by_phase_glob_ids(
            moris::uint aNumPhases,
            moris::mtk::Mesh &aBackgroundMeshData);

        // ----------------------------------------------------------------------------------

        /*
        * Get a child mesh reference
        */
        Child_Mesh const &
        get_child_mesh(
            moris::size_t const &aChildMeshIndex) const;

        // ----------------------------------------------------------------------------------

        /*
        * non const version
        */
        Child_Mesh &
        get_child_mesh(
            moris::size_t const &aChildMeshIndex);

        // ----------------------------------------------------------------------------------

        void
        set_num_subphases(
            moris::uint aNumSubPhases);

        // ----------------------------------------------------------------------------------

        moris::uint
        get_num_subphases();

        // ----------------------------------------------------------------------------------

        void
        setup_subphase_to_child_mesh_connectivity();

        // ----------------------------------------------------------------------------------

        Matrix<IndexMat> const &
        get_subphase_to_child_mesh_connectivity();

        // ----------------------------------------------------------------------------------

        moris_id
        get_subphase_id(
            moris_index aSubphaseIndex);

        // ----------------------------------------------------------------------------------

        uint
        get_bulk_phase_index(
            moris_index aSubPhaseIndex);

        // ----------------------------------------------------------------------------------

        void
        populate_subphase_vector(
            moris::Matrix<moris::IndexMat> &aSubphase);
        // ----------------------------------------------------------------------------------

        moris::Memory_Map
        get_memory_usage();

        // ----------------------------------------------------------------------------------

    protected:
        /*!
         * Tell cut mesh about groupings of its children meshes
         */
        void
        add_child_mesh_groups(
            Cell<Child_Mesh *> &tOwnedChildrenMeshes,
            Cell<Child_Mesh *> &tNotOwnedChildrenMeshes,
            Cell<moris_id> &tNotOwnedOwningProc);

        // ----------------------------------------------------------------------------------
        /*!
             * Get the child meshes which are owned and not shared
             */
        Cell<Child_Mesh *> &
        get_owned_child_meshes();

        // ----------------------------------------------------------------------------------

        /*!
             * Get the child mesh which are not owned
             */
        Cell<Child_Mesh *> &
        get_not_owned_child_meshes();

        // ----------------------------------------------------------------------------------

        /*
             * Get the not owned children mesh owners
             */
        Cell<moris_id> &
        get_not_owned_child_owners();

        // ----------------------------------------------------------------------------------

        void
        remove_all_child_meshes_but_selected(Cell<moris::uint> const & aMeshesToKeep,
                                             Cell<moris::uint> const & aMeshesToDelete,
                                             Cell<moris_index>       & aCellsToRemoveFromMesh);

        // ----------------------------------------------------------------------------------

        void
        reindex_cells(Cell<moris_index> & aOldIndexToNewCellIndex);

        // ----------------------------------------------------------------------------------

        friend Model;
        friend Mesh_Cleanup;

    private:
        Model *mModel;

        // ----------------------------------------------------------------------------------

        // spatial dimension
        moris::uint mSpatialDim;

        moris::size_t mNumberOfChildrenMesh;

        // All children meshes
        Cell<Child_Mesh*> mChildrenMeshes;

        // Groupings of children meshes determined by parent cell information (pointers to the mChildrenMeshes data)
        Cell<Child_Mesh *> mOwnedChildrenMeshes;    /* All children meshes which are fully owned by this processor and not shared with another processor */
        Cell<Child_Mesh *> mNotOwnedChildrenMeshes; /* All children meshes which are shared with another processors and not owned by this processor */
        Cell<moris_id>     mNotOwnedOwningProc;         /* For the mOwnedSharedChildrenMeshes, the other processes which share this mesh */

        // Interface elements
        Cell<Interface_Element> mInterfaceElements;

        // Number of entities total in child meshes and if current count is accurate
        // mutable because some const function need this information, and if the counts
        // are not consistent we need to be able to update these vars
        mutable bool mConsistentCounts;
        mutable Cell<moris::size_t> mNumEntities;

        // topology of child elements (i.e. TET4)
        enum CellTopology mChildElementTopo;

        // number of subphases
        moris::uint mNumSubPhases;
        Matrix<IndexMat> mSubPhaseIndexToChildMesh;
        Matrix<IndexMat> mSubPhaseIndexToChildMeshSubphaseIndex;

    private:
        // ----------------------------------------------------------------------------------

        void
        get_entity_counts() const;
    };
} // namespace xtk

#endif /* SRC_XTK_CL_XTK_CUT_MESH_HPP_ */
