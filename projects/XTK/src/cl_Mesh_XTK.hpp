/*
 * cl_Mesh_XTK.hpp
 *
 *  Created on: Nov 9, 2016
 *      Author: doble
 */

#ifndef SRC_XTK_CL_MESH_XTK_HPP_
#define SRC_XTK_CL_MESH_XTK_HPP_

#include "cl_XTK_Enums.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_Mesh_Enums.hpp" // MTK/src
#include "cl_Cell.hpp" // CON/src
#include "cl_Debug.hpp" // TOL/src
namespace moris
{
    class MeshXTK
    {

    public:

        // MeshXTK is simple mesh class for internal use within the XFEM toolkit.
        //

        /*
         * Constuctor
         * D - spatial dimension of an element/Cell
         * i.e for 3d D = 3;
         *     for 2d D = 2;
         */
        MeshXTK(moris::uint D);

        MeshXTK(){}

        ~MeshXTK();

        /*
         * set_template_type stores the template type used to create mesh
         *
         * @param[in] aTemplate    - Specifies the template used
         */

        void
        generate_templated_mesh(enum TemplateType    aTemplate);

        /* @ brief Sets up templated ancestry
         *
         * @ param[in] aTemplate       - specifies the template ancestry to use
         * @ param[in] aParentEntities - cell of row vectors of parent entity indices
         * Require: parents entities live in a full Mesh
         */
        void
        set_entity_ancestry(enum TemplateType                        aTemplate,
                            moris::Cell<moris::Mat<moris::uint>>   & aParentEntities);

        void
        modify_templated_mesh(enum TemplateType    aTemplate);

        /*
         * In the Auxiliary connectivity, there are local XTK indice. These can be converted to global ids or local indices if necessary
         * this allows for setting up aux connectivity prior to assigning global ids and process local indices. Assumes entire XTK mesh
         * is of the same type(i.e. TET4)
         *
         * mAuxConnectivity Structure:
         * rows             - Number of entities with dimension aD
         * col 0            - Counter
         * col 1-(n+1)      - aDprime index
         * col (n+1)-(2n+1) - aDprime index
         */
        void
        init_aux_connectivity(enum EntityRank aD,
                              enum EntityRank aDPrime1,
                              enum EntityRank aDPrime2);

        /*
         * aFlag - 0 means the provided aDPrime1Ind is appended to the end of existing nodes
         *       - 1 means the provided aDPrime1Ind is an XTK index
         *
         * aDPrime2Ind must be XTK local index
         */
        void
        add_entity_to_aux_connectivity(moris::uint aDPrime1Ind,
                                       moris::uint aDPrime2Ind,
                                       moris::uint aFlag);


        /*
         * For template generation purposes, should not be used in general.
         */
        void
        set_aux_connectivity(moris::Mat<moris::uint>  aAuxConn);

        moris::Mat<moris::uint>
        get_aux_connectivity();

        moris::Mat<moris::uint>
        sort_nodes(enum TemplateType       aTemplate,
                   moris::Mat<moris::uint> aAuxConnectivity);


        // Public Functions to build/alter Mesh ----------------------------------
        /*
         * set_nodes_index stores a set of ordered node Inds in the XTK Mesh
         */
        void
        set_node_index(moris::Mat<moris::uint>    & aNodes);

        /*
         * set_node_ids stores a set of ordered node inds in the XTK Mesh
         */
        void
        set_node_ids(moris::Mat<moris::uint>    & aNodeIds);

        /*
         * Tells the XTK mesh where a node index will be placed once it has been communicated
         * This function is used if you have all node pointers else use
         */
        void
        set_pending_node_index_pointers(moris::Cell<moris::uint*>  aNodeIndPtr);

        /*
         * XTK mesh retrieves pending node indices. Prior to this call the unique node assignments
         * must be complete via the request structure in XTK model
         */
        void
        get_pending_node_inds();

        /*
         * add_nodes adds node IDs to the existing list of nodes
         * If possible, add nodes in bulk
         */
        void
        add_node_ind(moris::Mat<moris::uint>    &aNodeInd);

        moris::uint
        get_node_ind(moris::uint  aIndex);

        moris::Mat<moris::uint>
        get_all_node_inds();
//        /*
//         * Gives the XTK mesh a reference to its nodal coordinates via a model. The XTK mesh does
//         * not own the coordinates
//         */
//         void
//         set_model_pointer(moris::xtk::Model  *  pModel);

        // Mesh Entity Accessing Public Functions --------------------------------

        moris::uint
        get_num_entities(enum EntityRank aEntity);

        moris::uint
        get_num_entities_connected_to_entity(enum EntityRank aEntity,
                                             enum EntityRank aEntityPrime,
                                             moris::uint     aEntityIndex);

        /*
         * returns the entities of rank aEntityPrime connected to entity of rank aEntity
         *
         * @param[in] aEntity      - Parent Entity rank
         * @param[in] aEntityPrime - Secondary Entity Rank
         * @param[in] aEntityIndex - Parent Entity Index
         * @param[in] aType        - Specifies if a xtk local (0) or processor local inds  (1) returned
         */
        moris::Mat<moris::uint>
        get_entities_connected_to_entity(enum EntityRank aEntity,
                                         enum EntityRank aEntityPrime,
                                         moris::uint     aEntityIndex,
                                         moris::uint     aType);

        /*
         * Returns the parent entity rank and index that an xtk entity was built on
         * @param[in] aEntityRank  - XTK entity rank
         * @param[in] aEntityIndex - XTK entity index
         *
         * Returns a row matrix
         */
        moris::Mat<moris::uint>
        get_parent_entity(enum EntityRank aEntityRank,
                          moris::uint     aEntityIndex);




        /*
         * Sets the parent element index that this mesh was created from
         */
        void
        set_parent_element_index(moris::uint    aEID);

        /*
         * returns the parent element index
         */
        moris::uint get_parent_element_index();



        /*
         * Accesses connectivity for use when asking the mesh for an entity
         * i.e. get_nodes_connected_to_element for use example.
         */

        moris::Mat<moris::uint>
        access_connectivity(enum EntityRank             aEntity,
                            enum EntityRank             aEntityPrime,
                            moris::uint                 aEntityIndex,
                            moris::uint                 aType)  const;

        moris::Mat<moris::uint>
        access_connectivity(moris::Mat<moris::uint>  &aInds,
                            moris::Mat<moris::uint>  &aOffs,
                            moris::uint                    aEntityIndex)const;

        /*
         * Returns a number of entities for intialization or loop sizes
         * i.e see get_num_nodes.
         */
        moris::uint
        get_num_entities(moris::uint    aIndex,
                         moris::uint    aEntityID) const;


        // Mesh Connectivity public functions-------------------------------------


        /*
         * set connectivity indices and offsets for a d->dprime connectivity
         * DO NOT USE WITH TEMPLATE TOPOLOGIES!
         */

        void
        set_connectivity(moris::uint                d,
                         moris::uint                dprime,
                         moris::Mat<moris::uint>    aIndices,
                         moris::Mat<moris::uint>    aOffsets);

        /*
         *Returns matrix of connectivity for entire xtk mesh (local Indexes)
         *Requires that the type of connectivity is uniform (allows for return in moris::Mat rather than a moris::cell)
         *@param[in] aType - specifies whether XTK local (0) or Processor Local Inds(1), or Global Ids(2) only valid for nodes
         */

        moris::Mat<moris::uint>
        get_full_connectivity(enum EntityRank    d,
                              enum EntityRank    dprime,
                              moris::uint        aType);

    private:
        moris::uint                    mD;                    // Spatial dimension of an element
        moris::uint                    mNumElems;             // Number of elements
        moris::uint                    mParentEID;            // Parent element Index
        moris::Mat<moris::uint>        mNodesInd;             // Ordered list of Node Indices (index is local to processor) (stored in a row)
        moris::Mat<moris::uint>        mNodesId;              // Ordered list of Node IDs (Needed for hierachal subdivision tempate selection) (id is global) (stored in a row)
        moris::Cell<moris::uint*>      mPtrPendingNodeIndex;

        // Auxilliary connectivity data structure (Might move into an actual struct or class)
        moris::Mat<moris::uint>        mAuxConnectivity;
        moris::uint                    mAuxConnectivityNum;
        enum EntityRank                mAuxD;
        enum EntityRank                mAuxDPrime1;
        enum EntityRank                mAuxDPrime2;
        // End


        // Mesh Topology - member variables --------------------------

        /* The topology data is a cell of connectivity information. Each connectivity has an indices vector and an offset vector.
         * The advantage of structuring connectivity in this format is eliminates the need for cells of cells;
         *
         *
         * 0   - Node
         * 1   - Edge
         * 2   - Face
         * 3   - Element
         * D-1 - Facet
         * D   - Cell
         * where D is the dimensionality of the element/cell
         *
         * in 2D
         * Edge = Facet
         * Face = Element
         *
         * in 3D
         * Face = facet
         * Element = Cell
         *
         * Connectivities are stored in the following format
         * while d<D
         * mIndices((mD+1) * d + dprime);
         * mOffsets((mD+1) * d + dprime);
         *
         *
         * For D = 3
         * Indices ---------------------------------------
         * OD -> dprime connectivity indices
         * mIndices(0)  = 0D to 0D connectivity
         * mIndices(1)  = 0D to 1D Connectivity
         * mIndices(2)  = 0D to 2D Connectivity
         * mIndices(3)  = 0D to 3D Connectivity
         *
         * 1D -> dprime connectivity indices
         * mIndices(4)  = 1D to 0D Connectivity
         * mIndices(5)  = 1D to 1D Connectivity
         * mIndices(6)  = 1D to 2D Connectivity
         * mIndices(7)  = 1D to 3D Connectivity
         *
         * 2D -> dprime connectivity indices
         * mIndices(8)  = 2D to 0D Connectivity
         * mIndices(9)  = 2D to 1D Connectivity
         * mIndices(10) = 2D to 2D Connectivity
         * mIndices(11) = 2D to 3D Connectivity
         *
         * 3D -> dprime connectivity indices
         * mIndices(12) = 3D to 0D Connectivity
         * mIndices(13) = 3D to 1D Connectivity
         * mIndices(14) = 3D to 2D Connectivity
         * mIndices(15) = 3D to 3D Connectivity
         *
         *
         * Offset -------------------------------------------
         *
         * OD -> dprime connectivity offsets
         * mOffsets(0)  = 0D to 0D connectivity
         * mOffsets(1)  = 0D to 1D Connectivity
         * mOffsets(2)  = 0D to 2D Connectivity
         * mOffsets(3)  = 0D to 3D Connectivity
         *
         * 1D -> dprime connectivity offsets
         * mOffsets(4)  = 1D to 0D Connectivity
         * mOffsets(5)  = 1D to 1D Connectivity
         * mOffsets(6)  = 1D to 2D Connectivity
         * mOffsets(7)  = 1D to 3D Connectivity
         *
         * 2D -> dprime connectivity offsets
         * mOffsets(8)  = 2D to 0D Connectivity
         * mOffsets(9)  = 2D to 1D Connectivity
         * mOffsets(10) = 2D to 2D Connectivity
         * mOffsets(11) = 2D to 3D Connectivity
         *
         * 3D -> dprime connectivity offsets
         * mOffsets(12) = 3D to 0D Connectivity
         * mOffsets(13) = 3D to 1D Connectivity
         * mOffsets(14) = 3D to 2D Connectivity
         * mOffsets(15) = 3D to 3D Connectivity
         *
         *
         * Ancestry -------------------------------------------
         * mAncestry(0) = 0D to 0D Ancestry
         * mAncestry(1) = 1D to 1D Ancestry
         * mAncestry(2) = 2D to 2D Ancestry
         * mAncestry(3) = 3D to 3D Ancestry
         */

        moris::Cell<moris::Mat<moris::uint>> mIndices;
        moris::Cell<moris::Mat<moris::uint>> mOffsets;
        moris::Cell<moris::Mat<moris::uint>> mAncestryInds;
        moris::Cell<moris::Mat<moris::uint>> mAncestryRank;

        /*
         * maps d and d prime to the index where this type of connectivity is stored in indices and offsets data
         *
         * @param[in] d            first connectivity dimension
         *
         * @param[in] dprime       second connectivity dimension
         *
         * @return location index
         *
         * i.e to access 3d->0d connectivity d = 3 and dprime = 0, return i = 12.
         */
        moris::uint
        index_map(moris::uint       d,
                  moris::uint       dprime) const;

        moris::Cell<moris::Mat<moris::uint>>
        access_template_indices(enum TemplateType  aTemplate);

        moris::Cell<moris::Mat<moris::uint>>
        access_template_offsets(enum TemplateType aTemplate);

        void
        replace_parent_entity(enum EntityRank  aParentRank,
                              moris::uint      aParentInd,
                              moris::Mat<moris::uint> tChildrenConn,
                              moris::Mat<moris::uint> tSortedNodes);
    };
}



#endif /* SRC_XTK_CL_MESH_XTK_HPP_ */
