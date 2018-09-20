//
// Created by messe on 3/7/18.
//

#ifndef MORIS_CL_GE_SDF_MESH_DATA_HPP
#define MORIS_CL_GE_SDF_MESH_DATA_HPP

#include <map>
#include "chronos.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Map.hpp" // CON/src
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"
#include "cl_Mesh_Enums.hpp" // MTK/src
#include "cl_Database.hpp" // MTK/src
//#include "cl_Hierarchical_Mesh.hpp" // STK/src/Heirarchical
#include "GeUtilities.hpp"

namespace ge {
// -----------------------------------------------------------------------------
    /**
      * @brief The struct SDF_Mesh_Data contains mesh informations.
      *        It is updated before the raycast algorithm is called
      */
    class SDF_Mesh_Data
    {
        const moris::database &mBackgroundMesh; //!< a wrapper of the mesh we are working with
        moris::Matrix< moris::DDRMat > mNodeCoords;     //!< contains the coordinates of the nodes
        moris::Matrix< moris::DDUMat > mNumberOfNodesPerElement;
        moris::Matrix< moris::DDUMat > mElementTopology;

        moris::uint mNumberOfNodes = 0;             //!< number of nodes in the mesh on current proc
        moris::uint mNumberOfElements = 0;          //!< number of cells in the mesh on current proc
        moris::Matrix< moris::DDRMat > mMinNodeCoordinate; //!< minimum coordinates of the mesh
        moris::Matrix< moris::DDRMat > mMaxNodeCoordinate; //!< minimum coordinates of the mesh

        moris::Matrix< moris::DDUMat > mNodesOnProc;
        moris::Matrix< moris::DDUMat > mElementsOnProc;
        //moris::Matrix< moris::DDUMat > mLocalNodesOnProc;
        //moris::Matrix< moris::DDUMat > mLocalElementsOnProc;
// =============================================================================
    public:
// =============================================================================

        /**
        * @brief             The constructor for the mesh data struct
        *
        */
        SDF_Mesh_Data (const moris::database &aBackgroundMesh) :
                mBackgroundMesh(aBackgroundMesh),
                mMinNodeCoordinate(3, 1),
                mMaxNodeCoordinate(3, 1)
        {
        }

// -----------------------------------------------------------------------------

        ~SDF_Mesh_Data () = default;

// -----------------------------------------------------------------------------
        /**
        * @brief             This function updates the  mMeshData data struct.
        *                    It must be called before calculate_raycast()
        *                    or calculate_sdf().
        *
        */
        void update ();
// ----------------------------------------------------------------------------

/**
         *
         * @brief             return the number of elements on proc
         *
         */
        moris::uint
        get_number_of_elements() const
        {
            return mNumberOfElements;
        }

// -----------------------------------------------------------------------------
        /**
         *
         * @brief             return the number of nodes on proc
         *
         */
        moris::uint
        get_number_of_nodes() const
        {
            return mNumberOfNodes;
        }

// -----------------------------------------------------------------------------
        /**
         *
         * @brief             return the ids of nodes on proc
         *
         */
        moris::Matrix< moris::DDUMat >
        get_node_ids() const
        {
            return mNodesOnProc;
        }
   // -----------------------------------------------------------------------------
           /**
           * @brief             return the coordinates of a local node Ind
           *
           * @param[in] aNode   local node number on proc
           *
           */
        moris::Matrix< moris::DDRMat >
        get_node_coords(const moris::uint aLocalNodeInd) const
        {
            return mNodeCoords.get_column(aLocalNodeInd);
        }

// -----------------------------------------------------------------------------
        /**
        * @brief             return the local nodes of a local element Ind
        *
        * @param[in] aNode   local node indexon proc
        *
        */
        moris::Matrix< moris::DDUMat >
        get_nodes_of_element(const moris::uint aLocalNodeInd) const
        {

            moris::Matrix< moris::DDUMat > aNodes(mNumberOfNodesPerElement(aLocalNodeInd), 1);

            for (moris::uint k=0; k<mNumberOfNodesPerElement(aLocalNodeInd); ++k)
            {
                aNodes(k) = mElementTopology(k, aLocalNodeInd);
            }
            return aNodes;

        }

// -----------------------------------------------------------------------------

        /**
        * @brief             return the global ID of a node
        *
        * @param[in] aNode   local node index on proc
        *
        */
        moris::uint
        get_global_node_id(const moris::uint aLocalNodeInd) const
        {
            return mNodesOnProc(aLocalNodeInd);
        }

// -----------------------------------------------------------------------------

        /**
        * @brief             return the maximum of all element IDs
        *
        */
        moris::uint
        get_max_element_id() const
        {
            return mElementsOnProc.max();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief             return the global ID of an element
         *
         * @param[in] aNode   local element index on proc
         *
         */
        moris::uint
        get_global_element_id(const moris::uint aLocalElementInd) const
        {
            return mElementsOnProc(aLocalElementInd);
        }

// -----------------------------------------------------------------------------

        /**
        * @brief             return the local ID of a node
        *
        * @param[in] aNode   local node index on proc
        *
        */
        /*moris::uint
        get_local_node_id(const moris::uint aLocalNodeInd) const
        {
            return mLocalNodesOnProc(aLocalNodeInd);
        }*/

// -----------------------------------------------------------------------------

        /**
        * @brief             return the local ID of an element
        *
        * @param[in] aNode   local node number on proc
        *
        */
        /*moris::uint
        get_local_element_id(const moris::uint aLocalElementId) const
        {
            return mLocalElementsOnProc(aLocalElementId);
        }*/
// -----------------------------------------------------------------------------

        /**
        *
        * @brief             return the minimum coordinate of the mesh
        *
        */
        moris::Matrix< moris::DDRMat >
        get_min_coord() const
        {
            return mMinNodeCoordinate;
        }

// -----------------------------------------------------------------------------

        /**
        *
        * @brief             return the maximum coordinate of the mesh
        *
        */
        moris::Matrix< moris::DDRMat >
        get_max_coord() const
        {
            return mMaxNodeCoordinate;
        }

    };
}
#endif //MORIS_CL_GE_SDF_MESH_DATA_HPP
