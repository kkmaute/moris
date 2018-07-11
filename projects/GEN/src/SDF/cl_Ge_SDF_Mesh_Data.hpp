//
// Created by messe on 3/7/18.
//

#ifndef MORIS_CL_GE_SDF_MESH_DATA_HPP
#define MORIS_CL_GE_SDF_MESH_DATA_HPP

#include <map>
#include "chronos.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Map.hpp" // CON/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Sp_Mat.hpp" // LNA/src
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
        moris::Mat< moris::real > mNodeCoords;     //!< contains the coordinates of the nodes
        moris::Mat< moris::uint > mNumberOfNodesPerElement;
        moris::Mat< moris::uint > mElementTopology;

        moris::uint mNumberOfNodes = 0;             //!< number of nodes in the mesh on current proc
        moris::uint mNumberOfElements = 0;          //!< number of cells in the mesh on current proc
        moris::Mat< moris::real > mMinNodeCoordinate; //!< minimum coordinates of the mesh
        moris::Mat< moris::real > mMaxNodeCoordinate; //!< minimum coordinates of the mesh

        moris::Mat< moris::uint > mNodesOnProc;
        moris::Mat< moris::uint > mElementsOnProc;
        //moris::Mat< moris::uint > mLocalNodesOnProc;
        //moris::Mat< moris::uint > mLocalElementsOnProc;
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
        moris::Mat< moris::uint >
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
        moris::Mat<moris::real>
           get_node_coords(const moris::uint aLocalNodeInd) const
        {
            return mNodeCoords.cols(aLocalNodeInd, aLocalNodeInd);
        }

// -----------------------------------------------------------------------------
        /**
        * @brief             return the local nodes of a local element Ind
        *
        * @param[in] aNode   local node indexon proc
        *
        */
        moris::Mat<moris::uint>
        get_nodes_of_element(const moris::uint aLocalNodeInd) const
        {

            moris::Mat<moris::uint> aNodes(mNumberOfNodesPerElement(aLocalNodeInd), 1);

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
        moris::Mat< moris::real >
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
        moris::Mat< moris::real >
        get_max_coord() const
        {
            return mMaxNodeCoordinate;
        }

    };
}
#endif //MORIS_CL_GE_SDF_MESH_DATA_HPP
