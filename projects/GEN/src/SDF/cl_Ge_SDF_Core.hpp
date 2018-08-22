/*
 * SDF_Core.hpp
 *
 *  Created on: Feb 28, 2018
 *      Author: messe
 *
 *
 *
 */

#ifndef SRC_GEOMENG_CL_GE_SDF_CORE_HPP_
#define SRC_GEOMENG_CL_GE_SDF_CORE_HPP_

#ifdef MORIS_USE_ARMA
#include <armadillo>
#else
#include <Eigen/Core>
#endif

#include "chronos.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "fn_find.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "GeUtilities.hpp"
#include "cl_Ge_SDF_Triangle.hpp"
#include "cl_Ge_SDF_Mesh_Data.hpp"
#include "cl_Ge_SDF_Data.hpp"

namespace ge {
// -----------------------------------------------------------------------------

    /**
      * @brief The core object which calculates the signed distance field
      */
    class SDF_Core
    {

// =============================================================================
    public:
// =============================================================================
        /**
        * @brief                       constructor for SDF generator.
        *
        */
        SDF_Core ()
        {
        }
// -----------------------------------------------------------------------------
        /**
         * @brief default destructor for SDF generator.
         */
        ~SDF_Core () = default;

// =============================================================================
//  Public Callers
// =============================================================================

        /**
        * @brief             run the ray cast algorithm
        */
        void
        calculate_raycast (
                const ge::SDF_Mesh_Data &aMeshData,
                ge::SDF_Data &aSDFData);
// -----------------------------------------------------------------------------
        /**
* @brief             run noth the the ray cast and the sdf algorithm
*/
        void
        calculate_raycast_and_sdf(
                const ge::SDF_Mesh_Data &aMeshData,
                ge::SDF_Data &aData);

// -----------------------------------------------------------------------------

        /**
         * @brief              Saves data of proc to VTK. Can also visualize HEX27
         */
        void
        save_to_vtk ( const SDF_Mesh_Data &aMeshData,
                      const SDF_Data &aData,
                      const std::string &aFilePath);
// =============================================================================
//  Private Functions
// =============================================================================

        /**
          * @brief  voxelizing subroutine
          *
          * @param[in] aAxis               projection axis
          * @param[in] aPreselectFunction  pointer to corresponding preselection function
          */
        void
        voxelize (const SDF_Mesh_Data &aMeshData,
                       SDF_Data &aOutputData,
                       const moris::uint aAxis);

// -----------------------------------------------------------------------------

        /**
          * @brief  subroutine for triangle preselection along x-axis
          *
          * @param[in,out] aData       struct containing calculation data
          * @param[in]     aNodeCoords coordinate of node to be projected
          *
          */
        void
        preselect_triangles_x (SDF_Data &aData,
                                    const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------

        /**
          * @brief  subroutine for triangle preselection along y-axis
          *
          * @param[in,out] aData       struct containing calculation data
          * @param[in]     aNodeCoords coordinate of node to be projected
          *
          */
        void
        preselect_triangles_y (SDF_Data &aData,
                                    const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------

        /**
          * @brief  subroutine for triangle preselection along z-axis
          *
          * @param[in,out] aData       struct containing calculation data
          * @param[in]     aNodeCoords coordinate of node to be projected
          *
          */
        void
        preselect_triangles_z (SDF_Data &aData,
                                    const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------
        /**
          * @brief  subroutine for triangle intersection in voxelize algorithm
          *
          * @param[in,out] aData       struct containing calculation data
          * @param[in]     aAxis       projection axis
          * @param[in]     aNodeCoords coordinate of node to be projected
          *
          */
        void
        intersect_triangles (
                ge::SDF_Data &aData,
                const moris::uint aAxis,
                const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------
        /**
         * @brief  subroutine for triangle intersection
         *
         * @param[in,out] aData       struct containing calculation data
         * @param[in]     aAxis       projection axis
         * @param[in]     aNodeCoords coordinate of node to be projected
         *
         */
        void
        intersect_ray_with_triangles (
                ge::SDF_Data &aData,
                const moris::uint aAxis,
                const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------

        /**
          * @brief  subroutine for voxelizing
          *
          * @param[in,out] aMeshData       struct containing mesh information
          * @param[in,out] aData       struct containing calculation data
          * @param[in]     aAxis       projection axis
          * @param[in]     aNode       local number of the Node
          * @param[in]     aNodeCoords coordinates of the node
          */
        void
        check_if_node_is_inside (
                const SDF_Mesh_Data &aMeshData,
                SDF_Data &aData,
                const moris::uint aAxis,
                const moris::uint aLocalNodeInd,
                const moris::Mat<moris::real> &aNodeCoords);

// ----------------------------------------------------------------------------

/**
        * @brief Compute candidate points based on intersection:
        *        Calculating the distance to the triangles is very expensive.
        *        We can save a lot of effort by preselecting points near the
        *        surface, and, if mSettings.mCandidateSearchDepth > 0, in the
        *        near neighborhood of the surface. Nodes near the surface are
        *        detected by element-wise checking for sign changes.
        */
        void calculate_candidate_points_and_buffer_diagonal (
                const SDF_Mesh_Data &aMeshData,
                SDF_Data &aData);

// ----------------------------------------------------------------------------

/**
        *
        * @brief subroutine for candidate point selection
        *
        */
        moris::real get_diagonal_length_of_element (
                const ge::SDF_Mesh_Data &aMeshData,
                const moris::Mat<moris::uint> &aNodesOfElement) const;

// ----------------------------------------------------------------------------
    /**
     *
     *
     * @brief subroutine for signed distance field calculation
     *
     */
    void calculate_udf(const SDF_Mesh_Data& aMeshData,
                             SDF_Data& aData);

// ----------------------------------------------------------------------------
        /**
         *
         * @brief subroutine for unsigned distance field calculation
         *
         */
        moris::Mat<moris::uint> get_nodes_within_triangle_bounding_box(
                const SDF_Mesh_Data& aMeshData,
                SDF_Data& aOutputData,
                const moris::uint aTriangle);

// ----------------------------------------------------------------------------
        /**
         *
         *
         * @brief subroutine for signed distance field calculation
         *
         */
        void sign_udf(const SDF_Mesh_Data& aMeshData,
                           SDF_Data& aData);

    };
} /* namespace ge */

#endif /* SRC_GEOMENG_CL_GE_SDF_CORE_HPP_ */
