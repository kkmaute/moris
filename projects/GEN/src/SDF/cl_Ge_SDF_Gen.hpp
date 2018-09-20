/*
 * cl+Ge_SDF_Gen.hpp
 *
 *  Created on: Mar 6, 2018
 *      Author: messe
 */

#ifndef SRC_GEOMENG_CL_GE_SDF_GEN_HPP_
#define SRC_GEOMENG_CL_GE_SDF_GEN_HPP_

#include <string>
#include "cl_Cell.hpp" // CON/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Database.hpp" // MTK/src
#include "cl_Ge_SDF_Mesh_Data.hpp"
#include "cl_Ge_SDF_Data.hpp"
#include "cl_Ge_SDF_Core.hpp"
#include "cl_Ge_SDF_Triangle_File.hpp"

namespace ge {
// -----------------------------------------------------------------------------
    class SDF_Gen
    {
        SDF_Core      mCore;
        SDF_Mesh_Data mMeshData;
        moris::Cell< SDF_Data > mObjects;
        moris::uint mNumberOfObjects;

// =============================================================================
    public:
// =============================================================================

    /**
      * @brief standard constructor for SDF Generator
      *
      * @param[in] aHMesh       reference to mesh object
      * @param[in] aFilePaths   cell containing paths to triangle files
      * @param[in] aObjectSDFs  cell of moris Mats containing the SDF for the object
      */
        SDF_Gen (
                const moris::database                    & aBackgroundMesh, // !< a wrapper of the mesh we are working with
                const moris::Cell<std::string>           & aFilePaths,      // !< file paths
                moris::Cell< moris::Matrix< moris::DDRMat > >  & aObjectSDFs,
                moris::Cell< moris::BoostBitset >        & aObjectSDFFlags      // !< individual SDFs for each object
        );      // !< individual SDFs for each object);    //!< Cell with file paths
// -----------------------------------------------------------------------------
        /**
          * @brief standard destructor for SDF Generator
          *
          */
        ~SDF_Gen () = default;

// -----------------------------------------------------------------------------

        /**
        * @brief             Calls the raycast algorithm for individual object.
        *
        */
        void calculate_raycast();
// -----------------------------------------------------------------------------
        /**
        * @brief             user setting
        *
        */
        void
        set_candidate_search_depth(const moris::uint aCandidateSearchDepth);

// -----------------------------------------------------------------------------

        /**
         * @brief             Internally calls the raycast algorithm and calculates
         *                    the signed distance field.
         *
         */
         void calculate_raycast_and_sdf();

// -----------------------------------------------------------------------------

         /**
          * @brief             Internally calls the raycast algorithm and calculates
          *                    the signed distance field for specified object.
          *
          */
         void
         calculate_raycast_and_sdf( const moris::uint aObject );

// -----------------------------------------------------------------------------

         /**
          * @brief             Internally calls the raycast algorithm for all objects
          *                    and calculates the signed distance field for the first
          *                    object only.
          *
          */
          void
          calculate_raycast_and_first_sdf();

// -----------------------------------------------------------------------------

/**
* @brief Returns moris::Mat<uint> containing all elements at the surface
*
*/
        moris::Matrix< moris::DDUMat >
        get_elements_at_surface();

        moris::Matrix< moris::DDUMat >
        get_elements_at_surface( const moris::uint aObject ) const;

// -----------------------------------------------------------------------------

/**
 * @brief Returns a moris::Mat<uint> containing all elements in the volume
 *
 */

        moris::Matrix< moris::DDUMat >
        get_node_ids();

// -----------------------------------------------------------------------------

        /**
         * @brief Returns the ids of all nmodes
         *
         */

        moris::Matrix< moris::DDUMat >
        get_elements_in_volume();

        moris::Matrix< moris::DDUMat >
        get_elements_in_volume( const uint aObject ) const;

// -----------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        get_inside_outside_sign_from_raycast( const uint aObject ) const;

// -----------------------------------------------------------------------------
        /**
         * @brief Returns moris::Mat<uint> containing all elements at the surface
         *                                 without those from first object
         *
         */
        moris::Matrix< moris::DDUMat >
        get_elements_at_surface_for_hmr_ref();

// -----------------------------------------------------------------------------

        /**
         * @brief Returns a moris::Mat<uint> containing all elements in the volume
         * *                                 without those from first object
         *
         */

        //moris::Matrix< moris::DDUMat >
        //get_elements_in_volume_for_HMR_ref();

// -----------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        get_node_signs();

    };
} /* namespace ge */

#endif /* SRC_GEOMENG_CL_GE_SDF_GEN_HPP_ */
