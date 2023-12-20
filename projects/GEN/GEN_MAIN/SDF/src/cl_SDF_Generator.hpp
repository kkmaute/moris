/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Generator.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_GENERATOR_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_GENERATOR_HPP_

#include <string>
#include <memory>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_SDF_Object.hpp"
namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class SDF_Generator
        {
            //! file containing object data
            Object mObject;

            //! verbosity flag
            bool          mVerboseFlag = false;
//-------------------------------------------------------------------------------
public:
//-------------------------------------------------------------------------------
            /**
             * constructor with pointer. Creates an SDF Object from aObjectPath.
             */

            SDF_Generator( const std::string & aObjectPath,
                           const bool aVerboseFlag = true );

            SDF_Generator( const std::string & aObjectPath,
                           Matrix< DDRMat >&   aObjectOffset,
                           const bool aVerboseFlag = true );

//-------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~SDF_Generator(){};

//-------------------------------------------------------------------------------

            /**
             * performs a raycast for an MTK mesh and returns Matrices with
             * element IDs and indices
             */
            void
            raycast(  mtk::Mesh          * aMesh,
                    Matrix< IndexMat > & aElementsAtSurface );

//-------------------------------------------------------------------------------

            /**
             * performs a raycast for an MTK mesh and returns Matrices with
             * element IDs and indices ( shared pointer version )
             */
            void
            raycast( std::shared_ptr< mtk::Mesh > aMesh,
                    Matrix< IndexMat >            & aElementsAtSurface );

//-------------------------------------------------------------------------------

            /**
             * performs a raycast for an MTK mesh and returns Matrices with
             * element IDs and indices
             */
            void
            raycast(  mtk::Mesh        * aMesh,
                    Matrix< IndexMat > & aElementsAtSurface,
                    Matrix< IndexMat > & aElementsInVolume );

//-------------------------------------------------------------------------------

            /**
             * performs a raycast for an MTK mesh and returns Matrices with
             * element IDs and indices ( shared pointer version )
             */
            void
            raycast( std::shared_ptr< mtk::Mesh > aMesh,
                    Matrix< IndexMat >           & aElementsAtSurface,
                    Matrix< IndexMat >           & aElementsInVolume );

//-------------------------------------------------------------------------------

            /**
             * calculates the SDF for a given mesh
             */
            void
            calculate_sdf(
                    mtk::Mesh          * aMesh,
                    Matrix< DDRMat>    & aSDF );

//-------------------------------------------------------------------------------

            /**
             * calculates the SDF for a given mesh ( shared pointer version )
             */
            void
            calculate_sdf(
                    std::shared_ptr< mtk::Mesh > aMesh,
                    Matrix< DDRMat>              & aSDF );

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_GENERATOR_HPP_ */

