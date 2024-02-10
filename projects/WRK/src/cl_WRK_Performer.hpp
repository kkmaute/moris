/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Performer.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

// This can be removed once we have the refinement interface
namespace moris::gen
{
    enum class Geometric_Region : signed char;
}

namespace moris::wrk
{
    class Performer
    {
        public:

            /**
             * Constructor
             */
            Performer()=default;

            /**
             * Destructor
             */
            ~Performer()=default;

            /**
             * Return the number of fields that can be used for refinement
             *
             * @return Number of fields for refinement
             */
            virtual uint get_num_refinement_fields() = 0;

            /**
             * Returns fields so that HMR can perform refinement based on the data from this performer
             *
             * @param aFieldIndex Index of the field
             * @param aNodeIndex Index of the node
             * @param aCoordinates Coordinates of the node
             */
            virtual gen::Geometric_Region get_geometric_region(
                    uint                  aFieldIndex,
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates) = 0;

            virtual const Vector< uint >& get_num_refinements(uint aFieldIndex ) = 0;

            virtual const Vector< uint >& get_refinement_mesh_indices(uint aFieldIndex ) = 0;

            /**
             * Gets the index of an HMR user-defined refinement function for the given field index
             *
             * @param aFieldIndex Index of the field
             * @param aRefinementIndex The current refinement step being performed
             * @return User-defined function index, or -1 to use default refinement
             */
            virtual sint get_refinement_function_index(
                    uint aFieldIndex,
                    uint aRefinementIndex) = 0;
    };
}
