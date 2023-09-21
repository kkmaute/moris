/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Multigeometry.hpp
 *
 */

#pragma once

#include "cl_GEN_Level_Set_Geometry.hpp"

namespace moris::ge
{
    class Multigeometry : public Field
    {

    private:
        Cell<std::shared_ptr< Field > > mGeometries;

    public:

        /**
         * Multigeometry constructor
         *
         * @param aGeometries Created geometries
         */
        Multigeometry(Cell<std::shared_ptr< Field > > aGeometries);

        /**
         * Given a node coordinate, the geometry needs to return the distance to the nearest function.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates vector of coordinate values
         * @return distance to nearest function
         */
        real get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinate, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities);

        /**
         * Adds a geometry to this multigeometry.
         *
         * @param aGeometry Geometry to add
         */
        void add_geometry(std::shared_ptr< Field > aGeometry);

    };
}
