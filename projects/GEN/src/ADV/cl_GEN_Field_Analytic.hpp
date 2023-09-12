/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Analytic.hpp
 *
 */

#ifndef MORIS_CL_GEN_FIELD_ANALYTIC_HPP
#define MORIS_CL_GEN_FIELD_ANALYTIC_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Analytic : virtual public Field
        {
          public:
            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Analytic();

            /**
             * Trivial destructor, necessary to avoid memory leaks
             */
            virtual ~Field_Analytic();

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aCoordinates vector of coordinate values
             * @return Field value
             */
            virtual real get_field_value( const Matrix< DDRMat >& aCoordinates ) = 0;

            /**
             * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Vector of sensitivities
             */
            const Matrix< DDRMat >& get_dfield_dadvs(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates );

            /**
             * Given a node coordinate, returns a vector of the field derivatives with respect to its ADVs.
             *
             * @param aCoordinates Vector of coordinate values
             * @return Vector of sensitivities
             */
            virtual const Matrix< DDRMat >& get_dfield_dadvs( const Matrix< DDRMat >& aCoordinates ) = 0;

            /**
             * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            void get_dfield_dcoordinates(
                    uint                    aNodeIndex,
                    const Matrix< DDRMat >& aCoordinates,
                    Matrix< DDRMat >&       aSensitivities );

            /**
             * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            virtual void get_dfield_dcoordinates(
                    const Matrix< DDRMat >& aCoordinates,
                    Matrix< DDRMat >&       aSensitivities ) = 0;
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_FIELD_ANALYTIC_HPP
