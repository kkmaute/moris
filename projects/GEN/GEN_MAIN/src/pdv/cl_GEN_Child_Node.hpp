/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Child_Node.hpp
 *
 */

#ifndef MORIS_CL_GEN_CHILD_NODE_HPP
#define MORIS_CL_GEN_CHILD_NODE_HPP

#include "cl_Matrix.hpp"
#include "cl_XTK_Basis_Function.hpp"
#include "cl_Cell.hpp"
#include "cl_GEN_Geometric_Proximity.hpp"
namespace moris
{
    namespace mtk
    {
        class Cell;
    }

    namespace ge
    {

        // Forward declaration of Field class
        class Field;

        //------------------------------------------------------------------------------

        enum class Element_Intersection_Type
        {
            Linear_1D,
            Linear_2D,
            Linear_3D
        };

        //------------------------------------------------------------------------------

        class Child_Node
        {
            friend class Geometry_Engine;

            //------------------------------------------------------------------------------

          protected:
            Matrix< DDUMat >         mAncestorNodeIndices;
            Cell< Matrix< DDRMat > > mAncestorNodeCoordinates;
            Matrix< DDRMat >         mBasisValues;

            //------------------------------------------------------------------------------

          private:
            Matrix< DDRMat > mJoinedSensitivities;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * Constructor
             *
             * @param aAncestorNodeIndices Node indices of the ancestors of this child node
             * @param aAncestorNodeCoordinates Global coordinates of the ancestors of this child node
             * @param aBasisFunction Basis function of the ancestor topology
             * @param aLocalCoordinatesInAncestor Local coordinate of this child inside of the ancestor element
             */
            Child_Node(
                    Matrix< DDUMat >                aAncestorNodeIndices,
                    Cell< Matrix< DDRMat > >        aAncestorNodeCoordinates,
                    const Element_Intersection_Type aBasisFunction,
                    Matrix< DDRMat >                aLocalCoordinatesInAncestor );

            //------------------------------------------------------------------------------

            Child_Node(
                    moris::mtk::Cell* aCell,
                    Matrix< DDRMat >* aLocalCoordinates,
                    bool              aEvaluateAsLinear );

            //------------------------------------------------------------------------------

            /**
             * Get the field value on the child node based on values from its ancestors.
             *
             * @param aField Field pointer, referenced during call from field class
             * @return Field value
             */
            virtual real interpolate_field_value( Field* aField );

            //------------------------------------------------------------------------------

            /**
             * Joins the field sensitivities on the child node based on its ancestors.
             *
             * @param aField Field pointer, referenced during call from field class
             * @return Field sensitivities
             */
            virtual const Matrix< DDRMat >& join_field_sensitivities( Field* aField );

            //------------------------------------------------------------------------------

            /**
             * Joins the depending ADV IDs on the child node based on its ancestors.
             *
             * @param aField Field pointer, referenced during call from field class
             * @return Field ADV IDs
             */
            virtual Matrix< DDSMat > join_determining_adv_ids( Field* aField );

            //------------------------------------------------------------------------------

            /**
             * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aField Field pointer, referenced during call from field class
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            virtual void get_dfield_dcoordinates(
                    Field*            aField,
                    Matrix< DDRMat >& aSensitivities );

            //------------------------------------------------------------------------------

        };    // class Child_Node

        //------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_CHILD_NODE_HPP
