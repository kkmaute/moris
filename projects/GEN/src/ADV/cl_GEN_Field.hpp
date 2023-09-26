/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

namespace moris::ge
{
    class Field : public ADV_Manager
    {
      public:
        mtk::Mesh_Pair mMeshPair = mtk::Mesh_Pair( nullptr, nullptr );

      protected:
        uint mNumOriginalNodes = 0;

        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         */
        template< typename Vector_Type >
        Field( Vector_Type&     aADVs,
               Matrix< DDUMat > aFieldVariableIndices,
               Matrix< DDUMat > aADVIndices,
               Matrix< DDRMat > aConstants )
                : ADV_Manager( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
        {
        }

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        explicit Field( Matrix< DDRMat > aConstants );

        /**
         * Constructor that sets all field variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aFieldVariableIndices Field variable indices for assigning the shared ADV IDs
         * @param aSharedADVIds Shared ADV IDs needed for this field
         */
        Field( const Matrix< DDUMat >&  aFieldVariableIndices,
               const Matrix< DDSMat >&  aSharedADVIds );

      public:

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        virtual real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return d(field value)/d(ADV_j)
         */
        virtual const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) = 0;

        /**
         * Sets the dependencies of this field after they have been found by the owning property. By default
         * does nothing.
         *
         * @param aDependencyFields Other fields that this field depends on.
         */
        virtual void set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields );

        /**
         * Add a new child node for evaluation, implemented for discrete integration fields.
         *
         * @param aNodeIndex Index of the child node
         * @param aChildNode Contains information about how the child node was created
         */
        virtual void add_child_node(
                uint                          aNodeIndex,
                std::shared_ptr< Child_Node > aChildNode );

        /**
         * In relevant derived classes, uses additional information from the given interpolation mesh to define
         * potentially new nodes on the field. Implemented for discrete interpolation fields.
         *
         * @param aMesh Interpolation mesh with additional nodes
         */
        virtual void add_nodal_data( mtk::Interpolation_Mesh* aMesh );

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        virtual void reset_nodal_data();

        void set_num_original_nodes( uint aNumOriginalNodes );

    };
}
