/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Analytic.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"
#include "cl_MTK_Mapper.hpp"

namespace moris::gen
{
    /**
     * Class for fields that are specified by an analytic function.
     *
     * @tparam N Number of coordinates required to specify a center (or other reference point) for translations and rotations.
     * These need to be the first N parameters to an analytic field in order for field arrays and rotations to work.
     */
    template< uint N >
    class Field_Analytic : public Field
    {
      public:
        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         */
        Field_Analytic( ADV_ARG_TYPES )
                : Field( ADV_ARGS )
        {
        }

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        explicit Field_Analytic( const Vector< real >& aConstants )
                : Field( aConstants, "" )
        {
        }

        /**
         * Copy constructor with replacement variables for new constants.
         *
         * @param aCopy Analytic field to copy
         * @param aReplaceVariables Variable indices to replace
         * @param aNewConstants New constants
         */
        Field_Analytic(
                const Field_Analytic< N >& aCopy,
                const Vector< uint >& aReplaceVariables,
                const Vector< real >& aNewConstants )
                : Field( aCopy, aReplaceVariables, aNewConstants )
        {
        }

        /**
         * Gets the number of reference coordinates this field has.
         *
         * @return Number of reference coordinates
         */
        uint get_number_of_reference_coordinates() override
        {
            return N;
        }

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) final
        {
            return this->get_field_value( aCoordinates );
        }

        /**
         * Gets a field value of a derived node.
         *
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         * @return Field value
         */
        real get_field_value(
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) final
        {
            return this->get_field_value( aDerivedNode.get_global_coordinates() );
        }

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
                const Matrix< DDRMat >& aCoordinates ) final
        {
            return this->get_dfield_dadvs( aCoordinates );
        }

        /**
         * Gets a vector of the field derivatives with respect to ADVs of a derived node.
         *
         * @param aSensitivities Sensitivities to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_dfield_dadvs(
                Matrix< DDRMat >&   aSensitivities,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) final
        {
            aSensitivities = this->get_dfield_dadvs( aDerivedNode.get_global_coordinates() );
        }

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
                Matrix< DDRMat >&       aSensitivities ) final
        {
            this->get_dfield_dcoordinates( aCoordinates, aSensitivities );
        }

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

        /**
         * Returns a nullptr, since all analytic fields by definition do not need to be remapped
         *
         * @return nullptr
         */
        std::shared_ptr< mtk::Field > get_mtk_field() final
        {
            // TODO make this just return a nullptr once the refinement interface is finished, as an MTK field won't be needed then
            mtk::Mapper tMapper;
            std::shared_ptr< mtk::Field > tField = this->create_mtk_field( mMeshPairForAnalytic );
            tMapper.perform_mapping( tField.get(), mtk::EntityRank::NODE, mtk::EntityRank::BSPLINE );
            return tField;
        }
    };
}

/**
 * \def ANALYTIC_FIELD_DECLARATION( class_name, num_dimensions, num_variables, ... )
 * Automatically creates a constructor that has ADV arguments for general field creation, and an additional copy function
 *
 * @param class_name Name of the specific field constructor to create
 * @param num_dimensions Number of dimensions this field is defined with
 * @param num_variables Number of variables this field takes in
 * @param ... Additional code to be inserted into the constructor body via __VA_ARGS__
 */
#define ANALYTIC_FIELD_ADV_CONSTRUCTOR( class_name, num_dimensions, num_variables, ... )                            \
    class_name( ADV_ARG_TYPES )                                                                                     \
            : Field_Analytic< num_dimensions >( ADV_ARGS )                                                          \
    {                                                                                                               \
        VARIABLE_CHECK( num_variables );                                                                            \
        __VA_ARGS__                                                                                                 \
    }                                                                                                               \
    class_name( const class_name& aCopy, const Vector< uint >& aReplaceVariables, const Vector< real >& aNewConstants ) \
            : Field_Analytic< num_dimensions >( aCopy, aReplaceVariables, aNewConstants )                           \
    {                                                                                                               \
    }                                                                                                               \
    std::shared_ptr< Field > copy( const Vector< uint >& aReplaceVariables, const Vector< real >& aNewConstants )       \
    {                                                                                                               \
        return std::make_shared< class_name >( *this, aReplaceVariables, aNewConstants );                           \
    }
