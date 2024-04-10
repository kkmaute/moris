/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Design_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{
    /**
     * This is a struct used to simplify \ref moris::gen::Design_Field constructors. It contains additional parameters that
     * are used by all fields, with given defaults.
     */
    struct Field_Parameters
    {
        sint         mDiscretizationIndex;         // Index of a mesh for discretization (-2 = none, -1 = store nodal values)
        real         mDiscretizationLowerBound;    // Lower bound for the B-spline coefficients in this field
        real         mDiscretizationUpperBound;    // Upper bound for the B-spline coefficients in this field
        bool         mUseMultilinearInterpolation; // Whether to use multilinear interpolation for all derived node field values

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Design field parameter list
         */
        explicit Field_Parameters( const Parameter_List& aParameterList );
    };

    class Design_Field
    {
      protected:
        std::shared_ptr< Field > mField;
        Node_Manager*            mNodeManager;

      private:
        Field_Parameters mParameters;
        Matrix< DDRMat > mInterpolatedSensitivities;
        Vector< sint > mInterpolatedADVIDs;

      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         * @param aNodeManager Node manager from the geometry engine, if available
         */
        Design_Field(
                std::shared_ptr< Field > aField,
                Field_Parameters         aParameters,
                Node_Manager&            aNodeManager );

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         */
        void discretize(
                mtk::Mesh_Pair        aMeshPair,
                sol::Dist_Vector*     aOwnedADVs,
                const Vector< sint >& aSharedADVIds,
                uint                  aADVOffsetID );

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMTKField Input MTK field to map based on
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         */
        void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs,
                const Vector< sint >&         aSharedADVIds,
                uint                          aADVOffsetID );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) const;

        /**
         * Gets the IDs of ADVs which this design component depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Vector< sint > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return d(field value)/d(ADV_j)
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADVs
         */
        template< typename Vector_Type >
        void set_advs( Vector_Type& aADVs )
        {
            mField->set_advs( aADVs );
        }

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs );

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        void reset_nodal_data( mtk::Interpolation_Mesh* aMesh );

        /**
         * Gets the name of this design field.
         *
         * @return Name
         */
        std::string get_name();

        /**
         * Gets whether this field will be using multilinear interpolation to get derived node field values.
         *
         * @return Multilinear interpolation flag
         */
        bool use_multilinear_interpolation() const;

        /**
         * Gets an MTK field, if this design field uses one that needs to be remapped to a new mesh
         *
         * @return
         */
        std::shared_ptr< mtk::Field > get_mtk_field();
    };
}    // namespace moris::gen
