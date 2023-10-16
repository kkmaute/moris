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
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Child_Node.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::ge
{
    /**
     * This is a struct used to simplify \ref moris::ge::Design_Field constructors. It contains additional parameters that
     * are used by all fields, with given defaults.
     */
    struct Field_Parameters
    {
        Cell< uint > mNumberOfRefinements;       //! The number of refinement steps to use for this field
        Cell< uint > mRefinementMeshIndices;     //! Indices of meshes to perform refinement on
        sint         mRefinementFunctionIndex;   //! Index of a user-defined refinement function (-1 = default)
        sint         mDiscretizationIndex;       //! Index of a mesh for discretization (-2 = none, -1 = store nodal values)
        real         mDiscretizationLowerBound;  //! Lower bound for the B-spline coefficients in this field
        real         mDiscretizationUpperBound ; //! Upper bound for the B-spline coefficients in this field

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Design field parameter list
         */
        explicit Field_Parameters( const ParameterList& aParameterList );
    };

    class Design_Field
    {
      private:
        std::shared_ptr< Field > mField;
        Field_Parameters         mParameters;

      public:

        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        Design_Field(
                std::shared_ptr< Field > aField,
                Field_Parameters         aParameters );

        /**
         * Gets if this field is to be used for seeding a B-spline field.
         *
         * @return Logic for B-spline creation
         */
        bool intended_discretization();

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aCoefficientIndices Coefficient indices to be mapped to
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         */
        void discretize(
                mtk::Mesh_Pair        aMeshPair,
                sol::Dist_Vector*     aOwnedADVs,
                const Matrix<DDUMat>& aCoefficientIndices,
                const Matrix<DDSMat>& aSharedADVIds,
                uint                  aADVOffsetID );

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
         * Gets the IDs of ADVs which this design component depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Returns if this design object has sensitivities, based on its ADV information.
         *
         * @return if this field depends on ADVs
         */
        bool depends_on_advs();

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
         * Add a new child node for evaluation, implemented for discrete integration fields.
         *
         * @param aNodeIndex Index of the child node
         * @param aChildNode Contains information about how the child node was created
         */
        void add_child_node(
                uint                          aNodeIndex,
                std::shared_ptr< Child_Node > aChildNode );

        /**
         * In relevant derived classes, uses additional information from the given interpolation mesh to define
         * potentially new nodes on the field. Implemented for discrete interpolation fields.
         *
         * @param aMesh Interpolation mesh with additional nodes
         */
        void add_nodal_data( mtk::Interpolation_Mesh* aMesh );

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        void reset_nodal_data();

        /**
         * Gets the name of this design's field
         *
         * @return Name
         */
        std::string get_name();

        /**
         * This function will return true when called less than the number of refinements set for this field,
         * and false otherwise. This is to determine for a given refinement call if this field needs refinement.
         *
         * @return if to perform an additional refinement with this field
         */
        const Cell< uint >& get_num_refinements();

        const Cell< uint >& get_refinement_mesh_indices();

        /**
         * Gets the index of a user-defined refinement function used within HMR.
         *
         * @return User-defined refinement function index
         */
        sint get_refinement_function_index();

        /**
         * Gets a discretization mesh index for a discretized field.
         *
         * @return Mesh index
         */
        moris_index get_discretization_mesh_index() const;

        /**
         * Gets the order of the discretization mesh.
         *
         * @return discretization order
         */
        uint get_discretization_order() const;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        real get_discretization_lower_bound();

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        real get_discretization_upper_bound();

        // TODO where should I really put this stuff?
        void set_num_original_nodes( uint aNumOriginalNodes );
        void set_mesh_pair( mtk::Mesh_Pair aMeshPair )
        {
            mField->mMeshPair = aMeshPair;
        }

        /**
         * Allows for access to the GEN field
         *
         * @return Underlying field
         */
        std::shared_ptr< Field > get_field()
        {
            return mField;
        }

    };
}
