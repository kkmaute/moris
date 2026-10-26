/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_DataBase.hpp
 *
 */

#ifndef SRC_cl_MTK_Vertex_Interpolation_DataBase_HPP_
#define SRC_cl_MTK_Vertex_Interpolation_DataBase_HPP_

#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris::mtk
{
    class Mesh;

    class Vertex_Interpolation_DataBase : public mtk::Vertex_Interpolation
    {
      private:
        moris_index mVertexIndex;
        moris_index mVertexOrder;
        mtk::Mesh*  mMesh = nullptr;
        Vector< Vector< moris_index > > mWeightNumerators;
        Vector< Vector< moris_index > > mWeightDenominators;
        Vector< Vector< moris::real > > mWeightsPerLeader;
        Vector< Vector< moris::uint > > mWeightNumeratorIndices;
        Vector< Vector< moris::uint > > mWeightDenominatorIndices;
        bool       mHasExtendedBasis = false;

      public:
        // default constrctor
        Vertex_Interpolation_DataBase() = default;

        //------------------------------------------------------------------------------

        Vertex_Interpolation_DataBase( moris_index aVertexIndex,
                moris_index                        aVertexOrder,
                mtk::Mesh*                         aMesh );

        //------------------------------------------------------------------------------

        ~Vertex_Interpolation_DataBase() override = default;

        //------------------------------------------------------------------------------
        /**
         * returns the IDs of the interpolation coefficients
         */
        Matrix< IdMat >
        get_ids() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the indices of the interpolation coefficients
         */
        Matrix< IndexMat >
        get_indices() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the proc owners of the IDs of this vertex
         */
        Matrix< IdMat >
        get_owners() const override;

        //------------------------------------------------------------------------------
        /**
         * set the interpolation weights
         */
        void
        set_weights( const moris::Matrix< DDRMat >& aWeights ) override;

        //------------------------------------------------------------------------------

        /**
         * set the basis extension flag
         */
        void
        set_basis_extension_flag( bool aHasExtendedBasis );

        //------------------------------------------------------------------------------

        bool
        contains_extended_basis() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the interpolation weights
         */
        const Matrix< DDRMat >*
        get_weights() const override;

        //------------------------------------------------------------------------------
        /**
         * returns the numerators associated with basis weights (when basis extension is being used)
         */
        const Vector< Vector< moris_index > >&
        get_basis_weight_numerators() const override;

        //------------------------------------------------------------------------------
        /**
         * returns the denominators associated with basis weights (when basis extension is being used)
         */
        const Vector< Vector< moris_index > >&
        get_basis_weight_denominators() const override;

        //------------------------------------------------------------------------------
        /**
         * returns the per-leader weights associated with basis weights (when basis extension is being used )
         */
        const Vector< Vector< real > >&
        get_basis_weight_per_leader_basis() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the numerator indices associated with basis weights (when basis extension is being used)
         */
        const Vector< Vector< moris::uint > >&
        get_basis_weight_numerator_count() const override;

        //------------------------------------------------------------------------------
        /**
         * returns the denominator indices associated with basis weights (when basis extension is being used)
         */
        const Vector< Vector< moris::uint > >&
        get_basis_weight_denominator_count() const override;

        //------------------------------------------------------------------------------
        /**
         * set the coefficient objects
         */
        void
        set_coefficients( Vector< Vertex* >& aCoefficients ) override;

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects
         */
        Vector< Vertex* >&
        get_coefficients() override;

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects (const version)
         */
        const Vector< Vertex* >&
        get_coefficients() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the number of coefficients attributed to this basis
         */
        uint
        get_number_of_coefficients() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief populate the data that is being passed as a refernce
         *
         */

        void
        set_outward_data();

        //------------------------------------------------------------------------------

        /**
         * @brief memory usage of the vertex interpolations
         *
         * @return size_t
         */

        size_t
        capacity();
    };
}    // namespace moris::mtk

#endif /* cl_MTK_Vertex_Interpolation_DataBase.hpp */
