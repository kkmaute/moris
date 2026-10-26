/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_VERTEX_INTERPOLATION_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_VERTEX_INTERPOLATION_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

//------------------------------------------------------------------------------

namespace moris::mtk
{
    //! forward declaration of vertex
    class Vertex;

    class Vertex_Interpolation
    {
        //------------------------------------------------------------------------------

      protected:
        // FIXME: MOVE THESE TO THE APPROPRIATE CHILD CLASSES
        // FIXME: BY HAVING MEMBER DATA HERE WE IMPOSE A DATA
        // FIXME: STRUCTURE ON ALL IMPLEMENTATIONS.
        Vector< mtk::Vertex* >         mCoefficients;
        moris::Matrix< moris::DDRMat > mWeights;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        Vertex_Interpolation(){};

        //------------------------------------------------------------------------------

        virtual ~Vertex_Interpolation(){};

        //------------------------------------------------------------------------------
        /**
         * returns the IDs of the interpolation coefficients
         */
        virtual Matrix< IdMat >
        get_ids() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the indices of the interpolation coefficients
         */
        virtual Matrix< IndexMat >
        get_indices() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the proc owners of the IDs of this vertex
         */
        virtual Matrix< IdMat >
        get_owners() const = 0;

        //------------------------------------------------------------------------------
        /**
         * set the interpolation weights
         */
        virtual void
        set_weights( const moris::Matrix< DDRMat >& aWeights ) = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the interpolation weights
         */
        virtual const Matrix< DDRMat >*
        get_weights() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns if interpolation contains extended basis
         */
        virtual bool
        contains_extended_basis() const 
        {
            return false;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the SPG for which a given extended basis is associated with. This is only relevant for extended basis functions
         */
        virtual Vector< Vector< moris_index > > const &
        get_basis_weight_numerators() const
        {
            MORIS_ERROR( 0, "get_basis_weight_numerators not implemented in base vertex interpolation" );
            static Vector< Vector< moris_index > > tVector;
            return tVector;
        }

        //------------------------------------------------------------------------------

        /**
         * returns all the SPGs a given extended basis' averaging weights are associated with. This is only relevant for extended basis functions
         */

        virtual Vector< Vector< moris_index > > const &
        get_basis_weight_denominators() const
        {
            MORIS_ERROR( 0, "get_basis_weight_denominators not implemented in base vertex interpolation" );
            static Vector< Vector< moris_index > >tVector;
            return tVector;
        }

        //------------------------------------------------------------------------------
        /**
         * returns the per-bad basis weights for a given extended basis. This is only relevant for extended basis functions
         */

        virtual Vector< Vector< real > > const &
        get_basis_weight_per_leader_basis() const
        {
            MORIS_ERROR( 0, "get_basis_weight_per_leader_basis not implemented in base vertex interpolation" );
            static Vector< Vector< real > > tVector;
            return tVector;
        }

        //------------------------------------------------------------------------------

        virtual Vector< Vector< uint > > const &
        get_basis_weight_denominator_count() const
        {
            MORIS_ERROR( 0, "get_basis_weight_denominator_count not implemented in base vertex interpolation" );
            static Vector< Vector< uint > > tVector;
            return tVector;
        }

        //------------------------------------------------------------------------------

        virtual Vector< Vector< uint > > const &
        get_basis_weight_numerator_count() const
        {
            MORIS_ERROR( 0, "get_basis_weight_numerator_count not implemented in base vertex interpolation" );
            static Vector< Vector< uint > > tVector;
            return tVector;
        }

        //------------------------------------------------------------------------------
        /**
         * set the coefficient objects
         */
        virtual void
        set_coefficients( Vector< Vertex* >& aCoefficients ) = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects
         */
        virtual Vector< Vertex* >&
        get_coefficients() = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects (const version)
         */
        virtual const Vector< Vertex* >&
        get_coefficients() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the number of coefficients attributed to this basis
         */
        virtual uint
        get_number_of_coefficients() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the ijk id for hm meshes
         */
        virtual Matrix< IdMat >
        get_ijkl_id() const
        {
            MORIS_ERROR( false, "Vertex_Interpolation::get_ijkl_id() not implemented for this mesh type." );
            return Matrix< IdMat >( 0, 0 );
        }

        //------------------------------------------------------------------------------

        virtual void
        print() const
        {
            MORIS_ERROR( false, "Vertex_Interpolation::print() - Not implemented for base class." );
        }

        //------------------------------------------------------------------------------

    };    // class Vertex_Interpolation
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_VERTEX_INTERPOLATION_HPP_ */
