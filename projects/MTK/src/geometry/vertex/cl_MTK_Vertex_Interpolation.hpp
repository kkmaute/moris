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

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CNT/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
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
            moris::Cell< mtk::Vertex* >    mCoefficients;
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
             * set the coefficient objects
             */
            virtual void
            set_coefficients( moris::Cell< Vertex* >& aCoefficients ) = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects
             */
            virtual moris::Cell< Vertex* >&
            get_coefficients() = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects (const version)
             */
            virtual const moris::Cell< Vertex* >&
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
    }         // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_VERTEX_INTERPOLATION_HPP_ */
