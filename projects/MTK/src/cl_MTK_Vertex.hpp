/*
 * cl_MTK_Vertex.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_VERTEX_HPP_
#define SRC_MESH_CL_MTK_VERTEX_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex(){};

//------------------------------------------------------------------------------

            /**
             * Destructor, virtual
             */
            virtual
            ~Vertex(){};

//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            virtual Matrix< DDRMat >
            get_coords() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return Matrix < DDRMat >(0,0);
            }

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_id
            get_id() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return 0;
            }

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_index
            get_index() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return 0;
            }

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this vertex
             */
            virtual Matrix< IdMat > const &
            get_adof_ids() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return mDummyMat;
            }

//------------------------------------------------------------------------------

            virtual Matrix< IndexMat >
            get_adof_indices() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return mDummyMat;
            }

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this vertex
             */
            virtual Matrix< IdMat >
            get_adof_owners() const
            {
                MORIS_ERROR(0,"Function not implemented in base vertex");
                return mDummyMat;
            }

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this vertex
             */
            virtual moris::Cell< Vertex* > &
            get_adof_pointers()
            {
                return mDummyAdofs;
            }

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this vertex
             */
            virtual const  moris::Cell< Vertex* > &
            get_adof_pointers()  const
            {
                return mDummyAdofs;
            }

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this vertex
             */
            virtual const Matrix< DDRMat > *
            get_t_matrix() const
            {
                return &mDummyRealMat;
            }

//------------------------------------------------------------------------------

            /**
             * returns the id of the proc that owns this vertex
             */
            virtual moris_id
            get_owner() const
            {
                return 0;
            }


            Matrix< IndexMat > mDummyMat;
            Matrix< DDRMat > mDummyRealMat;

            // TODO MOVE ADOF RELATED STUFF OUT OF VERTEX
            moris::Cell< Vertex* > mDummyAdofs;


//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_VERTEX_HPP_ */
