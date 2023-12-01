/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_VERTEX_HPP_
#define SRC_MESH_CL_MTK_VERTEX_HPP_

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CNT/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        class Vertex
        {
          protected:
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
            virtual ~Vertex() = default;

            //------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            virtual Matrix< DDRMat >
            get_coords() const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );
                return Matrix< DDRMat >( 0, 0 );
            }

            //------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            virtual moris_id
            get_id() const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return gNoID;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the processor unique index of this vertex
             */
            virtual moris_index
            get_index() const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return gNoIndex;
            }

            //------------------------------------------------------------------------------

            /**
             * @brief returns the base vertex of an mtk vertex. 
             * Specifically, in XTK enriched meshes this returns the vertex index in the background;
             * and in VIS it returns the MTK vertex a given VIS vertex was constructed from 
             */
            virtual mtk::Vertex const *
            get_base_vertex() const
            {
                return this;
            }

            //------------------------------------------------------------------------------

            // fixme: change this into moris_id
            virtual moris_index
            get_owner() const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return 0;
            }

            //------------------------------------------------------------------------------

            virtual Vertex_Interpolation*
            get_interpolation( const uint aBSplineMeshIndex )
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return nullptr;
            }

            //------------------------------------------------------------------------------

            virtual const Vertex_Interpolation*
            get_interpolation( const uint aBSplineMeshIndex ) const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return nullptr;
            }

            virtual bool
            has_interpolation( const uint aBSplineMeshIndex )
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );

                return false;
            }

            //------------------------------------------------------------------------------

            virtual uint
            get_level() const
            {
                return 0;
            }

            //------------------------------------------------------------------------------

            virtual void
            flag()
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );
            }

            //------------------------------------------------------------------------------

            virtual void
            unflag()
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );
            }

            //------------------------------------------------------------------------------

            bool
            is_flagged() const
            {
                MORIS_ERROR( false, "Function not implemented in base vertex" );
                return false;
            }

            bool
            operator==( const mtk::Vertex& aVertex ) const
            {
                return this->get_id() == aVertex.get_id();
            }

            //------------------------------------------------------------------------------

        }; // class mtk::Vertex

        //------------------------------------------------------------------------------

    } /* namespace mtk */

    //------------------------------------------------------------------------------

    // operators for printing
    inline std::ostream&
    operator<<( std::ostream& os, const mtk::Vertex& dt )
    {
        os << "Vertex Id: " << dt.get_id() << " | Vertex Index: " << dt.get_index();

        return os;
    }

    //------------------------------------------------------------------------------
    
    inline std::ostream&
    operator<<( std::ostream& os, mtk::Vertex const * const & dt )
    {
        os << "Vertex Id: " << dt->get_id() << " | Vertex Index: " << dt->get_index();

        return os;
    }

    //------------------------------------------------------------------------------
    
    inline bool
    comparePtrToVertexIdBased(
            moris::mtk::Vertex* a,
            moris::mtk::Vertex* b )
    {
        return ( a->get_id() < b->get_id() );
    }

    //------------------------------------------------------------------------------

} /* namespace moris */

//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_VERTEX_HPP_ */
