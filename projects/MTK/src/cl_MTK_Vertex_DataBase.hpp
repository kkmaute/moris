/*
 * cl_MTK_Vertex_DataBase.hpp
 *
 *  Created on: Dec  11, 2021
 *      Author: momo
 */
#ifndef SRC_cl_MTK_Vertex_DataBase_HPP_
#define SRC_cl_MTK_Vertex_DataBase_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
namespace mtk
{
    // class Mesh;
    // class Interpolation_Mesh_Analysis;

    class Vertex_DataBase : public mtk::Vertex
    {
      private:
         moris_index mVertexIndex;
        mtk::Mesh*  mMesh = nullptr;


      public:
        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */

        Vertex_DataBase() = default;

        //------------------------------------------------------------------------------

        /**
         * @brief Construct a new Vertex_DataBase object // this is used for lagrange interpolation vertices
         *
         * @param aVertex an mtk vertex to get id,index,owner
         * @param aCoordinatePtr  coordinate pointer
         * @param mVertexInterpolationPtr pointer to the index of the
         * @param aDim spatial dimension  of the coords
         */
        Vertex_DataBase( moris_index aVertexIndex,
            mtk::Mesh* const&        aMesh );


        //------------------------------------------------------------------------------

        /**
         * Destructor, virtual
         */

        ~Vertex_DataBase() = default ;

        //------------------------------------------------------------------------------


        /**
         * @brief Get the coords of the vertex
         *
         * @return Matrix< DDRMat > coordinate matrix
         */
        virtual Matrix< DDRMat >
        get_coords() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the id of the object
         *
         * @return moris_id id of the vertex
         */
        virtual moris_id
        get_id() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the index of the object
         *
         * @return moris_index index of the vertex
         */
        virtual moris_index
        get_index() const override;


        //------------------------------------------------------------------------------

        // FIXME: change this into moris_id

        /**
         * @brief Get the owner of the object
         *
         * @return moris_index owner of the vertex
         */
        virtual moris_index
        get_owner() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the interpolation vertex of the lagrange vertex
         *
         * @param aBSplineMeshIndex  a bspline mesh index
         * @return Vertex_Interpolation* a bspline vertex
         */
        Vertex_Interpolation*
        get_interpolation( const uint aBSplineMeshIndex ) override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the constant interpolation vertex of the lagrange vertex
         *
         * @param aBSplineMeshIndex  a bspline mesh index
         * @return const Vertex_Interpolation* a bspline vertex
         */

        const Vertex_Interpolation*
        get_interpolation( const uint aBSplineMeshIndex ) const override;

        //------------------------------------------------------------------------------

        /**
         * @brief to indicate if it has an associated bsplie veretx or not
         *
         * @param aBSplineMeshIndex a bspline mesh index
         * @return true
         * @return false
         */

        virtual bool
        has_interpolation( const uint aBSplineMeshIndex ) override;


        //------------------------------------------------------------------------------

        /**
         * @brief overlead of operator == to see if two vertices have the same id
         *
         * @param aVertex
         * @return true
         * @return false
         */
        bool
        operator==( const mtk::Vertex& aVertex ) const
        {
            return this->get_id() == aVertex.get_id();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief memory usage of the mesh
         * 
         * @return size_t 
         */

        size_t
        capacity();
    };

}// namespace mtk
}// namespace moris


#endif /* cl_MTK_Vertex_DataBase.hpp */