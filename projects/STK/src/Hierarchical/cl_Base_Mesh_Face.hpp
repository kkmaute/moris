/*
 * cl_Base_Mesh_Face.hpp
 *
 *  Created on: Feb 27, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_BASE_MESH_FACE_HPP_
#define SRC_MESH_CL_BASE_MESH_FACE_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src

namespace moris
{

    class Base_Mesh_Face
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Base_Mesh_Face()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Base_Mesh_Face() = default;

        /**
          * Provides the number of faces in x-direction
          *
          * @param[in] aLevel              Level of the element.
          * @param[in] aModelDim                Number of dimensions.
          * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
          *
          * @param[out] Face_number        Number of faces in x-direction.
          *
          */
         uint
         give_number_of_faces_x(
                 uint const & aLevel,
                 uint const & aModelDim,
                 Mat<uint> const & aNumberOfElementsPerDirection) const;

         /**
           * Provides the number of faces in y-direction
           *
           * @param[in] aLevel              Level of the element.
           * @param[in] aModelDim                Number of dimensions.
           * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
           *
           * @param[out] Face_number        Number of faces in y-direction.
           *
           */
          uint
          give_number_of_faces_y(
                  uint const & aLevel,
                  uint const & aModelDim,
                  Mat<uint> const & aNumberOfElementsPerDirection) const;

          /**
            * Provides the number of faces in z-direction
            *
            * @param[in] aLevel              Level of the element.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
            *
            * @param[out] Face_number        Number of faces in z-direction.
            *
            */
           uint
           give_number_of_faces_z(
                   uint const & aLevel,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
             * Provides the number of faces
             *
             * @param[in] aLevel              Level of the element.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             *
             * @param[out] Face_number        Number of faces in z-direction.
             *
             */
            uint
            give_number_of_faces(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the level of the face from the hierarchical mesh.
             *
             * @param[in] aFaceId            Face number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             *
             * @param[out] level        Level of the face from the hierarchical mesh.
             *
             */
            uint
            give_face_level(
                    uint const & aFaceId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
             * Provides the face number in x-direction of the position i,j,k
             *
             * @param[in] aLevel              Level of the element.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             * @param[in] aIJKPosition        Position I,J,K.
             *
             * @param[out] Face_number        Number of face in x-direction.
             *
             */
           uint
           give_face_x_of_position(
                   uint const & aLevel,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<uint> const & aIJKPosition) const;

           /**
             * Provides the face number in y-direction of the position i,j,k
             *
             * @param[in] aLevel              Level of the element.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             * @param[in] aIJKPosition        Position I,J,K.
             *
             * @param[out] Face_number        Number of face in y-direction.
             *
             */
           uint
           give_face_y_of_position(
                   uint const & aLevel,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<uint> const & aIJKPosition) const;

           /**
             * Provides the face number in z-direction of the position i,j,k
             *
             * @param[in] aLevel              Level of the element.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             * @param[in] aIJKPosition        Position I,J,K.
             *
             * @param[out] Face_number        Number of face in z-direction.
             *
             */
           uint
           give_face_z_of_position(
                   uint const & aLevel,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<uint> const & aIJKPosition) const;

           /**
            * Provides the faces, which are connected to an element.
            *
            * @param[in] aElement            Element number.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
            *
            * @param[out] Faces        Faces, which are connected to an element. In 2D: Faces(4,1) = [left face_X, right face_X, bottom face_Y, top_face_Y], In 2D: faces(6,1) = [left face_X, right face_X,  bottom face_Y, top_face_Y, front face_z, back face_z]
            *
            */
           Mat<uint>
           give_element_faces(
                   uint const & aElement,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
            * Provides the elements, which are connected to an edge.
            *
            * @param[in] aEdge            Edge Id.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
            *
            * @param[out] Element        Elements, which are connected to an edge.
            *
            */
           Mat<uint>
           give_elements_of_face(
                   uint const & aFaceId,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
            * Provides the position of an face
            *
            * @param[in] aFaceId            Face number.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
            *
            * @param[out] Position         Position I,J,K.
            *
            */
           Mat<uint>
           give_face_position(
                   uint const & aFaceId,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
            * Provides the ownership of an edge
            *
            * @param[in] aFaceId            Face number.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
            * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
            *
            * @param[out] Rank        Proc rank
            *
            */
           uint
           give_face_owner(
                   uint const & aFaceId,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<uint> const & aProcRange,
                   Mat<uint> const & aProcNeighbors) const;

           /**
            * Provides a list of procs which share the face
            *
            * @param[in] aFaceId            Face number.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
            * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
            *
            * @param[out] Share        Vector with procs, which share the face
            *
            */
           Mat<uint>
           give_face_share(
                   uint const & aFaceId,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<uint> const & aProcRange,
                   Mat<uint> const & aProcNeighbors) const;

    };
}

#endif /* SRC_MESH_CL_BASE_MESH_FACE_HPP_ */
