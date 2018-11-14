/*
 * cl_MTK_Face_Cluster.hpp
 *
 *  Created on: Nov 5, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FACET_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_FACET_CLUSTER_HPP_

#include "typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp"
namespace moris
{
namespace mtk
{

class Facet_Cluster
{
public:
    /*
     * Default constructor
     */
    Facet_Cluster(){}

    //------------------------------------------------------------------------------
    /*
     * Constructor which takes the parent face and number of child faces
     */
    Facet_Cluster(moris_index aParentFaceIndex,
                 uint        aNumChildFaces):
                     mParentFaceIndex(aParentFaceIndex),
                     mChildFaceIndices(aNumChildFaces,1,MORIS_INDEX_MAX),
                     mParametricBounds(aNumChildFaces)
    {

    }
    //------------------------------------------------------------------------------
    /*
     * set child face with ordinal in cluster
     */
    void
    set_child_face(moris_index              aChildOrd,
                   moris_index              aChildFaceIndex,
                   Matrix< DDRMat > const & aParametricBounds)
    {
        MORIS_ASSERT(mChildFaceIndices(aChildOrd) == MORIS_INDEX_MAX,
                     "Child ordinal information already set in face cluster");

        // set the face index and parametric bounds
        mChildFaceIndices(aChildOrd) = aChildFaceIndex;
        mParametricBounds(aChildOrd) = aParametricBounds;
    }
    //------------------------------------------------------------------------------
    /*
     * Get child faces
     */
    Matrix<IndexMat> const &
    get_child_faces() const
    {
        return mChildFaceIndices;
    }
    //------------------------------------------------------------------------------
    /*
     * Get child faces
     */
    Matrix<DDRMat> const &
    get_child_parametric_bounds(moris_index aChildFaceOrd) const
    {
        return mParametricBounds(aChildFaceOrd);
    }
    //------------------------------------------------------------------------------


private:
    // The face which all parametric bounds are relative to.
    moris_index                   mParentFaceIndex;

    // Child face within a parent face
    Matrix<IndexMat>              mChildFaceIndices;

    // Bounding box of the child faces. Cell index corresponds to loc in mChildFaceIndex.
    moris::Cell<Matrix< DDRMat >> mParametricBounds;
};
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_FACET_CLUSTER_HPP_ */
