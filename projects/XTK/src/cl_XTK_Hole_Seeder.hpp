/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Hole_Seeder.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_HOLE_SEEDER_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_HOLE_SEEDER_HPP_

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"
namespace moris
{
    namespace ge
    {
        class Superellipsoid;
    }
}

namespace xtk
{

class Hole_Seeder
{
public:

    // Sphere constructors
    /*
     * constructor for if you have the mesh
     */
    Hole_Seeder( moris::mtk::Mesh* aMTKMesh,
                 moris::real       aRadiusX,
                 moris::real       aRadiusY,
                 moris::real       aRadiusZ,
                 moris::real       aNexp,
                 moris::uint       aNumInX,
                 moris::uint       aNumInY,
                 moris::uint       aNumInZ);

    /*
     * Constructor if the mesh hasnt been created
     */
    Hole_Seeder(moris::real       aRadiusX,
                moris::real       aRadiusY,
                moris::real       aRadiusZ,
                moris::real       aNexp,
                moris::uint       aNumInX,
                moris::uint       aNumInY,
                moris::uint       aNumInZ);

    // box constructors

    void
    set_mesh( moris::mtk::Mesh* aMTKMesh);

    void
    seed_field();

    moris::Matrix<moris::DDRMat> const &
    get_seeded_field();

    moris::Cell<std::shared_ptr<moris::ge::Superellipsoid>> &
    get_seeded_geometies();

private:
    moris::mtk::Mesh* mMTKMesh;
    moris::real mRadiusX;
    moris::real mRadiusY;
    moris::real mRadiusZ;
    moris::real mNexp;
    moris::uint mNumSpheresInX;
    moris::uint mNumSpheresInY;
    moris::uint mNumSpheresInZ;
    moris::Matrix<moris::DDRMat> mSeededField;
    moris::Cell<std::shared_ptr<moris::ge::Superellipsoid>> mSpheres;

};

}

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_HOLE_SEEDER_HPP_ */

