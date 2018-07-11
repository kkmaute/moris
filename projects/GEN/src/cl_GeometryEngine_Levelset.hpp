/*
 * cl_geomeng_Level_Set_Interpreter.hpp
 *
 *  Created on: Sep 12, 2016
 *      Author: Keenan Doble
 */

#ifndef CL_GEOMETRY_INTERPRETER_LEVELSET_HPP_
#define CL_GEOMETRY_INTERPRETER_LEVELSET_HPP_

#include "core.hpp"
#include "assert.hpp"
#include "tools.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_GeometryEngine.hpp"
#include "cl_Cell.hpp" // CON/src
#include "cl_Mesh_XTK.hpp" // XTK specific mesh header // XTK/src
#include "cl_Mesh.hpp"    // Mesh header // MTK/src
#include "cl_Interpolation.hpp" // TOL/src
/*
 * The Level_Set_Interpreter class is a derived class of the Geometry_Interpreter class
 * which performs the operations specific to a level set field geometry object. Currently,
 * it is assumed that level set field equation has the same coordinate basis as the mesh.
 *
 */

namespace moris
{

    class LevelsetEngine : public GeometryEngine
    {

    public:

        LevelsetEngine(moris::tools::Function*    aGeometryFunction);
        ~LevelsetEngine();

    private:
        moris::tools::Function* mLevelsetFunction; // Pointer to a parent function

        /**
         * @brief This function checks to see if the element has an interface and
         * depending on the aCheckType may create interface pts and parts.
         *
         * @param[in] aNodeCoords       - Node coordinate
         * @param[in] aNodeToEntityConn - Connectivity between nodes and parent entity
         * @param[in] aCheckType - Specifies what type of intersection check is to be performed
         *                         0 - No information on interface requred
         *                         1 - information on interface required
         */
        moris::Cell<moris::GeometryEngine::Object>
        intersection_check(moris::Mat<moris::real>              const &   aNodeCoords,
                           moris::Mat<moris::uint>              const &   aNodetoEntityConn,
                           moris::uint                                    aCheckType);

        /*
         * @brief compute_levelset_value for each node in the mesh.
         *
         */
        moris::Mat<moris::real>
        compute_levelset_value(moris::Mat<moris::real> const  & aNodeCoords);

        /*
         * computes the levelset value for an MeshXTK
         */
//        void
//        compute_levelset_value(moris::MeshXTK                   & aSimpleMesh,
//                               moris::Mat<moris::real>          & aLevelSetValues);

        /*
         * compute_intersection_info, calculates the relevant intersection information placed in the geometry object
         * @param[in]  aEntityNodeInds - node to entity connectivity
         * @param[in]  aNodeVars      - node level set values
         * @param[in]  aCheckType     - if a local location is necessary 1, else 0.
         * @param[out] Returns an intersection flag and local coordinates if aCheckType 1 in cell 1 and ndoe sensitivity information in cell 2 if intersection point located
         */
        moris::Cell<moris::Mat<moris::real>>
        compute_intersection_info(moris::Mat<moris::uint>  const  & aEntityNodeInds,
                                  moris::Mat<moris::real>  const  & aNodeVars,
                                  moris::Mat<moris::real>  const & aNodeCoords,
                                  moris::uint                      aCheckType);

        /*
         * REQUIREMENT each node connects to all other nodes (as is the case for tet4)
         * This function would need to be modified to support more complicated geometries
         */
        moris::Mat<moris::real>
        locate_interface(moris::Mat<moris::uint>  const & aElementNodeIds,
                         moris::Mat<moris::real>  const & aNodeVars);

    };
}

#endif /* CL_GEOMETRY_INTERPRETER_LEVELSET_HPP_ */
