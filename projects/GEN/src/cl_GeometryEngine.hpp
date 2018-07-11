/*
 * cl_geomeng_Geometry_Interpreter.hpp
 *
 *  Created on: Sep 12, 2016
 *      Author: doble
 */

#ifndef CL_GEOMETRY_INTERPRETER_HPP_
#define CL_GEOMETRY_INTERPRETER_HPP_
#include "cl_XTK_Enums.hpp" // XTK/src
#include "ios.hpp"
#include "core.hpp"
#include "chronos.hpp"
//#include "linalg.hpp"
#include "cl_Cell.hpp" // CON/src
#include "assert.hpp"

#include "tools.hpp"
#include "cl_Mesh_XTK.hpp" // XTK/src
#include "cl_Mesh.hpp" // MTK/src

/*
 * The geometry interpreter class establishes interface discretization
 * and controls evaluation of sensitivities. The class receives a mesh
 * and geometry description (Level Set Field) and outputs a surface description and
 * sensitivities of surface description with respect to geometric object.
 */
namespace moris
{
    class GeometryEngine
    {

    public:
        /*
         * Constructor
         */
        GeometryEngine();

        /*
         * Destructor
         */
        virtual ~GeometryEngine() = default;

        struct Object
        {
            Object();

            Object(moris::uint     aParentEntityIndex);

            ~Object();

            void
            set_parent_entity_index(moris::uint aEntityIndex);

            moris::uint
            get_parent_entity_index();

            /*
             * Set the coordinates where the geometry intersects the provided node connectivity
             *
             * @param[in] aInterfaceCoords - x y z coordinates of interface
             * @param[in] aRelevantNodes   - Nodes that the above coordinate is attached to
             */
            void
            set_interface_location(moris::Mat<moris::real>   aInterfaceCoords,
                                   moris::Mat<moris::uint>   aRelevantNodes);

            /* Currently set_interface_lcl_coord is only needed for an edge and requires only 1 value,
             * In future  want to extend to different dimension entities
             */
            void
            set_interface_lcl_coord(moris::real aLclCoord);
            /*
             * Global coordinate of interface point
             */
            void
            set_interface_glb_coord(moris::Mat<moris::real>  & aGlbCoord);
            /*
             * Sensitivity with respect to design relevant design variables
             */
            void
            set_dx_dp(moris::real & adxdp);
            moris::Mat<moris::real>
            get_interface_lcl_coord();

        private:
            moris::uint                           mParentEntityIndex;
            moris::real                           mdxdp;
            moris::Mat<moris::real>               mInterfaceGlbCoords;
            moris::Mat<moris::real>               mInterfaceLclCoords;
            moris::Cell<moris::Mat<moris::uint>>  mRelevantNodes;

        };
        /*is_intersected The XTK passes a mesh ptr to the Geometry Engine. The geometry engine then checks whether the object
         * is intersected using the intersection_check private function
         *
         * @param[in] aNodeCoords       - Node coordinate
         * @param[in] aNodeToEntityConn - Connectivity between nodes and parent entity
         * @param[in] aCheckType        - Specifies what type of intersection check is to be performed
         *                                   0 - No information on interface requred
         *                                   1 - information on interface required
         */
        moris::Cell<moris::GeometryEngine::Object>
        is_intersected(moris::Mat<moris::real>              const &   aNodeCoords,
                       moris::Mat<moris::uint>              const &   aNodetoEntityConn,
                       moris::uint                                    aCheckType);
//        moris::Cell<moris::geomeng::GeometryObject> const &
//        is_intersected(moris::Cell<moris::MeshXTK>               & aSimpleMesh,
//                       enum IntersectionCheck                      aCheckType);

    protected:

        virtual moris::Cell<moris::GeometryEngine::Object>
        intersection_check(moris::Mat<moris::real> const &   aNodeCoords,
                           moris::Mat<moris::uint> const &   aNodetoEntityConn,
                           moris::uint                       aCheckType) = 0;

    };

}
#endif /* CL_GEOMETRY_INTERPRETER_HPP_ */
