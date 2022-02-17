/*
 * cl_MIG.hpp
 *
 *  Created on: Jan  31, 2022
 *      Author: momo
 */
#ifndef SRC_cl_MIG
#define SRC_cl_MIG

#include "cl_Param_List.hpp"

namespace moris::mtk
{
    class Mesh_Manager;
} // namespace moris::mtk

namespace moris::ge
{
    class Geometry_Engine;
}

namespace moris::mig
{
    class MIG
    {
      private:
        // mesh manager containing ip and ig meshes
        std::shared_ptr< moris::mtk::Mesh_Manager > mMeshManager;

        //parameter list
        moris::ParameterList                       mParameterList;

        //geomtry engine in order to link newly created nodes
        moris::ge::Geometry_Engine*                  mGeometryEngine;

      public:

        // ----------------------------------------------------------------------------

        /**
         * @brief Construct a new MIG object
         * 
         * @param aMeshManager 
         * @param aParameterList 
         * @param aGeometryEngine 
         */

        MIG( std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
            moris::ParameterList&                        aParameterList,
            moris::ge::Geometry_Engine*                  aGeometryEngine );

        // ----------------------------------------------------------------------------

        ~MIG()=default;

        // ----------------------------------------------------------------------------

        void
        perform();
    };
}// namespace moris::mig


#endif /* cl_MIG.hpp */