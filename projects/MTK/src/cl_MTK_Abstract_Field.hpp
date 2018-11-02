/*
 * cl_MTK_Abstract_Field.hpp
 *
 *  Created on: Oct 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_ABSTRACT_FIELD_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_ABSTRACT_FIELD_HPP_

#include <string>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Mesh;
//------------------------------------------------------------------------------
        /**
         * An abstract field is defined by a Mesh, an entity rank,
         * a label, and data
         */
        class Abstract_Field
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            std::string        mLabel;
            const EntityRank   mEntityRank;
            const Mesh       * mMesh;
            Matrix< DDRMat >   mData;
            const uint         mDimension;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            Abstract_Field(
                    const std::string & aLabel,
                    const EntityRank  & aRank,
                    const Mesh        * aMesh,
                    const uint          aDimension=1 );

//------------------------------------------------------------------------------

            /**
             * Altertative constructor that passes a shared pointer
             * The share pointer is unpacked to a raw pointer during construction
             */
            Abstract_Field(
                    const std::string             & aLabel,
                    const EntityRank              & aRank,
                    const std::shared_ptr< Mesh >   aMesh,
                    const uint                      aDimension=1 );

//------------------------------------------------------------------------------
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* PROJECTS_MTK_SRC_CL_MTK_ABSTRACT_FIELD_HPP_ */
