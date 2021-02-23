/*
 * cl_FEM_Field.hpp
 *
 *  Created on: Jan 19, 2021
 *      Author: schmidt
 */

#ifndef PROJECTS_HMR_SRC_CL_FEM_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_FEM_FIELD_HPP_

#include <memory>

#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Field.hpp"

namespace moris
{
    namespace mtk
    {
        class Field;
    }
    namespace fem
    {
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

        class Field : public mtk::Field
        {
            protected:


                //------------------------------------------------------------------------------
            public :
                //------------------------------------------------------------------------------

                Field( mtk::Mesh_Pair * aMeshPair,
                       uint const                         & aDiscretizationMeshIndex =0 );

                //------------------------------------------------------------------------------

                ~Field();

                //------------------------------------------------------------------------------

                void set_field_type( const uint & aType );

                //-----------------------------------------------------------------------------

                void set_field_from_file( const std::string & aString );

                //-----------------------------------------------------------------------------

                //-----------------------------------------------------------------------------

        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_FEM_FIELD_HPP_ */
