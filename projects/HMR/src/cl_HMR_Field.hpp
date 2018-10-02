/*
 * cl_HMR_Field.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_

#include <memory>

#include "typedefs.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        class Field : public mtk::Field
        {
            //! pointer to database
            std::shared_ptr< Database > mDatabase;

            //! mesh that holds data
            Lagrange_Mesh_Base * mLagrangeMesh;

            // index of field in mesh
            const uint mFieldIndex;
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            Field(  const std::string             & aLabel,
                    std::shared_ptr< mtk::Mesh >    aMesh,
                    std::shared_ptr< Database >     aDatabase,
                    Lagrange_Mesh_Base *            aLagrangeMesh );

//------------------------------------------------------------------------------

            ~Field();

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_node_values();

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_node_values() const;

        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_ */
