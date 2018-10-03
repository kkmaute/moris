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
            uint mFieldIndex;
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

            const std::string &
            get_label() const;

//------------------------------------------------------------------------------

            void
            set_label( const std::string & aLabel );

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_node_values();

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_node_values() const;

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_coefficients();
//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_coefficients() const;

//------------------------------------------------------------------------------

            /**
             * sets the pointer of the mesh to another mesh
             * this is needed by the mapper
             */
            void
            change_mesh( Lagrange_Mesh_Base * aMesh, const uint aFieldIndex );

//------------------------------------------------------------------------------

            /**
             * returns the pointer of the underlying mesh
             */
            const Lagrange_Mesh_Base *
            get_mesh() const
            {
                return mLagrangeMesh;
            }

//------------------------------------------------------------------------------

            /**
             * return the field index on the linked mesh
             */
            uint
            get_field_index() const
            {
                return mFieldIndex;
            }

//------------------------------------------------------------------------------

            void
            save_field_to_hdf5( const std::string & aFilePath );

//------------------------------------------------------------------------------

            void
            save_bspline_coeffs_to_binary( const std::string & aFilePath );

//------------------------------------------------------------------------------

            void
            save_node_values_to_binary( const std::string & aFilePath );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_CL_HMR_FIELD_HPP_ */
