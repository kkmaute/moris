/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field.hpp
 *
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
          private:
            //! Field type. Identifier for a field interpolator. Similar to dof type or dv type
            Cell< enum mtk::Field_Type > mFieldType;

            //! Name of IQI which is used to populate field
            std::string mIQIName;

            //! output file name and path
            std::string mOutputFilePath;

            // ! bool indicating if field shall be populated with the help of an IQI
            bool mPopulateFieldWithIQI = false;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            Field(
                    mtk::Mesh_Pair              aMeshPair,
                    enum mtk::Field_Entity_Type tFieldEntityType         = mtk::Field_Entity_Type::NODAL,
                    uint                        aDiscretizationMeshIndex = 0 );

            //------------------------------------------------------------------------------

            virtual ~Field();

            //------------------------------------------------------------------------------

            void set_field_type( const moris::Cell< mtk::Field_Type >& aType );

            //-----------------------------------------------------------------------------

            void set_field_from_file( const std::string& aString );

            //-----------------------------------------------------------------------------

            void set_field_to_file( const std::string& aString );

            //-----------------------------------------------------------------------------

            void output_field_to_file();

            //-----------------------------------------------------------------------------

            void set_IQI_name( const std::string& aString );

            //-----------------------------------------------------------------------------

            bool
            get_populate_field_with_IQI()
            {
                return mPopulateFieldWithIQI;
            }

            //-----------------------------------------------------------------------------

            const std::string& get_IQI_name();

            //-----------------------------------------------------------------------------

            void get_values(
                    Matrix< IndexMat > const &      aIndex,
                    Matrix< DDRMat >&               aValues,
                    Cell< mtk::Field_Type > const & aFieldTypes );

            //-----------------------------------------------------------------------------

            // FIXME replace this
            void
            set_field_value(
                    const moris_index              tIndex,
                    const moris::Matrix< DDRMat >& aValue )
            {
                mValues.set_row( tIndex, aValue );
            }

            //-----------------------------------------------------------------------------

            /**
             * @brief child class implementation: computes and stores nodal values
             */
            virtual void
            compute_nodal_values()
            {
                MORIS_ERROR( false, "fem::Field::compute_nodal_values - not implemented.\n" );
            }

            // ----------------------------------------------------------------------------------------------

            /**
             * @brief child class implementation: computes derivatives of nodal values
             */
            virtual void
            compute_derivatives_of_field_value(
                    Matrix< DDRMat >& aDerivatives,
                    Matrix< DDUMat >& aCoefIndices,
                    uint const &      aNodeIndex,
                    uint const &      aFieldIndex )
            {
                MORIS_ERROR( false, "fem::Field::compute_derivatives_of_field_value - not implemented.\n" );
            }

            //-----------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    }    // namespace fem
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_FEM_FIELD_HPP_ */
