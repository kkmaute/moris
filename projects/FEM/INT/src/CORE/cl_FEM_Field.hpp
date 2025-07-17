/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
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
            Vector< mtk::Field_Type > mFieldType;

            //! Name of IQI which is used to populate field
            std::string mIQIName;

            //! input/output file name and path
            std::string mInputFilePath;
            std::string mOutputFilePath;
            uint mTimeIndex = 0;
            bool mUpdateFromFile = false;

            // ! bool indicating if field shall be populated with the help of an IQI
            bool mPopulateFieldWithIQI = false;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            Field(
                    const mtk::Mesh_Pair&  aMeshPair,
                    mtk::Field_Entity_Type tFieldEntityType         = mtk::Field_Entity_Type::NODAL,
                    uint                   aDiscretizationMeshIndex = 0 );

            //------------------------------------------------------------------------------

            ~Field() override;

            /**
             * Updates this field from file, if necessary.
             */
            void update_field() override;

            //------------------------------------------------------------------------------

            void set_field_type( const Vector< mtk::Field_Type >& aType );

            //-----------------------------------------------------------------------------

            void set_field_from_file(
                    const std::string& aString,
                    uint               aTimeIndex,
                    uint               aFieldIndex,
                    bool               aUpdateTimeIndex = false );

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
                    Vector< mtk::Field_Type > const & aFieldTypes );

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
            void
            compute_nodal_values() override
            {
                MORIS_ERROR( false, "fem::Field::compute_nodal_values - not implemented.\n" );
            }

            //-----------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    }    // namespace fem
} /* namespace moris */
