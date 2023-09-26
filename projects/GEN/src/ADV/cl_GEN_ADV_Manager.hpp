/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris::ge
{
    class ADV_Manager
    {

      protected:
        Cell< real* >    mVariables;
        Matrix< DDRMat > mSensitivities;

      private:
        Matrix< DDRMat >  mConstants;
        Matrix< DDSMat >  mDeterminingADVIds;
        bool              mHasADVs;
        sol::Dist_Vector* mSharedADVs = nullptr;

      protected:
        /**
         * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aParameters Additional parameters
         */
        template< typename Vector_Type >
        ADV_Manager(
                Vector_Type&            aADVs,
                const Matrix< DDUMat >& aFieldVariableIndices,
                const Matrix< DDUMat >& aADVIndices,
                const Matrix< DDRMat >& aConstants );

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        ADV_Manager( const Matrix< DDRMat >& aConstants );

        /**
         * Constructor, sets variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aFieldVariableIndices Variable indices for assigning the shared ADV IDs
         * @param aSharedADVIds Shared ADV IDs needed
         */
        ADV_Manager(
                const Matrix< DDUMat >& aFieldVariableIndices,
                const Matrix< DDSMat >& aSharedADVIds );

      public:
        /**
         * Destructor
         */
        virtual ~ADV_Manager();

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADVs
         */
        template< typename Vector_Type >
        void set_advs( Vector_Type& aADVs );

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        virtual void import_advs( sol::Dist_Vector* aOwnedADVs );

        /**
         * Gets the IDs of ADVs which this manager depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        virtual Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets if this manager has ADVs (non-constant parameters)
         *
         * @return if this manager has ADVs
         */
        bool has_advs();

      private:
        /**
         * Checks variable inputs and resizes the internal field variables based these inputs.
         *
         * @param aFieldVariableIndices Indices of Field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the Field variables
         */
        void assign_adv_dependencies(
                const Matrix< DDUMat >& aFieldVariableIndices,
                const Matrix< DDUMat >& aADVIndices );

        /**
         * Fills the remaining field variables with constant parameters.
         */
        void fill_constant_parameters();
    };
}
