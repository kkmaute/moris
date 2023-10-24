/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_ADV_Manager.hpp
 *
 */

#pragma once

#include "cl_GEN_ADV.hpp"
#include "cl_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris::ge
{
    class ADV_Manager
    {

      private:
        Cell< ADV >       mADVs;
        Matrix< DDSMat >  mDeterminingADVIds;
        bool              mHasADVs;
        sol::Dist_Vector* mSharedADVs = nullptr;

      public:
        /**
         * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aParameters Additional parameters
         */
        ADV_Manager(
                Matrix< DDRMat >&       aADVs,
                const Matrix< DDUMat >& aFieldVariableIndices,
                const Matrix< DDUMat >& aADVIndices,
                const Matrix< DDRMat >& aConstants );

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        explicit ADV_Manager( const Matrix< DDRMat >& aConstants );

        /**
         * Constructor, sets variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aFieldVariableIndices Variable indices for assigning the shared ADV IDs
         * @param aSharedADVIds Shared ADV IDs needed
         */
        ADV_Manager(
                const Matrix< DDUMat >& aFieldVariableIndices,
                const Matrix< DDSMat >& aSharedADVIds );

        /**
         * Destructor
         */
        ~ADV_Manager();

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADVs
         */
        template< typename Vector_Type >
        void set_advs( Vector_Type& aADVs );
        
        /**
         * Gets the value of a specific design variable so it can be used as a part of a design discretization.
         * 
         * @param aVariableIndex Index of the variable in this manager to reference
         * @return Design variable value
         */
        real get_variable( uint aVariableIndex );

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs );

        /**
         * Gets the IDs of ADVs that this manager depends on for evaluations.
         *
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids();

        /**
         * Gets if this manager has ADVs (at least one non-constant parameter)
         *
         * @return if this manager has ADVs
         */
        bool has_advs();

      private:

        /**
         * Creates the ADVs managed by this object.
         *
         * @param aADVs ADV vector
         * @param aConstants Constants to fill in other values
         */
        void create_advs(
                Matrix< DDRMat >& aADVs,
                const Matrix< DDRMat >& aConstants );
    };
}
