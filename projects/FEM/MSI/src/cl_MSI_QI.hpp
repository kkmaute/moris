/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_QI.hpp
 *
 */

#include "moris_typedefs.hpp"
#include <functional>
#include <optional>
#include "cl_Matrix.hpp"
#include "cl_Library_Enums.hpp"

namespace moris::MSI
{
    class QI
    {

      public:
        // Public and const so it can be accessed directly without a getter
        const Module_Type mModule;

      private:
        // Forward solve variables
        mutable real            mValue         = MORIS_REAL_MAX;    // Stores the value of the QI
        mutable bool            mIsEvaluated   = false;             // Flag if the QI value has been computed
        std::function< real() > mValueFunction = nullptr;           // Function to compute the QI value (may never be set if the QI is evaulated eagerly)

        // Sensitivity analysis variables
        mutable std::optional< Matrix< DDRMat > > mdQIdADV;                      // Stores dQi/dADV sensitivity
        std::function< Matrix< DDRMat >() >       mDQIdADVFunction = nullptr;    // Function to compute the dQi/dADV sensitivity (may never be set if the sensitivity is evaluated eagerly)

      public:
        // // Forward problem only constructor - eagerly computed QI
        // explicit QI( const std::string& aName, const real& aValue );

        // // Forward problem only constructor - lazily computed QI
        // explicit QI( const std::string& aName, std::function< real() > aValueFunction );

        // Optimization problem constructor - eagerly computed QI and sensitivity
        explicit QI(
                Module_Type             aModule,
                const real&             aValue,
                const Matrix< DDRMat >& aDQIdADV = Matrix< DDRMat >() );

        // Optimization problem constructor - lazily computed QI and sensitivity
        explicit QI(
                Module_Type                          aModule,
                std::function< real() >              aValueFunction,
                std::function< Matrix< DDRMat >&() > aDQIdADVFunction );

        // Optimization problem constructor - eagerly computed QI and lazily computed sensitivity
        explicit QI(
                Module_Type                          aModule,
                const real&                          aValue,
                std::function< Matrix< DDRMat >&() > aDQIdADVFunction );

        // Optimization problem constructor - lazily computed QI and eagerly computed sensitivity
        explicit QI(
                Module_Type             aModule,
                std::function< real() > aValueFunction,
                const Matrix< DDRMat >& aDQIdADV = Matrix< DDRMat >() );

        /**
         * Gets the value of the QI, computes the value if it hasn't been computed already
         */
        real val() const;

        /*
         * Sets the value of the QI and marks it as evaluated
         */
        void set_val( real aValue );


        /**
         * Gets the dQi/dADV sensitivity, computes the sensitivity if it hasn't been computed already
         */
        const Matrix< DDRMat >& dADV() const;

        /**
         * Sets the dQi/dADV sensitivity and marks it as evaluated
         */
        void set_dADV( const Matrix< DDRMat >& aDQIdADV );

        /**
         * Resets the QI to signify the value must be recomputed
         */
        void reset();
    };
}    // namespace moris::MSI
