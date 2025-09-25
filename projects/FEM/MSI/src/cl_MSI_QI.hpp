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
#include "cl_SOL_Dist_Vector.hpp"

namespace moris::MSI
{
    class QI
    {
      public:
        const Module_Type mModule;

      private:
        // Forward solve variables
        mutable real            mValue         = MORIS_REAL_MAX;    // Stores the value of the QI
        mutable bool            mIsEvaluated   = false;             // Flag if the QI value has been computed
        std::function< real() > mValueFunction = nullptr;           // Function to compute the QI value (may never be set if the QI is evaulated eagerly)

        // Sensitivity analysis variables
        mutable std::optional< sol::Dist_Vector* > mdQI;                      // Stores dQi/dPDV unless the module is a GQI. Then it is dQi/dADV
        std::function< sol::Dist_Vector*() >       mdQIFunction = nullptr;    // Function to compute the dQi/dADV sensitivity (may never be set if the sensitivity is evaluated eagerly)

      public:
        // Optimization problem constructor - eagerly computed QI and sensitivity
        explicit QI(
                Module_Type       aModule,
                const real        aValue = MORIS_REAL_MAX,
                sol::Dist_Vector* adQI   = nullptr );

        // Optimization problem constructor - lazily computed QI and sensitivity
        explicit QI(
                Module_Type                          aModule,
                std::function< real() >              aValueFunction,
                std::function< sol::Dist_Vector*() > aDQIdADVFunction );

        // Optimization problem constructor - eagerly computed QI and lazily computed sensitivity
        explicit QI(
                Module_Type                          aModule,
                const real                           aValue,
                std::function< sol::Dist_Vector*() > aDQIdADVFunction );

        // Optimization problem constructor - lazily computed QI and eagerly computed sensitivity
        explicit QI(
                Module_Type             aModule,
                std::function< real() > aValueFunction,
                sol::Dist_Vector*       adQIdADV = nullptr );

        /**
         * Gets the value of the QI, computes the value if it hasn't been computed already
         */
        real val() const;

        /*
         * Sets the value of the QI and marks it as evaluated
         */
        void set_val( real aValue );

        bool has_sensitivities() const;

        /**
         * Gets the dQi/PDV sensitivity (dQI/dADV for GEN QIs), computes the sensitivity if it hasn't been computed already
         */
        sol::Dist_Vector* sensitivity() const;

        /**
         * Sets the dQi/dADV sensitivity and marks it as evaluated
         */
        void set_sensitivity( sol::Dist_Vector* adQIdADV );

        void set_sensitivity( uint aIndex, real aValue );

        /**
         * Resets the QI to signify the value must be recomputed
         */
        void reset();
    };
}    // namespace moris::MSI
