/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_GQI.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_GQI_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_GQI_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_QI.hpp"    //FEM/INT/src
#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class GQI : public QI
    {
      protected:
        //! initialization flag
        bool mIsInitialized = false;

        std::string                                             mGeometryName            = "";
        std::shared_ptr< const MSI::Design_Variable_Interface > mDesignVariableInterface = nullptr;

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        GQI();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~GQI() {};

        //------------------------------------------------------------------------------

        void set_geometry_name( const std::string& aGeometryName )
        {
            mGeometryName = aGeometryName;
        }

        //------------------------------------------------------------------------------

        void set_design_variable_interface( std::shared_ptr< const MSI::Design_Variable_Interface > aDesignVariableInterface )
        {
            mDesignVariableInterface = aDesignVariableInterface;
        }

        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void compute_QI() = 0;

        //------------------------------------------------------------------------------

        /**
         * BRENDAN documentation
         */
        virtual real get_QI() = 0;

      private:
        //------------------------------------------------------------------------------
        /**
         * initialize parameters
         */
        virtual void initialize();

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest and fill aQI with value
         * @param[ in ] aQI IQI value at evaluation point
         */
        virtual void compute_QI( Matrix< DDRMat >& aQI ) = 0;

        //------------------------------------------------------------------------------
        /**
         * evaluate the quantity of interest
         * @param[ in ] Matrix
         *\
         * @return void
         */
        virtual void evaluate_QI( Matrix< DDRMat >& aMat ) = 0;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void compute_dQIdu( real aWStar ) = 0;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        virtual void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) = 0;

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_GQI_HPP_ */
