#ifndef PROJECTS_FEM_INT_SRC_IQI_CL_FEM_IQI_JUMP_DOF_HPP_
#define PROJECTS_FEM_INT_SRC_IQI_CL_FEM_IQI_JUMP_DOF_HPP_

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CON/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Jump_Dof : public IQI
        {
            //! initialization flag
            bool mIsInitialized = false;

            //! spatial and time derivative information
            uint mSpatialDerivativeDirection = 0;
            uint mSpatialDerivativeOrder     = 0;
            uint mTimeDerivativeOrder        = 0;
            //------------------------------------------------------------------------------

          public:

            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Jump_Dof();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Jump_Dof(){};

            //------------------------------------------------------------------------------

          private:

            //------------------------------------------------------------------------------
            
            /**
             * @brief initialize the parameters 
             * 
             */
            void initialize();

            //------------------------------------------------------------------------------
            /**
             * compute the quantity of interest
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_QI( Matrix< DDRMat >& aQI );

            //------------------------------------------------------------------------------
            /**
             * Evaluate the quantity of interest and fill aQI with value
             * @param[ in ] aQI IQI value at evaluation point
             */
            void compute_QI( real aWStar );

            //------------------------------------------------------------------------------

            /**
             * @brief evaluates IQI
             * 
             * @param aMat IQI matrix
             */

            void evaluate_QI( Matrix< DDRMat > & aMat );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dQIdu( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
             * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
             */
            void compute_dQIdu(
                    moris::Cell< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&             adQIdu );

            //------------------------------------------------------------------------------
        };
    } /* end namespace fem */
} /* end namespace moris */


#endif /* PROJECTS_FEM_INT_SRC_IQI_CL_FEM_IQI_Jump_Dof_HPP_ */
