/*
 * cl_FEM_IQI_Volume_Fraction.hpp
 *
 *  Created on: Nov 20, 2019
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_FRACTION_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_FRACTION_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Volume_Fraction : public IQI
        {

            enum class IQI_Stabilization_Type
            {
                RECIPROCAL_TOTAL_VOLUME,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IQI_Stabilization_Type > mStabilizationMap;

            //------------------------------------------------------------------------------
        public:
            //------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IQI_Volume_Fraction();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Volume_Fraction(){};

            //------------------------------------------------------------------------------
            /**
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                              std::string                                aStabilizationString )
            {
                // FIXME check that stabilization string makes sense?
                std::cout<<static_cast< uint >( mStabilizationMap[ aStabilizationString ] )<<" string"<<std::endl;
                std::cout<<this->get_stabilization_parameters().size()<<" size"<<std::endl;

                // set the stabilization parameter in the stabilization parameter cell
                this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
            }

            //------------------------------------------------------------------------------

        private:

            /**
             * compute the quantity of interest
             * @param[ in ] aQI quantity of interest matrix to fill
             */
            void compute_QI( Matrix< DDRMat > & aQI );

            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] adQIdu derivative of quantity of interest matrix to fill
             */
            void compute_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu );

        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_ */
