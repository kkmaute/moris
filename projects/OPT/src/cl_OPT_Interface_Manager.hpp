//
// Created by christopherson on 3/20/20.
//

#ifndef MORIS_CL_OPT_INTERFACE_MANAGER_HPP
#define MORIS_CL_OPT_INTERFACE_MANAGER_HPP

#include "cl_OPT_Interface.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Interface_Manager : public Interface
        {
        private:
            Cell<std::shared_ptr<Interface>> mInterfaces;
            Matrix<DDUMat> mNumADVsPerInterface;
            Matrix<DDUMat> mNumCriteriaPerInterface;
            Matrix<DDSMat> mProcessorBoundaries;

            bool mSharedADVs;
            bool mParallel;

            uint mNumInterfaces;

        public:

            /**
             * Constructor
             */
            Interface_Manager(ParameterList aParameterList, Cell<std::shared_ptr<Interface>> aInterfaces);

            /**
             * Destructor
             */
            ~Interface_Manager()
            {
            }

            /**
             * Sets the individual interfaces based on a cell of parameter lists
             */
            void set_interfaces();

            /**
             * Initializes the vector of ADV values
             */
            Matrix<DDRMat> initialize_advs();

            /**
             * Gets the lower bound values for the advs
             *
             * @return vector of lower bounds
             */
            Matrix<DDRMat> get_lower_adv_bounds();

            /**
             * Gets the upper bound values for the advs
             *
             * @return vector of upper bounds
             */
            Matrix<DDRMat> get_upper_adv_bounds();

            /**
             * Gets the criteria values
             *
             * @return vector of criteria
             */
            Matrix<DDRMat> get_criteria(Matrix<DDRMat> aNewADVs);

            /**
             * Gets the derivative of the criteria with respect to the advs
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            Matrix<DDRMat> get_dcriteria_dadv();

            /**
             * Gets the local advs based on whether or not they are shared
             */
            Matrix<DDRMat> get_local_advs(Matrix<DDRMat> aGlobalADVs, uint tInterfaceIndex);

        };
    }   // namespace opt
}   // namespace moris

#endif //MORIS_CL_OPT_INTERFACE_MANAGER_HPP
