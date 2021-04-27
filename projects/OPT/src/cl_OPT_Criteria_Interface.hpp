#ifndef MORIS_CL_OPT_CRITERIA_INTERFACE_HPP
#define MORIS_CL_OPT_CRITERIA_INTERFACE_HPP

#include "cl_Matrix.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Criteria_Interface
        {
            private:

                bool mEvaluated = false;
                Matrix<DDRMat> mSensitivities;

            protected:
                bool mInitializeOptimizationRestart = false;

            public:

                /**
                 * Constructor
                 */
                Criteria_Interface()
                {
                }

                /**
                 * Destructor
                 */
                virtual ~Criteria_Interface()
                {
                }

                /**
                 * Initializes the vectors of ADV values, lower bounds, and upper bounds
                 *
                 * @param aADVs Initial ADVs to be filled.
                 * @param aLowerBounds Lower ADV bounds to be filled.
                 * @param aUpperBounds Upper ADV bounds to be filled.
                 */
                virtual void initialize(
                        Matrix<DDRMat>& aADVs,
                        Matrix<DDRMat>& aLowerBounds,
                        Matrix<DDRMat>& aUpperBounds) = 0;

                /**
                 * Gets the criteria values given a new set of ADVs
                 *
                 * @return vector of criteria
                 */
                Matrix<DDRMat> get_criteria(const Matrix<DDRMat> & aNewADVs);

                /**
                 * Gets the derivatives of the criteria with respect to the advs, and computes if not already done
                 *
                 * @return matrix d(criteria)_i/d(adv)_j
                 */
                Matrix<DDRMat> get_dcriteria_dadv();

                /**
                 * indicates a request to restart the optimization
                 *
                 */
                bool get_restart_optimization()
                {
                    return mInitializeOptimizationRestart;
                }

            private:

                /**
                 * Gets the criteria values given a new set of ADVs
                 *
                 * @return vector of criteria
                 */
                virtual Matrix<DDRMat> perform(const Matrix<DDRMat> & aNewADVs) = 0;

                /**
                 * Computes the derivatives of the criteria with respect to the advs
                 *
                 * @return matrix d(criteria)_i/d(adv)_j
                 */
                virtual Matrix<DDRMat> compute_dcriteria_dadv() = 0;
        };
    }
}

#endif //MORIS_CL_OPT_CRITERIA_INTERFACE_HPP
