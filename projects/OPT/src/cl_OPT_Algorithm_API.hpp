#ifndef MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_

// C++ header files
#include <string>

// MORIS project header files.
#include "typedefs.hpp" // COR/src
#include "cl_OPT_Algorithm.hpp" // OPT/src
#include "cl_OPT_Problem.hpp" // OPT/src
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Algorithm_API
        {
            private:

                Algorithm* mAlgorithm;

            public:

                /**
                 * Constructor
                 *
                 * @param[in] aAlgId Optimization Algorithm
                 */
                Algorithm_API(ParameterList aParameterList);

                Algorithm_API()
                {}

                /**
                 * Copy Constructor through cloning
                 */
                Algorithm_API(const Algorithm_API& aAlgAPI )
                    : mAlgorithm(nullptr )
                {
                    mAlgorithm = aAlgAPI.mAlgorithm->clone();
                }

                /**
                 * Destructor
                 */
                ~Algorithm_API();

                /**
                 * @brief Accessor to set a value in the parameter list of Algorithm
                 *
                 * @param[in] aKey Key corresponding to the mapped value that
                 *            needs to be accessed
                 */
                ParameterListTypes&
                set_param( const char* aKey )
                {
                    return mAlgorithm->set_param(aKey);
                }

                /**
                 * @brief Wrapper around Algorithm::solve
                 *
                 * @param[in] aOptProb Object of type Problem containing relevant
                 *            data regarding ADVs, the objective and constraints
                 */
                void solve(std::shared_ptr<Problem> aOptProb )
                {
                    mAlgorithm->solve(aOptProb);
                }
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_ */
