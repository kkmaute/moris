#ifndef MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_

// C++ header files
#include <string>

// MORIS project header files.
#include "typedefs.hpp" // COR/src
#include "cl_Opt_Alg.hpp" // OPT/src

// Class forward declarations
class OptProb;

namespace moris
{
    namespace opt
    {
        class OptAlgAPI
        {
            private:

                OptAlg* mOptAlg;

            public:

                /**
                 * Constructor
                 *
                 * @param[in] aAlgId Optimization Algorithm
                 */
                OptAlgAPI( const std::string aAlgId );

                /**
                 * Copy Constructor through cloning
                 */
                OptAlgAPI( const OptAlgAPI& aAlgAPI )
                    : mOptAlg( nullptr )
                {
                    mOptAlg = aAlgAPI.mOptAlg->clone();
                }

                /**
                 * Destructor
                 */
                ~OptAlgAPI();

                /**
                 * @brief Accessor for the parameter list of OptAlg
                 */
                auto
                params()->decltype( mOptAlg->params() )
                {
                    return mOptAlg->params();
                }

                /**
                 * @brief Accessor to set a value in the parameter list of OptAlg
                 *
                 * @param[in] aKey Key corresponding to the mapped value that
                 *            needs to be accessed
                 */
                boost::variant< bool, sint, real, const char* >&
                set_param( const char* aKey )
                {
                    return mOptAlg->set_param( aKey );
                }

                /**
                 * @brief Wrapper around OptAlg::solve
                 *
                 * @param[in] aOptProb Object of type OptProb containing relevant
                 *            data regarding ADVs, the objective and constraints
                 */
                void solve( OptProb & aOptProb )
                {
                    mOptAlg->solve( aOptProb );
                }
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGAPI_HPP_ */
