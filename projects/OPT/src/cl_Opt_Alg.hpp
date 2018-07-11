#ifndef MORIS_OPTIMIZATION_CL_OPTALG_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALG_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_Param_List.hpp" // CON/src
#include "cl_Opt_Prob.hpp" // OPT/src

namespace moris
{
    namespace opt
    {
        class OptAlg
        {
        protected:

            // Declare member variables

            OptProb mOptProb; // Object of type optimization problem

            uint mNumAdv;   // number of advs
            uint mNumCon;   // number of constraints
            uint mNumEqCon; // number of equality constraints

            Mat< real > mAdvVec;    // Abstract Design Variable Vector
            Mat< real > mAdvUpVec;  // Upper bounds on Abstract Design Variable Vector
            Mat< real > mAdvLowVec; // Lower bounds on  Abstract Design Variable Vector

            real         mObjVal; // objective
            Mat < real > mConVal; // matrix of constraints
            Mat < sint > mConType;// flags for type of constraint
            Mat < real > mDObj;   // derivative of objective w.r.t advs
            Mat < real > mDCon;   // derivative of constraints w.r.t advs
            Mat < sint > mActive; // flag for active/inactive constraints

            Param_List< boost::variant< bool, sint, real, const char* > > mParameterList; // The Algorithm specific parameter list

        public:

            /**
             * Constructor
             */
            OptAlg();

            /**
             * @brief virtual copy constructor through cloning
             */
            virtual OptAlg*
            clone() const = 0;

            /**
             * Destructor
             */
            virtual ~OptAlg();

            /**
             * @brief Get the objective value
             */
            real & get_obj()
            {
                return mObjVal;
            }

            /**
             * @brief Get the vector of constraints
             */
            Mat< real > & get_con()
            {
                return mConVal;
            }

            /**
             * @brief Get the gradient of the objective
             */
            Mat< real > & get_gradobj()
            {
                return mDObj;
            }

            /**
             * @brief Get the gradient of the constraints
             */
            Mat< real > & get_gradcon()
            {
                return mDCon;
            }

            /**
             * @brief Calls the derived optimization algorithm
             *
             * @param[in] aOptProb Object of type OptProb containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            virtual void solve( OptProb & aOptProb ) = 0;

            /**
             * @brief Initialize the member variables
             */
            void initialize();

            /**
             * @brief Evaluates the objectives and the constraints
             *
             * @param[in] aIter Optimization iteration number
             */
            void func( uint aIter );

            /**
             *@brief Evaluates the derivative of objectives and the constraints
             */
            void grad();

            /**
             *@brief Orders the constraints such that the equality constraints
             *       appear before the inequality constraints
             */
            void order_con();

            /**
             *@brief Orders the gradients of the constraints in an order similar
             *       to the order of constraints
             */
            void order_grad_con();

            /**
             * @brief Accessor for the parameter list of OptAlg
             */
            template< typename Variant = boost::variant< bool, sint, real, const char* > >
            Param_List< Variant > &
            params()
            {
                return mParameterList;
            }

            /**
             * @brief Accessor to set a value in the parameter list of OptAlg
             *
             * @param[in] aKey Key corresponding to the mapped value that
             *            needs to be accessed
             */
            boost::variant< bool, sint, real, const char* >&
            set_param( char const* aKey)
            {
                return mParameterList(aKey);
            }
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALG_HPP_ */
