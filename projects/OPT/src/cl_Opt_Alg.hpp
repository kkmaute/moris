#ifndef MORIS_OPTIMIZATION_CL_OPTALG_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALG_HPP_

// MORIS project header files.
#include "core.hpp"
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

            Matrix< DDRMat >  mAdvVec;    // Abstract Design Variable Vector
            Matrix< DDRMat >  mAdvUpVec;  // Upper bounds on Abstract Design Variable Vector
            Matrix< DDRMat >  mAdvLowVec; // Lower bounds on  Abstract Design Variable Vector

            real         mObjVal; // objective
            Matrix< DDRMat >  mConVal; // matrix of constraints
            Matrix< DDSMat >  mConType;// flags for type of constraint
            Matrix< DDRMat >  mDObj;   // derivative of objective w.r.t advs
            Matrix< DDRMat >  mDCon;   // derivative of constraints w.r.t advs
            Matrix< DDSMat >  mActive; // flag for active/inactive constraints

            ParameterList mParameterList; // The Algorithm specific parameter list

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
            Matrix< DDRMat >  & get_con()
            {
                return mConVal;
            }

            /**
             * @brief Get the gradient of the objective
             */
            Matrix< DDRMat >  & get_gradobj()
            {
                return mDObj;
            }

            /**
             * @brief Get the gradient of the constraints
             */
            Matrix< DDRMat >  & get_gradcon()
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
            template< typename Variant = ParameterListTypes >
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
            ParameterListTypes&
            set_param( char const* aKey)
            {
                return mParameterList(aKey);
            }
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALG_HPP_ */
