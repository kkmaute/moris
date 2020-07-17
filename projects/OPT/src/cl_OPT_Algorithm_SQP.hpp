#ifndef MORIS_CL_OPT_ALGORITHM_SQP_HPP_
#define MORIS_CL_OPT_ALGORITHM_SQP_HPP_

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"

namespace moris
{
    namespace opt
    {
        class Algorithm_SQP : public Algorithm
        {
        public:
            friend void sqp_user_function(
                    int* Status, int* n, double* x,
                    int* needf, int* nF, double* f,
                    int* needG, int* lenG, double* G,
                    char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru );

        private:

            uint mOptIter; // optimization iteration counter

            const int    MAXINT    = std::numeric_limits<int>::max();     // (int) (std::pow( 2, 63 )) - 1;
            const double MAXDOUBLE = std::numeric_limits<double>::max();  //        std::pow( 2, 63 )  - 1;

            // declare algorithm specific parameters
            int    mMinWLen;
            double mObjAdd;
            int    mObjRow;
            char*  mProb;

        public:

            /**
             * Constructor
             */
            Algorithm_SQP(ParameterList aParameterList);

            /**
             * Destructor
             */
            ~Algorithm_SQP();

            /**
             * @brief copy constructor through cloning
             */
            Algorithm*
            clone() const
            {
                return new Algorithm_SQP(*this );
            }

            /**
             * @brief Set the default parameter values
             */
            void set_default();

            /**
             * @brief Set the parameters for the algorithm
             */
            void set_params( char* cw, int lencw, int* iw, int leniw, double* rw, int lenrw );

            /**
             * @brief MORIS interface for solving of optimization problem using
             *        GCMMA
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            void solve(std::shared_ptr<Problem> aOptProb );

            /**
             * @brief MORIS-SQP interface for evaluation of objective,
             *        constraints and their derivative w.r.t. design variables
             *
             * @param[in] n Number of design variables
             * @param[in] x Array of design variables
             * @param[in] needG Flag from Snopt determining whether gradients
             *            need to be computed
             */
            void func_grad( int n, double* x, int needG );
        };

        /**
         * External function call for computing objective, constraints, and their sensitivities, to interface with
         * SNOPT library.
         */
        void sqp_user_function(
                int* Status, int* n, double* x,
                int* needf, int* nF, double* f,
                int* needG, int* lenG, double* G,
                char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru);
    }
}

#endif /* MORIS_CL_OPT_ALGORITHM_SQP_HPP_ */
