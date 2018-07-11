#ifndef MORIS_OPTIMIZATION_CL_OPTALGSQP_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALGSQP_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Opt_Alg.hpp" // Base class // OPT/src

//----------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        class OptAlgSQP : public OptAlg
        {
        private:

            uint mOptIter; // optimization iteration counter

            const int    MAXINT    = std::numeric_limits<int>::max();     // (int) (std::pow( 2, 63 )) - 1;
            const double MAXDOUBLE = std::numeric_limits<double>::max();  //        std::pow( 2, 63 )  - 1;

            // declare algorithm specific parameters
            int    mMinWLen;
            double mObjAdd;
            int    mObjRow;
            char*  mProb;

            /**
             * @brief External function call for computing objective, constraints, and
             *        their sensitivities, to interface with SNOPT library
             */
            friend void OptalgSQP_usrfun(
                    int* Status, int* n, double* x,
                    int* needf, int* nF, double* f,
                    int* needG, int* lenG, double* G,
                    char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru);

        public:

            /**
             * Constructor
             */
            OptAlgSQP();

            /**
             * Destructor
             */
            ~OptAlgSQP();

            /**
             * @brief copy constructor through cloning
             */
            OptAlg*
            clone() const
            {
                return new OptAlgSQP( *this );
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
             * @param[in] aOptProb Object of type OptProb containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            void solve( OptProb & aOptProb );

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
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALGSQP_HPP_ */
