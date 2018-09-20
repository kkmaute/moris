#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Integrator::Integrator(
                const Integration_Rule & aIntegrationRule )
        {

            // check if this is a combined rule
            if ( aIntegrationRule.has_two_rules() )
            {
                // create space rule
                mSpaceCoeffs = aIntegrationRule.create_space_coeffs();

                // create time rule
                mTimeCoeffs  = aIntegrationRule.create_time_coeffs();

                // make sure that time rule is 1D
                MORIS_ERROR( mTimeCoeffs->get_number_of_dimensions() == 1,
                        "time rule must be one-dimensional" );

                // assign function pointer to number of dimensions
                mGetNumberOfDimensions
                    = &Integrator::get_number_of_dimensions_space_and_time;

                // assign function pointer for number of points
                mGetNumberOfPoints
                    = &Integrator::get_number_of_points_space_and_time;

                // assign function pointer for points
                mGetPoints
                    = &Integrator::get_points_space_and_time;

                // assign function pointer for weights
                mGetWeights
                    = &Integrator::get_weights_space_and_time;
            }
            else
            {
                // create spacetime rule
                mSpaceTimeCoeffs
                    = aIntegrationRule.create_space_time_coeffs();

                // assign function pointer to number of dimensions
                mGetNumberOfDimensions
                    = &Integrator::get_number_of_dimensions_spacetime;

                // assign function pointer for number of points
                mGetNumberOfPoints
                    = &Integrator::get_number_of_points_spacetime;

                // assign function pointer for points
                mGetPoints
                    = &Integrator::get_points_spacetime;

                // assign function pointer for weights
                mGetWeights
                    = &Integrator::get_weights_spacetime;
            }
        }

//------------------------------------------------------------------------------

        Integrator::~Integrator()
        {
            // delete space coeffs if they exist
            if ( mSpaceCoeffs != NULL )
            {
                delete mSpaceCoeffs;
            }

            // delete time coeffs if they exist
            if( mTimeCoeffs != NULL )
            {
                delete mTimeCoeffs;
            }

            // delete combinned coeffs if they exist
            if( mSpaceTimeCoeffs != NULL )
            {
                delete mSpaceTimeCoeffs;
            }
        }

//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_points()
        {
            return ( this->*( mGetNumberOfPoints ) )();
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_points()
        {
            return ( this->*( mGetPoints ) ) () ;
        }

//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_dimensions()
        {
            return ( this->*( mGetNumberOfDimensions ) ) () ;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_weights()
        {
            return ( this->*( mGetWeights ) ) () ;
        }

//------------------------------------------------------------------------------
//      private :
//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_dimensions_space_and_time()
        {
            return    this->mSpaceCoeffs->get_number_of_dimensions()
                    + this->mTimeCoeffs->get_number_of_dimensions();
        }

//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_points_space_and_time()
        {
            return    this->mSpaceCoeffs->get_number_of_points()
                    * this->mTimeCoeffs->get_number_of_points();
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_points_space_and_time()
        {
            // matrix with space points
            auto tSpacePoints  = this->mSpaceCoeffs->get_points();

            // matrix with time points
            auto tTimePoints = this->mTimeCoeffs->get_points();

            // get number of points in space
            uint tNumberOfSpacePoints
                =  this->mSpaceCoeffs->get_number_of_points();

            // get number of points in time
            uint tNumberOfTimePoints
                =  this->mTimeCoeffs->get_number_of_points();

            // get number of dimensions in space
            uint tNumberOfDimensions
                =  this->mSpaceCoeffs->get_number_of_dimensions();

            // create output matrix for spacetime
            Matrix< DDRMat > aPoints(
                    tNumberOfDimensions + 1,
                    tNumberOfSpacePoints*tNumberOfTimePoints );

            // initialize counter
            uint tCount = 0;

            // loop over time
            for( uint k=0; k<tNumberOfTimePoints; ++k )
            {
                // loop over space
                for( uint j=0; j<tNumberOfSpacePoints; ++j )
                {
                    // loop over dimensions in space
                    for( uint i=0; i<tNumberOfDimensions; ++i )
                    {
                        aPoints( i, tCount ) = tSpacePoints( i, j );
                    }

                    // time
                    aPoints( tNumberOfDimensions, tCount++ )
                    = tTimePoints( k );
                }
            }

            return aPoints;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_weights_space_and_time()
        {
            // matrix with space weights
            auto tSpaceWeights = this->mSpaceCoeffs->get_weights();

            // matrix with time weights
            auto tTimeWeights = this->mTimeCoeffs->get_weights();


            // get number of points in space
            uint tNumberOfSpacePoints
            =  this->mSpaceCoeffs->get_number_of_points();

            // get number of points in time
            uint tNumberOfTimePoints
            =  this->mTimeCoeffs->get_number_of_points();

            // create output matrix for weights
            Matrix< DDRMat > aWeights( 1,
                    tNumberOfSpacePoints*tNumberOfTimePoints );

            // initialize counter
            uint tCount = 0;

            // loop over time
            for( uint k=0; k<tNumberOfTimePoints; ++k )
            {
                // loop over space
                for( uint j=0; j<tNumberOfSpacePoints; ++j )
                {
                    // weight
                    aWeights( tCount++ )
                        = tSpaceWeights( j ) * tTimeWeights( k );
                }
            }

            return aWeights;
        }

//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_dimensions_spacetime()
        {
            return this->mSpaceTimeCoeffs->get_number_of_dimensions();
        }

//------------------------------------------------------------------------------

        uint
        Integrator::get_number_of_points_spacetime()
        {
            return this->mSpaceTimeCoeffs->get_number_of_points();
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_points_spacetime()
        {
            return this->mSpaceTimeCoeffs->get_points();
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat >
        Integrator::get_weights_spacetime()
        {
            return this->mSpaceTimeCoeffs->get_weights();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
