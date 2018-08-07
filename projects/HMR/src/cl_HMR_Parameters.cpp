/*
 * cl_HMR_Parameters.cpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#include "assert.hpp"

#include "cl_HMR_Parameters.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        void
        Parameters::print() const
        {
            if ( mVerbose && par_rank() == 0 )
            {
                std::fprintf( stdout, "\n" );
                std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                std::fprintf( stdout, "  user defined settings\n" ) ;
                std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                std::fprintf( stdout, "\n" );
                if ( mNumberOfElementsPerDimension.length() == 1 )
                {
                    std::fprintf( stdout, "  elements per dimension ....... : %lu\n",
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 0 ) );
                }
                else if (  mNumberOfElementsPerDimension.length() == 2 )
                {
                    std::fprintf( stdout, "  elements per dimension ....... : %lu x %lu\n",
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 0 ),
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 1 ) );
                }
                else if (  mNumberOfElementsPerDimension.length() == 3 )
                {
                    std::fprintf( stdout, "  elements per dimension ....... : %lu x %lu x %lu\n",
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 0 ),
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 1 ),
                            ( long unsigned int ) mNumberOfElementsPerDimension ( 3 )
                    );
                }
                std::fprintf( stdout,     "  buffer size .................. : %lu\n", ( long unsigned int ) mBufferSize );
                std::fprintf( stdout,     "  max polynomial ............... : %lu\n", ( long unsigned int ) mMaxPolynomial );
                std::fprintf( stdout, "\n" );
                std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                std::fprintf( stdout, "  automatically defined settings\n" ) ;
                std::fprintf( stdout, "--------------------------------------------------------------------------------\n" ) ;
                std::fprintf( stdout, "\n" );
                std::fprintf( stdout,     "  dimension .................... : %u\n", ( unsigned int ) this->get_number_of_dimensions() );
                std::fprintf( stdout,     "  padding size ................. : %lu\n", ( long unsigned int ) this->get_padding_size() );
                std::fprintf( stdout, "\n");
            }
        }

//--------------------------------------------------------------------------------

        auto
        Parameters::get_padding_size() const -> decltype ( mBufferSize )
        {
            // returns the larger value of max polynomial abd buffer size.
            // in the future, filter with will be regarded here
            return ( mBufferSize > mMaxPolynomial )
                    ? ( mBufferSize ) : ( mMaxPolynomial );
        }

//--------------------------------------------------------------------------------

        void
        Parameters::set_number_of_elements_per_dimension(
                const Mat<luint> & aNumberOfElementsPerDimension )
        {
            if( ! mParametersAreLocked )
            {
                mNumberOfElementsPerDimension = aNumberOfElementsPerDimension;
            }
            else
            {
                if( par_rank() == 0 )
                {
                    std::fprintf( stdout, "Error: calling Parameters->set_number_of_elements_per_dimension() is forbidden in this context." );
                    exit( -1 );
                }
            }
        }

//--------------------------------------------------------------------------------

        Mat<luint>
        Parameters::get_domain_ijk() const
        {
            // ask settings for number of dimensions
            auto tNumberOfDimensions = get_number_of_dimensions();

            // calculate padding size
            auto tPaddingSize = get_padding_size();

            // allocate output matrix
            Mat<luint> aDomain( 2, tNumberOfDimensions );

            // write beginning and ending of ijk domain in output matrix
            for ( uint k=0; k<tNumberOfDimensions; ++k )
            {
                aDomain( 0, k ) = tPaddingSize;
                aDomain( 1, k ) = aDomain( 0, k )
                        + mNumberOfElementsPerDimension( k ) - 1;
            }

            return aDomain;
        }

//--------------------------------------------------------------------------------

        /**
         * returns with, height and length of specified domain
         *
         * @return Mat<real>
         */
        Mat< real >
        Parameters::get_domain_dimensions() const
        {
            // see if dimensions have been set
            if( mDomainDimensions.length() != 0 )
            {
                // return user defined dimensions
                return mDomainDimensions;
            }
            else
            {
                // use default setting:

                // dimensions
                uint tNumberOfDimensions
                = mNumberOfElementsPerDimension.length();

                // return defalult values
                Mat< real > aDimensions( tNumberOfDimensions, 1 );

                // loop over all dimensions
                for( uint k=0; k<tNumberOfDimensions; ++k )
                {
                    // cast element number to real
                    aDimensions(k) = ( real ) mNumberOfElementsPerDimension( k );
                }

                // return domain so taht element length equals unity
                return aDimensions;
            }
        }

//--------------------------------------------------------------------------------

        /**
         * sets the maximum polynomial degree to given value
         *
         * @param[in] aMaxPolynomial
         *
         * @return void
         */
        void
        Parameters::set_max_polynomial( const luint & aMaxPolynomial )
        {
            if( ! mParametersAreLocked )
            {
                mMaxPolynomial = aMaxPolynomial;
            }
            else
            {
                if( par_rank() == 0 )
                {
                    std::fprintf( stdout, "Error: calling Parameters->set_max_polynomial() is forbidden in this context." );
                    exit( -1 );
                }
            }
        }

//--------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
