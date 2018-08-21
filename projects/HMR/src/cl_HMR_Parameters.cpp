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
        Parameters::error( const std::string & aMessage ) const
        {
            if( par_rank() == 0 )
            {

                std::fprintf( stdout, aMessage.c_str() );
                exit( -1 );
            }
        }

//--------------------------------------------------------------------------------

        void
        Parameters::error_if_locked( const std::string & aFunctionName ) const
        {
            if( mParametersAreLocked )
            {
                std::string tMessage
                    = "Error: calling function Parameters->" + aFunctionName +
                    "() is forbidden since parameters are locked.";

                this->error( tMessage );
            }
        }

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

        /**
         * sets the mesh orders according to given matrix
         */
        void
        Parameters::set_lagrange_orders( const Mat< uint > & aMeshOrders )
        {
            // test if calling this function is allowed
            this->error_if_locked("set_lagrange_orders" );

            MORIS_ERROR( aMeshOrders.max() <= 3, "Polynomial degree must be between 1 and 3" );
            MORIS_ERROR( 1 <= aMeshOrders.min() , "Polynomial degree must be between 1 and 3" );

            mLagrangeOrders = aMeshOrders;

            // make sure that max polynomial is up to date
            this->update_max_polynomial_and_truncated_buffer();
        }

//--------------------------------------------------------------------------------

        /**
         * sets the mesh orders according to given matrix
         */
        void
        Parameters::set_bspline_orders( const Mat< uint > & aMeshOrders )
        {
            // test if calling this function is allowed
            this->error_if_locked("set_lagrange_orders" );

            MORIS_ERROR( aMeshOrders.max() <= 3, "Polynomial degree must be between 1 and 3" );
            MORIS_ERROR( 1 <= aMeshOrders.min() , "Polynomial degree must be between 1 and 3" );

            mBSplineOrders = aMeshOrders;

            // make sure that max polynomial is up to date
            this->update_max_polynomial_and_truncated_buffer();

        }

//--------------------------------------------------------------------------------

        void
        Parameters::update_max_polynomial_and_truncated_buffer()
        {
            mMaxPolynomial =
                    ( mLagrangeOrders.max() > mBSplineOrders.max() ) ?
                            (mLagrangeOrders.max() ) : ( mBSplineOrders.max() );

            if ( mBSplineTruncationFlag )
            {
                mBufferSize = mMaxPolynomial;
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
            // test if calling this function is allowed
            this->error_if_locked("set_number_of_elements_per_dimension");

            mNumberOfElementsPerDimension = aNumberOfElementsPerDimension;

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

//-------------------------------------------------------------------------------

        /**
         * sets the patterns for the Lagrange Meshes
         */
        void
        Parameters::set_lagrange_patterns( const Mat< uint > & aPatterns )
        {
            // test if calling this function is allowed
            this->error_if_locked("set_lagrange_patterns");

            // test sanity of input
            MORIS_ERROR( aPatterns.length() == mLagrangeOrders.length(),
                    "set_lagrange_patterns() : referred refinement pattern does not exist. Call set_lagrange_orders() first." );

            mLagrangePatterns = aPatterns;
        }

//-------------------------------------------------------------------------------

        /**
         * sets the patterns for the Lagrange Meshes
         */
        void
        Parameters::set_bspline_patterns( const Mat< uint > & aPatterns )
        {
            // test if calling this function is allowed
            this->error_if_locked("set_bspline_patterns");

            // test sanity of input
            MORIS_ERROR( aPatterns.length() == mBSplineOrders.length(),
                    "set_bspline_patterns() : referred refinement pattern does not exist. Call set_bspline_orders() first." );

            mBSplinePatterns = aPatterns;
        }

//--------------------------------------------------------------------------------

        /**
         * define which Lagrange mesh is linked to which B-Spline mesh
         */
        void
        Parameters::set_lagrange_to_bspline( const Mat< uint > & aBSplineMeshIndices )
        {
            // test if calling this function is allowed
            this->error_if_locked("set_lagrange_to_bspline");

            // test sanity of input
            MORIS_ERROR( aBSplineMeshIndices.max() < mBSplineOrders.length(),
                "set_lagrange_to_bspline() : referred refinement pattern does not exist. Call set_bspline_orders() first." );


            MORIS_ERROR( aBSplineMeshIndices.length() == mLagrangeOrders.length(),
                    "set_lagrange_to_bspline() : referred refinement pattern does not exist. Call set_lagrange_orders() first." );

            mLagrangeToBSpline = aBSplineMeshIndices;
        }

//--------------------------------------------------------------------------------

        void
        Parameters::set_default_patterns()
        {
            if( mLagrangePatterns.length() == 0 )
            {
                // ask parameters for number of meshes
                uint tNumberOfMeshes = this->get_number_of_lagrange_meshes();

                // set everthing to first pattern
                mLagrangePatterns.set_size( tNumberOfMeshes, 1, 0 );
            }

            if ( mBSplinePatterns.length() == 0 )
            {
                // ask parameters for number of meshes
                uint tNumberOfMeshes = this->get_number_of_bspline_meshes();

                mBSplinePatterns.set_size( tNumberOfMeshes, 1, 0 );
            }

            // set default links
            if ( mLagrangeToBSpline.length() == 0 )
            {
                uint tNumberOfMeshes = this->get_number_of_lagrange_meshes();

                if ( this->get_number_of_bspline_meshes() >= tNumberOfMeshes )
                {
                    mLagrangeToBSpline.set_size( tNumberOfMeshes, 1, 0 );

                    for( uint k=0; k<tNumberOfMeshes; ++k )
                    {
                        mLagrangeToBSpline( k ) = k;
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        /**
         * sets the maximum polynomial degree to given value
         *
         * @param[in] aMaxPolynomial
         *
         * @return void
         */
        /*void
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
        } */

//--------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
