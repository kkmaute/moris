/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Parameters.cpp
 *
 */

#include "cl_HMR_Parameters.hpp"    //HMR/src

#include "assert.hpp"
#include "fn_Parsing_Tools.hpp"

#include "fn_unique.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::hmr
{

    //--------------------------------------------------------------------------------

    /*
     * parameter list constructor
     */
    Parameters::Parameters(
            Parameter_List&                             aParameterList,
            const std::shared_ptr< moris::Library_IO >& aLibrary )
    {
        // Clear default vectors
        mProcessorDimensions.clear();
        mLagrangeOrders.clear();
        mLagrangePatterns.clear();
        mBSplineOrders.clear();
        mBSplinePatterns.clear();

        // Number of elements per dimension
        string_to_matrix( aParameterList.get< std::string >( "number_of_elements_per_dimension" ), mNumberOfElementsPerDimension );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == 2 || mNumberOfElementsPerDimension.length() == 3,
                "Number of elements must be a matrix of length 2 or 3." );

        // get processor decomposition method
        this->set_processor_decomp_method( aParameterList.get< sint >( "processor_decomposition_method" ) );

        // get user defined processor dimensions. Only matters if decomp method == 3.
        string_to_vector( aParameterList.get< std::string >( "processor_dimensions" ), mProcessorDimensions );

        // get domain dimensions
        string_to_matrix( aParameterList.get< std::string >( "domain_dimensions" ), mDomainDimensions );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == mDomainDimensions.length(),
                "length of domain_dimensions must be equal to number_of_elements_per_dimension." );

        // get domain offset
        string_to_matrix( aParameterList.get< std::string >( "domain_offset" ), mDomainOffset );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == mDomainOffset.length(),
                "length of domain_offset must be equal to number_of_elements_per_dimension." );

        // set buffer sizes
        this->set_refinement_buffer( aParameterList.get< sint >( "refinement_buffer" ) );
        this->set_staircase_buffer( aParameterList.get< sint >( "staircase_buffer" ) );

        string_to_matrix( aParameterList.get< std::string >( "domain_sidesets" ), mSideSets );

        string_to_vector_of_vectors( aParameterList.get< std::string >( "lagrange_output_meshes" ), mOutputMeshes );

        MORIS_ERROR( mOutputMeshes.size() <= 2,
                "Output mesh list can only have one list of main output meshes and one list of secondary output meshes" );

        if ( aParameterList.get< std::string >( "lagrange_output_mesh_names" ).empty() )
        {
            uint tOutputMeshSize = mOutputMeshes.size();

            if ( tOutputMeshSize == 1 )
            {
                mOutputMeshNames.resize( 1, "HMR_Mesh_Main" );
                mOutputNameToIndexMap[ "HMR_Mesh_Main" ] = mOutputMeshes( 0 )( 0 );
            }
            else if ( tOutputMeshSize == 2 )
            {
                uint tNumSecOutputMeshes = mOutputMeshes( 1 ).size();
                mOutputMeshNames.resize( 1 + tNumSecOutputMeshes );
                mOutputNameToIndexMap[ "HMR_Mesh_Main" ] = mOutputMeshes( 0 )( 0 );

                mOutputMeshNames( 1 ) = "HMR_Mesh_Main";
                for ( uint Ik = 1; Ik < tNumSecOutputMeshes + 1; Ik++ )
                {
                    std::string tName              = "HMR_Sec_Mesh" + std::to_string( Ik );
                    mOutputMeshNames( Ik )        = tName;
                    mOutputNameToIndexMap[ tName ] = mOutputMeshes( 1 )( Ik - 1 );
                }
            }
        }
        else
        {
            string_to_vector( aParameterList.get< std::string >( "lagrange_output_mesh_names" ), mOutputMeshNames );

            uint tOutputMeshSize = mOutputMeshes.size();

            if ( tOutputMeshSize == 1 )
            {
                MORIS_ERROR( mOutputMeshNames.size() == 1,
                        "Number of output mesh names must be the same than number of output meshes" );

                mOutputNameToIndexMap[ mOutputMeshNames( 0 ) ] = mOutputMeshes( 0 )( 0 );
            }
            else if ( tOutputMeshSize == 2 )
            {
                uint tNumSecOutputMeshes = mOutputMeshes( 1 ).size();
                MORIS_ERROR( mOutputMeshNames.size() == ( 1 + tNumSecOutputMeshes ),
                        "Number of output mesh names must be the same than number of output meshes" );

                mOutputNameToIndexMap[ mOutputMeshNames( 0 ) ] = mOutputMeshes( 0 )( 0 );

                for ( uint Ik = 0; Ik < mOutputMeshes( 1 ).size(); Ik++ )
                {
                    mOutputNameToIndexMap[ mOutputMeshNames( Ik + 1 ) ] = mOutputMeshes( 1 )( Ik );
                }
            }
        }

        string_to_matrix( aParameterList.get< std::string >( "lagrange_input_meshes" ), mLagrangeInputMeshes );

        // Set orders/patterns
        string_to_vector( aParameterList.get< std::string >( "lagrange_orders" ), mLagrangeOrders );
        string_to_vector( aParameterList.get< std::string >( "lagrange_pattern" ), mLagrangePatterns );
        string_to_vector( aParameterList.get< std::string >( "bspline_orders" ), mBSplineOrders );
        string_to_vector( aParameterList.get< std::string >( "bspline_pattern" ), mBSplinePatterns );

        this->set_union_pattern( aParameterList.get< sint >( "union_pattern" ) );
        this->set_working_pattern( aParameterList.get< sint >( "working_pattern" ) );

        string_to_vector_of_vectors( aParameterList.get< std::string >( "lagrange_to_bspline" ), mLagrangeToBSplineMesh );

        if ( aParameterList.get< sint >( "severity_level" ) != 1 )
        {
            this->set_severity_level( aParameterList.get< sint >( "severity_level" ) );
        }

        // set truncation flag
        this->set_bspline_truncation( (bool)aParameterList.get< sint >( "truncate_bsplines" ) );

        //        // set minimum initial refinement
        string_to_matrix( aParameterList.get< std::string >( "initial_refinement" ), mInitialRefinementLevel );
        string_to_matrix( aParameterList.get< std::string >( "initial_refinement_pattern" ), mInitialRefinementPattern );

        MORIS_ERROR( mInitialRefinementLevel.numel() == mInitialRefinementPattern.numel(),
                "length of mInitialRefinementLevel must be equal to mInitialRefinementPattern." );

        this->set_max_refinement_level( aParameterList.get< sint >( "max_refinement_level" ) );

        // get multigrid parameter
        this->set_multigrid( aParameterList.get< sint >( "use_multigrid" ) == 1 );

        // get refinement interrelation parameter
        this->set_refinement_interrelation( aParameterList.get< sint >( "use_refinement_interrelation" ) == 1 );

        // get renumber lagrange nodes
        this->set_renumber_lagrange_nodes( aParameterList.get< sint >( "renumber_lagrange_nodes" ) == 1 );

        // get multigrid parameter
        this->set_number_aura( aParameterList.get< sint >( "use_number_aura" ) == 1 );

        this->set_use_advanced_t_matrices( aParameterList.get< sint >( "use_advanced_T_matrix_scheme" ) == 1 );

        this->set_refinement_for_low_level_elements( aParameterList.get< bool >( "use_refine_low_level_elements" ) );

        this->set_write_background_mesh( aParameterList.get< std::string >( "write_background_mesh" ) );

        this->set_write_output_lagrange_mesh( aParameterList.get< std::string >( "write_lagrange_output_mesh" ) );

        this->set_write_output_lagrange_mesh_to_exodus( aParameterList.get< std::string >( "write_lagrange_output_mesh_to_exodus" ) );

        this->set_write_refinement_pattern_file_flag( aParameterList.get< bool >( "write_refinement_pattern_file" ) );

        this->set_restart_refinement_pattern_file( aParameterList.get< std::string >( "restart_refinement_pattern_file" ) );

        this->set_basis_fuction_vtk_file_name( aParameterList.get< std::string >( "basis_function_vtk_file" ) );

        // get user-defined refinement functions
        Vector< std::string > tFunctionNames = string_to_vector< std::string >( aParameterList.get< std::string >( "refinement_function_names" ) );

        MORIS_ERROR( ( aLibrary != nullptr ) or ( tFunctionNames.size() == 0 ),
                "User-defined refinement function names were provided without a library to load them from." );

        for ( uint tFunctionIndex = 0; tFunctionIndex < tFunctionNames.size(); tFunctionIndex++ )
        {
            Refinement_Function tRefineFunc = aLibrary->load_function< Refinement_Function >( tFunctionNames( tFunctionIndex ) );
            mRefinementFunctions.push_back( tRefineFunc );
        }

        this->update_max_polynomial_and_truncated_buffer();
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::error_if_locked( const std::string& aFunctionName ) const
    {
        MORIS_ERROR( not mParametersAreLocked,
                "Error: calling function Parameters->%s() is forbidden since parameters are locked.", aFunctionName.c_str() );
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::lock()
    {
        mParametersAreLocked = true;
    }

    //--------------------------------------------------------------------------------
    void
    Parameters::print() const
    {
        if ( par_rank() == 0 )
        {
            MORIS_LOG_INFO( " " );
            MORIS_LOG_INFO( "-------------------------------------------------------------------------------- " );
            MORIS_LOG_INFO( "user defined settings " );
            MORIS_LOG_INFO( "-------------------------------------------------------------------------------- " );
            MORIS_LOG_INFO( " " );
            if ( mNumberOfElementsPerDimension.length() == 1 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %lu",
                        (long unsigned int)mNumberOfElementsPerDimension( 0 ) );
            }
            else if ( mNumberOfElementsPerDimension.length() == 2 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %lu x %lu ",
                        (long unsigned int)mNumberOfElementsPerDimension( 0 ),
                        (long unsigned int)mNumberOfElementsPerDimension( 1 ) );
            }
            else if ( mNumberOfElementsPerDimension.length() == 3 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %lu x %lu x %lu ",
                        (long unsigned int)mNumberOfElementsPerDimension( 0 ),
                        (long unsigned int)mNumberOfElementsPerDimension( 1 ),
                        (long unsigned int)mNumberOfElementsPerDimension( 2 ) );
            }

            MORIS_LOG_INFO( "refinement buffer............. : %lu", (long unsigned int)mRefinementBuffer );
            MORIS_LOG_INFO( "staircase buffer.............. : %lu", (long unsigned int)mStaircaseBuffer );
            MORIS_LOG_INFO( "max polynomial ............... : %lu", (long unsigned int)mMaxPolynomial );
            MORIS_LOG_INFO( " " );
            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
            MORIS_LOG_INFO( "automatically defined settings" );
            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
            MORIS_LOG_INFO( " " );
            MORIS_LOG_INFO( "dimension .................... : %u", (unsigned int)this->get_number_of_dimensions() );
            MORIS_LOG_INFO( "padding size ................. : %lu", (long unsigned int)this->get_padding_size() );
            MORIS_LOG_INFO( " " );
        }
    }

    //--------------------------------------------------------------------------------

    /**
     * Sets processor decomposition method.  Options 1, 2, or 3 as of 10/29/2020
     * 0=UserDefined. 1=Min Proc Interface (Original) 2=Min Mesh Interface
     */
    void
    Parameters::set_processor_decomp_method( uint aProcDecompMethod )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_processor_decomp_method" );

        mProcDecompMethod = aProcDecompMethod;
    }

    //--------------------------------------------------------------------------------

    /**
     * sets the mesh orders according to given matrix
     */
    void
    Parameters::set_lagrange_orders( const Vector< uint >& aMeshOrders )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_lagrange_orders" );

        MORIS_ERROR( aMeshOrders.max() <= 3, "Polynomial degree must be between 1 and 3" );
        MORIS_ERROR( 1 <= aMeshOrders.min(), "Polynomial degree must be between 1 and 3" );

        mLagrangeOrders = aMeshOrders;

        // make sure that max polynomial is up to date
        this->update_max_polynomial_and_truncated_buffer();
    }

    //--------------------------------------------------------------------------------

    /**
     * sets the mesh orders according to given matrix
     */
    void
    Parameters::set_bspline_orders( const Vector< uint >& aMeshOrders )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_lagrange_orders" );

        MORIS_ERROR( aMeshOrders.max() <= 3, "B-spline polynomial degree must be between 1 and 3" );
        MORIS_ERROR( 1 <= aMeshOrders.min(), "B-spline polynomial degree must be between 1 and 3" );

        // Assign B-spline orders
        mBSplineOrders = aMeshOrders;

        // make sure that max polynomial is up to date
        this->update_max_polynomial_and_truncated_buffer();
    }

    //--------------------------------------------------------------------------------

    uint Parameters::get_max_bspline_order() const
    {
        return mBSplineOrders.max();
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::update_max_polynomial_and_truncated_buffer()
    {
        mMaxPolynomial = ( mLagrangeOrders.max() > mBSplineOrders.max() ) ? ( mLagrangeOrders.max() ) : ( mBSplineOrders.max() );
    }

    //--------------------------------------------------------------------------------

    auto
    Parameters::get_padding_size() const -> decltype( mStaircaseBuffer )
    {
        // returns the larger value of max polynomial and buffer size.
        // in the future, filter with will be regarded here
        return std::max( std::max( mStaircaseBuffer, mMaxPolynomial ), mRefinementBuffer );
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_number_of_elements_per_dimension( const Matrix< DDLUMat >& aNumberOfElementsPerDimension )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_number_of_elements_per_dimension" );

        // check sanity of input
        MORIS_ERROR( aNumberOfElementsPerDimension.length() == 2 || aNumberOfElementsPerDimension.length() == 3,
                "Number of elements must be a matrix of length 2 or 3." );

        mNumberOfElementsPerDimension = aNumberOfElementsPerDimension;

        // auto setting for dimensions and offset
        this->set_default_dimensions_and_offset();
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_number_of_elements_per_dimension(
            luint aElementsX,
            luint aElementsY )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_number_of_elements_per_dimension" );

        mNumberOfElementsPerDimension.set_size( 2, 1 );
        mNumberOfElementsPerDimension( 0 ) = aElementsX;
        mNumberOfElementsPerDimension( 1 ) = aElementsY;

        // auto setting for dimensions and offset
        this->set_default_dimensions_and_offset();
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_number_of_elements_per_dimension(
            luint aElementsX,
            luint aElementsY,
            luint aElementsZ )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_number_of_elements_per_dimension" );

        mNumberOfElementsPerDimension.set_size( 3, 1 );
        mNumberOfElementsPerDimension( 0 ) = aElementsX;
        mNumberOfElementsPerDimension( 1 ) = aElementsY;
        mNumberOfElementsPerDimension( 2 ) = aElementsZ;

        // auto setting for dimensions and offset
        this->set_default_dimensions_and_offset();
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_default_dimensions_and_offset()
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_default_dimensions_and_offset" );

        auto tNumberOfDimensions = mNumberOfElementsPerDimension.length();

        // auto set for domain dimensions
        if ( mDomainDimensions.length() == 0 )
        {
            mDomainDimensions.set_size( tNumberOfDimensions, 1, 1.0 );
        }

        // auto set offset
        if ( mDomainOffset.length() == 0 )
        {
            mDomainOffset.set_size( tNumberOfDimensions, 1, 0.0 );
        }
    }
    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_dimensions( const Matrix< DDRMat >& aDomainDimensions )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( aDomainDimensions.length() == 2 || aDomainDimensions.length() == 3,
                "Domain Dimensions must be a matrix of length 2 or 3." );

        MORIS_ERROR( aDomainDimensions.max() > 0.0, "Domain Dimensions be greater than zero" );

        mDomainDimensions = aDomainDimensions;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_dimensions(
            real aDomainDimensionsX,
            real aDomainDimensionsY )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( aDomainDimensionsX > 0.0, "aDomainDimensionsX must be greater than zero" );

        MORIS_ERROR( aDomainDimensionsY > 0.0, "aDomainDimensionsY must be greater than zero" );

        mDomainDimensions.set_size( 2, 1 );
        mDomainDimensions( 0 ) = aDomainDimensionsX;
        mDomainDimensions( 1 ) = aDomainDimensionsY;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_dimensions(
            real aDomainDimensionsX,
            real aDomainDimensionsY,
            real aDomainDimensionsZ )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( aDomainDimensionsX > 0.0, "aDomainDimensionsX must be greater than zero" );

        MORIS_ERROR( aDomainDimensionsY > 0.0, "aDomainDimensionsY must be greater than zero" );

        MORIS_ERROR( aDomainDimensionsZ > 0.0, "aDomainDimensionsZ must be greater than zero" );

        mDomainDimensions.set_size( 3, 1 );
        mDomainDimensions( 0 ) = aDomainDimensionsX;
        mDomainDimensions( 1 ) = aDomainDimensionsY;
        mDomainDimensions( 2 ) = aDomainDimensionsZ;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_offset( const Matrix< DDRMat >& aDomainOffset )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_offset" );

        // check sanity of input
        MORIS_ERROR( aDomainOffset.length() == 2 || aDomainOffset.length() == 3,
                "Domain Offset must be a matrix of length 2 or 3." );

        mDomainOffset = aDomainOffset;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_offset(
            real aDomainOffsetX,
            real aDomainOffsetY )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_offset" );

        mDomainOffset.set_size( 2, 1 );
        mDomainOffset( 0 ) = aDomainOffsetX;
        mDomainOffset( 1 ) = aDomainOffsetY;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_offset(
            real aDomainOffsetX,
            real aDomainOffsetY,
            real aDomainOffsetZ )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_offset" );

        mDomainOffset.set_size( 3, 1 );
        mDomainOffset( 0 ) = aDomainOffsetX;
        mDomainOffset( 1 ) = aDomainOffsetY;
        mDomainOffset( 2 ) = aDomainOffsetZ;
    }

    //--------------------------------------------------------------------------------

    Matrix< DDLUMat >
    Parameters::get_domain_ijk() const
    {
        // ask settings for number of dimensions
        auto tNumberOfDimensions = get_number_of_dimensions();

        // calculate padding size
        auto tPaddingSize = get_padding_size();

        // allocate output matrix
        Matrix< DDLUMat > aDomain( 2, tNumberOfDimensions );

        // write beginning and ending of ijk domain in output matrix
        for ( uint k = 0; k < tNumberOfDimensions; ++k )
        {
            aDomain( 0, k ) = tPaddingSize;
            aDomain( 1, k ) = aDomain( 0, k ) + mNumberOfElementsPerDimension( k ) - 1;
        }

        return aDomain;
    }

    //--------------------------------------------------------------------------------

    /**
     * returns with, height and length of specified domain
     *
     * @return Matrix< DDRMat >
     */
    Matrix< DDRMat >
    Parameters::get_domain_dimensions() const
    {
        // see if dimensions have been set
        if ( mDomainDimensions.length() != 0 )
        {
            // return user defined dimensions
            return mDomainDimensions;
        }
        else
        {
            // use default setting:

            // dimensions
            uint tNumberOfDimensions = mNumberOfElementsPerDimension.length();

            // return defalult values
            Matrix< DDRMat > aDimensions( tNumberOfDimensions, 1 );

            // loop over all dimensions
            for ( uint k = 0; k < tNumberOfDimensions; ++k )
            {
                // cast element number to real
                aDimensions( k ) = (real)mNumberOfElementsPerDimension( k );
            }

            // return domain so that element length equals unity
            return aDimensions;
        }
    }

    //-------------------------------------------------------------------------------

    /**
     * sets the patterns for the Lagrange Meshes
     */
    void
    Parameters::set_lagrange_patterns( const Vector< uint >& aPatterns )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_lagrange_patterns" );

        // test sanity of input
        MORIS_ERROR( aPatterns.size() == mLagrangeOrders.size(),
                "set_lagrange_patterns() : referred refinement pattern does not exist. Call set_lagrange_orders() first." );

        MORIS_ERROR( aPatterns.max() < gNumberOfPatterns - 2,
                "set_lagrange_patterns() : Pattern %-5i and %-5i are reserved for union and working pattern. Choose different pattern or increase gNumberOfPatterns.",
                gNumberOfPatterns - 2,
                gNumberOfPatterns - 1 );

        mLagrangePatterns = aPatterns;
    }

    //-------------------------------------------------------------------------------

    /**
     * sets the patterns for the Lagrange Meshes
     */
    void
    Parameters::set_bspline_patterns( const Vector< uint >& aPatterns )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_bspline_patterns" );

        // test sanity of input
        MORIS_ERROR( aPatterns.size() == mBSplineOrders.size(),
                "set_bspline_patterns() : referred refinement pattern does not exist. Call set_bspline_orders() first." );

        MORIS_ERROR( aPatterns.max() < gNumberOfPatterns - 2,
                "set_bspline_patterns() : Pattern %-5i and %-5i are reserved for union and working pattern. Choose different pattern or increase gNumberOfPatterns.",
                gNumberOfPatterns - 2,
                gNumberOfPatterns - 1 );

        mBSplinePatterns = aPatterns;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::check_sanity() const
    {
        if ( par_rank() == 0 )
        {
            // get dimensions
            auto tNumberOfDimensions = this->get_number_of_dimensions();

            // check dimensions
            MORIS_ERROR( mNumberOfElementsPerDimension.length() == tNumberOfDimensions,
                    "Number of Elements Per Dimension does not match" );

            MORIS_ERROR( mDomainDimensions.length() == tNumberOfDimensions,
                    "Domain dimensions and Number of Elements per dimension do not match" );

            MORIS_ERROR( mDomainOffset.length() == tNumberOfDimensions,
                    "Domain offset and Number of Elements per dimension do not match" );

            MORIS_ERROR( mBSplinePatterns.size() == mBSplineOrders.size(),
                    "Number of B-spline meshes given by patterns (%lu) doesn't match number of B-spline meshes given by orders (%lu)",
                    mBSplinePatterns.size(),
                    mBSplineOrders.size() );

            MORIS_ERROR( mLagrangePatterns.size() ==  mLagrangeOrders.size(),
                    "Number of Lagrange meshes given by patterns (%lu) doesn't match number of Lagrange meshes given by orders (%lu)",
                    mLagrangePatterns.size(),
                    mLagrangeOrders.size() );
        }
    }

    //--------------------------------------------------------------------------------

    std::string
    Parameters::get_side_sets_as_string() const
    {
        std::string aString;
        matrix_to_string( mSideSets, aString );
        return aString;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_refinement_functions( const Vector< Refinement_Function >& aRefinementFunctions )
    {
        mRefinementFunctions = aRefinementFunctions;
    }

    //--------------------------------------------------------------------------------

    Refinement_Function
    Parameters::get_refinement_function( uint aFunctionIndex )
    {
        MORIS_ASSERT(
                aFunctionIndex < mRefinementFunctions.size(),
                "hmr::Parameters::get_refinement_function() - "
                "A user-defined refinement function with index %i was requested for use, "
                "but only %li user-defined refinement functions were provided to HMR.",
                aFunctionIndex,
                mRefinementFunctions.size() );

        return mRefinementFunctions( aFunctionIndex );
    }

    //--------------------------------------------------------------------------------

    bool
    Parameters::is_output_mesh( const uint aMeshIndex ) const
    {
        const Vector< Vector< uint > >& tOutputMeshes = this->get_output_mesh();
        bool tIsOutputMesh = false;

        // Loop over output meshes
        for ( const auto& iMeshList : tOutputMeshes )
        {
            for ( auto iOutputIndex : iMeshList )
            {
                if ( aMeshIndex == iOutputIndex )
                {
                    tIsOutputMesh = true;
                    break;
                }
            }
        }

        return tIsOutputMesh;
    }

    //--------------------------------------------------------------------------------

}    // namespace moris::hmr
