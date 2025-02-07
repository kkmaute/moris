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
        mNumberOfElementsPerDimension.clear();
        mDomainDimensions.clear();
        mLagrangeOrders.clear();
        mLagrangePatterns.clear();
        mBSplineOrders.clear();
        mBSplinePatterns.clear();

        // Number of elements per dimension
        mNumberOfElementsPerDimension = aParameterList.get< Vector< uint > >( "number_of_elements_per_dimension" );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.size() == 2 or mNumberOfElementsPerDimension.size() == 3,
                "Number of elements must be a matrix of length 2 or 3." );

        // get processor decomposition method
        this->set_processor_decomp_method( aParameterList.get< sint >( "processor_decomposition_method" ) );

        // get user defined processor dimensions. Only matters if decomp method == 3.
        mProcessorDimensions = aParameterList.get< Vector< uint > >( "processor_dimensions" );

        // get domain dimensions
        mDomainDimensions = aParameterList.get< Vector< real > >( "domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.size() == mDomainDimensions.size(),
                "length of domain_dimensions must be equal to number_of_elements_per_dimension." );

        // get domain offset
        string_to_vector( aParameterList.get< std::string >( "domain_offset" ), mDomainOffset );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.size() == mDomainOffset.size(),
                "length of domain_offset must be equal to number_of_elements_per_dimension." );

        // set buffer sizes
        this->set_refinement_buffer( aParameterList.get< uint >( "refinement_buffer" ) );
        this->set_staircase_buffer( aParameterList.get< uint >( "staircase_buffer" ) );

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
        string_to_vector( aParameterList.get< std::string >( "initial_refinement" ), mInitialRefinementLevel );

        this->set_max_refinement_level( aParameterList.get< sint >( "max_refinement_level" ) );

        // get multigrid parameter
        this->set_multigrid( aParameterList.get< sint >( "use_multigrid" ) == 1 );

        // get renumber lagrange nodes
        this->set_renumber_lagrange_nodes( aParameterList.get< sint >( "renumber_lagrange_nodes" ) == 1 );

        // get multigrid parameter
        this->set_number_aura( aParameterList.get< sint >( "use_number_aura" ) == 1 );

        this->set_use_advanced_t_matrices( aParameterList.get< sint >( "use_advanced_T_matrix_scheme" ) == 1 );

        this->set_refinement_for_low_level_elements( aParameterList.get< bool >( "use_refine_low_level_elements" ) );

        this->set_background_mesh_output_file_name( aParameterList.get< std::string >( "write_background_mesh" ) );

        this->set_lagrange_mesh_output_file_name( aParameterList.get< std::string >( "lagrange_mesh_output_file_name" ) );

        this->set_write_refinement_pattern_file_flag( aParameterList.get< bool >( "write_refinement_pattern_file" ) );

        this->set_restart_refinement_pattern_file( aParameterList.get< std::string >( "restart_refinement_pattern_file" ) );

        this->set_basis_fuction_vtk_file_name( aParameterList.get< std::string >( "basis_function_vtk_file" ) );

        // Always create side sets when using parameter list
        mCreateSideSets = true;

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
            if ( mNumberOfElementsPerDimension.size() == 1 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %u",
                        mNumberOfElementsPerDimension( 0 ) );
            }
            else if ( mNumberOfElementsPerDimension.size() == 2 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %u x %u ",
                        mNumberOfElementsPerDimension( 0 ),
                        mNumberOfElementsPerDimension( 1 ) );
            }
            else if ( mNumberOfElementsPerDimension.size() == 3 )
            {
                MORIS_LOG_INFO( "elements per dimension ....... : %u x %u x %u ",
                        mNumberOfElementsPerDimension( 0 ),
                        mNumberOfElementsPerDimension( 1 ),
                        mNumberOfElementsPerDimension( 2 ) );
            }

            MORIS_LOG_INFO( "refinement buffer............. : %u", mRefinementBuffer );
            MORIS_LOG_INFO( "staircase buffer.............. : %u", mStaircaseBuffer );
            MORIS_LOG_INFO( "max polynomial ............... : %u", mMaxPolynomial );
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
    Parameters::set_number_of_elements_per_dimension( const Vector< uint >& aNumberOfElementsPerDimension )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_number_of_elements_per_dimension" );

        // check sanity of input
        MORIS_ERROR( aNumberOfElementsPerDimension.size() == 2 or aNumberOfElementsPerDimension.size() == 3,
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

        mNumberOfElementsPerDimension.resize( 2 );
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

        mNumberOfElementsPerDimension.resize( 3 );
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

        auto tNumberOfDimensions = mNumberOfElementsPerDimension.size();

        // auto set for domain dimensions
        if ( mDomainDimensions.size() == 0 )
        {
            mDomainDimensions.resize( tNumberOfDimensions, 1.0 );
        }

        // auto set offset
        if ( mDomainOffset.size() == 0 )
        {
            mDomainOffset.resize( tNumberOfDimensions, 0.0 );
        }
    }
    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_dimensions( const Vector< real >& aDomainDimensions )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( aDomainDimensions.size() == 2 || aDomainDimensions.size() == 3,
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

        mDomainDimensions.resize( 2 );
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

        mDomainDimensions.resize( 3 );
        mDomainDimensions( 0 ) = aDomainDimensionsX;
        mDomainDimensions( 1 ) = aDomainDimensionsY;
        mDomainDimensions( 2 ) = aDomainDimensionsZ;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_domain_offset( const Vector< real >& aDomainOffset )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_domain_offset" );

        // check sanity of input
        MORIS_ERROR( aDomainOffset.size() == 2 || aDomainOffset.size() == 3,
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

        mDomainOffset.resize( 2 );
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

        mDomainOffset.resize( 3 );
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
    const Vector< real >&
    Parameters::get_domain_dimensions() const
    {
        return mDomainDimensions;
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
            MORIS_ERROR( mNumberOfElementsPerDimension.size() == tNumberOfDimensions,
                    "Number of Elements Per Dimension does not match" );

            MORIS_ERROR( mDomainDimensions.size() == tNumberOfDimensions,
                    "Domain dimensions and Number of Elements per dimension do not match" );

            MORIS_ERROR( mDomainOffset.size() == tNumberOfDimensions,
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
    Parameters::is_output_mesh( uint aMeshIndex ) const
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
