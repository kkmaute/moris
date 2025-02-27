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
            Module_Parameter_Lists&                     aParameterLists,
            const std::shared_ptr< moris::Library_IO >& aLibrary )
    {
        // Clear default vectors
        mProcessorDimensions.clear();
        mNumberOfElementsPerDimension.clear();
        mDomainDimensions.clear();
        mLagrangeOrders.clear();
        mLagrangePatterns.clear();
        mBSplineOrdersX.clear();
        mBSplineOrdersY.clear();
        mBSplineOrdersZ.clear();
        mBSplinePatterns.clear();

        // Get main HMR parameter list
        Parameter_List tHMRParameterList = aParameterLists( HMR::GENERAL )( 0 );

        // Number of elements per dimension
        mNumberOfElementsPerDimension = tHMRParameterList.get_vector< uint >( "number_of_elements_per_dimension" );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.size() == 2 or mNumberOfElementsPerDimension.size() == 3,
                "Number of elements must be a matrix of length 2 or 3." );

        // get processor decomposition method
        this->set_processor_decomp_method( tHMRParameterList.get< sint >( "processor_decomposition_method" ) );

        // get user defined processor dimensions. Only matters if decomp method == 3.
        mProcessorDimensions = tHMRParameterList.get_vector< uint >( "processor_dimensions" );

        // get domain dimensions
        mDomainDimensions = tHMRParameterList.get_vector< real >( "domain_dimensions" );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.size() == mDomainDimensions.size(),
                "length of domain_dimensions must be equal to number_of_elements_per_dimension." );

        // get domain offset
        mDomainOffset = tHMRParameterList.get_vector< real >( "domain_offset" );

        // Check size of offset
        if ( mDomainOffset.size() == 0 )
        {
            // If no offset specified, assume zero offset
            mDomainOffset.resize( mNumberOfElementsPerDimension.size(), 0.0 );
        }
        else
        {
            // Otherwise make sure size matches up
            MORIS_ERROR( mNumberOfElementsPerDimension.size() == mDomainOffset.size(),
                    "length of domain_offset must be equal to number_of_elements_per_dimension." );
        }

        // set buffer sizes
        this->set_refinement_buffer( tHMRParameterList.get< uint >( "refinement_buffer" ) );
        this->set_staircase_buffer( tHMRParameterList.get< uint >( "staircase_buffer" ) );
        this->set_union_pattern( tHMRParameterList.get< sint >( "union_pattern" ) );
        this->set_working_pattern( tHMRParameterList.get< sint >( "working_pattern" ) );

        string_to_vector_of_vectors( tHMRParameterList.get< std::string >( "lagrange_to_bspline" ), mLagrangeToBSplineMesh );

        if ( tHMRParameterList.get< sint >( "severity_level" ) != 1 )
        {
            this->set_severity_level( tHMRParameterList.get< sint >( "severity_level" ) );
        }

        // set truncation flag
        this->set_bspline_truncation( tHMRParameterList.get< bool >( "truncate_bsplines" ) );

        // set minimum initial refinement
        mInitialRefinementLevel = tHMRParameterList.get_vector< uint >( "initial_refinement" );

        // get multigrid parameter
        this->set_multigrid( tHMRParameterList.get< bool >( "use_multigrid" ) );

        // get renumber lagrange nodes
        this->set_renumber_lagrange_nodes( tHMRParameterList.get< sint >( "renumber_lagrange_nodes" ) == 1 );

        // get multigrid parameter
        this->set_number_aura( tHMRParameterList.get< bool >( "use_number_aura" ) );

        this->set_use_advanced_t_matrices( tHMRParameterList.get< sint >( "use_advanced_T_matrix_scheme" ) == 1 );

        this->set_refinement_for_low_level_elements( tHMRParameterList.get< bool >( "use_refine_low_level_elements" ) );

        this->set_background_mesh_output_file_name( tHMRParameterList.get< std::string >( "write_background_mesh" ) );

        this->set_lagrange_mesh_output_file_name( tHMRParameterList.get< std::string >( "lagrange_mesh_output_file_name" ) );

        this->set_write_refinement_pattern_file_flag( tHMRParameterList.get< bool >( "write_refinement_pattern_file" ) );

        this->set_restart_refinement_pattern_file( tHMRParameterList.get< std::string >( "restart_refinement_pattern_file" ) );

        this->set_basis_fuction_vtk_file_name( tHMRParameterList.get< std::string >( "basis_function_vtk_file" ) );

        // Always create side sets when using parameter list
        mCreateSideSets = true;

        // get user-defined refinement functions
        Vector< std::string > tFunctionNames = string_to_vector< std::string >( tHMRParameterList.get< std::string >( "refinement_function_names" ) );

        MORIS_ERROR( ( aLibrary != nullptr ) or ( tFunctionNames.size() == 0 ),
                "User-defined refinement function names were provided without a library to load them from." );

        for ( uint tFunctionIndex = 0; tFunctionIndex < tFunctionNames.size(); tFunctionIndex++ )
        {
            Refinement_Function tRefineFunc = aLibrary->load_function< Refinement_Function >( tFunctionNames( tFunctionIndex ) );
            mRefinementFunctions.push_back( tRefineFunc );
        }

        // Read from Lagrange and B-spline mesh parameter lists
        if ( aParameterLists( HMR::LAGRANGE_MESHES ).size() > 0 )
        {
            // Get number of Lagrange meshes and set sizes
            uint tNumberOfLagrangeMeshes = aParameterLists( HMR::LAGRANGE_MESHES ).size();
            mLagrangePatterns.resize( tNumberOfLagrangeMeshes );
            mLagrangeOrders.resize( tNumberOfLagrangeMeshes );
            mLagrangeToBSplineMesh.resize( tNumberOfLagrangeMeshes );

            // Loop over Lagrange meshes
            for ( uint iLagrangeMeshIndex = 0; iLagrangeMeshIndex < tNumberOfLagrangeMeshes; iLagrangeMeshIndex++ )
            {
                // Pattern and order
                mLagrangePatterns( iLagrangeMeshIndex ) = aParameterLists( HMR::LAGRANGE_MESHES )( iLagrangeMeshIndex ).get< uint >( "pattern_index" );
                mLagrangeOrders( iLagrangeMeshIndex ) = aParameterLists( HMR::LAGRANGE_MESHES )( iLagrangeMeshIndex ).get< uint >( "order" );

                // Output mesh parameters
                if ( aParameterLists( HMR::LAGRANGE_MESHES )( iLagrangeMeshIndex ).get< bool >( "is_output_mesh" ) )
                {
                    mOutputMeshNames.push_back( aParameterLists( HMR::LAGRANGE_MESHES )( iLagrangeMeshIndex ).get< std::string >( "output_mesh_name" ) );
                }
            }

            // Get number of B-spline meshes and set sizes
            uint tNumberOfBSplineMeshes = aParameterLists( HMR::BSPLINE_MESHES ).size();
            mBSplinePatterns.resize( tNumberOfBSplineMeshes );
            mBSplineOrdersX.resize( tNumberOfBSplineMeshes, 0 );
            mBSplineOrdersY.resize( tNumberOfBSplineMeshes, 0 );
            mBSplineOrdersZ.resize( tNumberOfBSplineMeshes, 0 );

            // Loop over B-spline meshes
            for ( uint iBSplineMeshIndex = 0; iBSplineMeshIndex < tNumberOfBSplineMeshes; iBSplineMeshIndex++ )
            {
                // B-spline pattern
                mBSplinePatterns( iBSplineMeshIndex ) = aParameterLists( HMR::BSPLINE_MESHES )( iBSplineMeshIndex ).get< uint >( "pattern_index" );

                // Dimension for order checking
                uint tNumberOfDimensions = mDomainDimensions.size();

                // Get raw orders, resized to 3
                Vector< uint > tBSplineOrders = aParameterLists( HMR::BSPLINE_MESHES )( iBSplineMeshIndex ).get_vector< uint >( "orders" );
                MORIS_ASSERT( tBSplineOrders.size() > 0, "At least one order must be given for each HMR B-spline mesh." );

                // Assign order in x
                mBSplineOrdersX( iBSplineMeshIndex ) = tBSplineOrders( 0 );

                // Assign y and z orders depending on input
                if ( tBSplineOrders.size() == 1 )
                {
                    mBSplineOrdersY( iBSplineMeshIndex ) = tBSplineOrders( 0 );
                    mBSplineOrdersZ( iBSplineMeshIndex ) = tBSplineOrders( 0 ) * ( tNumberOfDimensions > 2 );
                }
                else
                {
                    MORIS_ASSERT( tBSplineOrders.size() == tNumberOfDimensions,
                            "If specifying more than one B-spline order, the number of orders must match the number of dimensions." );
                    mBSplineOrdersY( iBSplineMeshIndex ) = tBSplineOrders( 1 );
                    if ( tNumberOfDimensions > 2 )
                    {
                        mBSplineOrdersZ( iBSplineMeshIndex ) = tBSplineOrders( 2 );
                    }
                }

                // Pair with Lagrange mesh
                uint tLagrangeMeshIndex = aParameterLists( HMR::BSPLINE_MESHES )( iBSplineMeshIndex ).get< uint >( "paired_lagrange_mesh_index" );
                mLagrangeToBSplineMesh( tLagrangeMeshIndex ).push_back( iBSplineMeshIndex );
            }
        }
        else
        {
            // Warning - deprecated section ahead!
            MORIS_LOG_WARNING( "Only a single HMR parameter list was found. It is recommended to use parameter lists for each Lagrange and B-spline mesh." );

            // Set orders/patterns
            string_to_vector( tHMRParameterList.get< std::string >( "lagrange_orders" ), mLagrangeOrders );
            string_to_vector( tHMRParameterList.get< std::string >( "lagrange_pattern" ), mLagrangePatterns );
            string_to_vector( tHMRParameterList.get< std::string >( "bspline_pattern" ), mBSplinePatterns );

            // B-spline orders
            Vector< uint > tBSplineOrders;
            string_to_vector( tHMRParameterList.get< std::string >( "bspline_orders" ), tBSplineOrders );
            mBSplineOrdersX = tBSplineOrders;
            mBSplineOrdersY = tBSplineOrders;
            if ( mDomainDimensions.size() > 2 )
            {
                mBSplineOrdersZ = tBSplineOrders;
            }
            else
            {
                mBSplineOrdersZ.resize( tBSplineOrders.size(), 0 );
            }

            // Set output mesh information
            string_to_vector_of_vectors( tHMRParameterList.get< std::string >( "lagrange_output_meshes" ), mOutputMeshes );

            MORIS_ERROR( mOutputMeshes.size() <= 2,
                    "Output mesh list can only have one list of main output meshes and one list of secondary output meshes" );
        }

        // Output mesh size TODO look at this further, should be set as a default in parameter lists
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
        mBSplineOrdersX = aMeshOrders;
        mBSplineOrdersY = aMeshOrders;
        if ( mDomainDimensions.size() > 2 )
        {
            mBSplineOrdersZ = aMeshOrders;
        }
        else
        {
            mBSplineOrdersZ.resize( aMeshOrders.size(), 0 );
        }

        // make sure that max polynomial is up to date
        this->update_max_polynomial_and_truncated_buffer();
    }

    //--------------------------------------------------------------------------------

    uint Parameters::get_max_bspline_order() const
    {
        return std::max( std::max( mBSplineOrdersX.max(), mBSplineOrdersY.max() ), mBSplineOrdersZ.max() );
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::update_max_polynomial_and_truncated_buffer()
    {
        mMaxPolynomial = std::max( mLagrangeOrders.max(), get_max_bspline_order() );
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
        MORIS_ERROR( aPatterns.size() == mBSplineOrdersX.size(),
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

            MORIS_ERROR( mBSplinePatterns.size() == mBSplineOrdersX.size(),
                    "Number of B-spline meshes given by patterns (%lu) doesn't match number of B-spline meshes given by orders (%lu)",
                    mBSplinePatterns.size(),
                    mBSplineOrdersX.size() );

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
