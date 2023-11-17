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

    // creates a parameter list with default inputs
    void
    load_hmr_parameter_list_from_xml( const std::string& aFilePath, ParameterList& aParameterList )
    {
        // create temporary Parser object
        XML_Parser          tParser( aFilePath );
        Cell< std::string > tFirst;
        Cell< std::string > tSecond;

        tParser.get_keys_from_subtree( "moris.hmr", "parameters", 0, tFirst, tSecond );

        for ( uint k = 0; k < tFirst.size(); ++k )
        {
            std::string tKey = tFirst( k );

            if ( tKey == "number_of_elements_per_dimension" )
            {
                aParameterList.set( "number_of_elements_per_dimension", tSecond( k ) );
            }
            else if ( tKey == "processor_decomposition_method" )
            {
                aParameterList.set( "processor_decomposition_method", (sint)std::stoi( tSecond( k ) ) );
            }
            if ( tKey == "processor_dimensions" )
            {
                aParameterList.set( "processor_dimensions", tSecond( k ) );
            }
            else if ( tKey == "domain_dimensions" )
            {
                aParameterList.set( "domain_dimensions", tSecond( k ) );
            }
            else if ( tKey == "domain_offset" )
            {
                aParameterList.set( "domain_offset", tSecond( k ) );
            }
            else if ( tKey == "domain_sidesets" )
            {
                aParameterList.set( "domain_sidesets", tSecond( k ) );
            }
            else if ( tKey == "refinement_buffer" )
            {
                aParameterList.set( "refinement_buffer", (sint)std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "staircase_buffer" )
            {
                aParameterList.set( "staircase_buffer", (sint)std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "bspline_orders" )
            {
                aParameterList.set( "bspline_orders", tSecond( k ) );
            }
            else if ( tKey == "lagrange_orders" )
            {
                aParameterList.set( "lagrange_orders", tSecond( k ) );
            }
            //            else if( tKey == "initial_refinement" )
            //            {
            //                aParameterList.set( "initial_refinement", ( sint ) std::stoi( tSecond( k ) ) );
            //            }
            else if ( tKey == "severity_level" )
            {
                aParameterList.set( "severity_level", (sint)std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "truncate_bsplines" )
            {
                aParameterList.set( "truncate_bsplines", (sint)string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "additional_lagrange_refinement" )
            {
                aParameterList.set( "additional_lagrange_refinement", (sint)std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "max_refinement_level" )
            {
                aParameterList.set( "max_refinement_level", (sint)std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "use_multigrid" )
            {
                aParameterList.set( "use_multigrid", (sint)string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "use_refinement_interrelation" )
            {
                aParameterList.set( "use_refinement_interrelation", (sint)string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "renumber_lagrange_nodes" )
            {
                aParameterList.set( "renumber_lagrange_nodes", (sint)string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "use_number_aura" )
            {
                aParameterList.set( "use_number_aura", (sint)string_to_bool( tSecond( k ) ) );
            }
        }
    }

    //--------------------------------------------------------------------------------

    /*
     * parameter list constructor
     */
    Parameters::Parameters(
            ParameterList&                       aParameterList,
            std::shared_ptr< moris::Library_IO > aLibrary )
    {
        string_to_mat( aParameterList.get< std::string >( "number_of_elements_per_dimension" ), mNumberOfElementsPerDimension );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == 2 || mNumberOfElementsPerDimension.length() == 3,
                "Number of elements must be a matrix of length 2 or 3." );

        // get processor decomposition method
        this->set_processor_decomp_method( aParameterList.get< sint >( "processor_decomposition_method" ) );

        // get user defined processor dimensions. Only matters if decomp method == 3.
        string_to_mat( aParameterList.get< std::string >( "processor_dimensions" ), mProcessorDimensions );

        // get domain dimensions
        string_to_mat( aParameterList.get< std::string >( "domain_dimensions" ), mDomainDimensions );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == mDomainDimensions.length(),
                "length of domain_dimensions must be equal to number_of_elements_per_dimension." );

        // get domain offset
        string_to_mat( aParameterList.get< std::string >( "domain_offset" ), mDomainOffset );

        // check sanity of input
        MORIS_ERROR( mNumberOfElementsPerDimension.length() == mDomainOffset.length(),
                "length of domain_offset must be equal to number_of_elements_per_dimension." );

        // set buffer sizes
        this->set_refinement_buffer( aParameterList.get< sint >( "refinement_buffer" ) );
        this->set_staircase_buffer( aParameterList.get< sint >( "staircase_buffer" ) );

        string_to_mat( aParameterList.get< std::string >( "domain_sidesets" ), mSideSets );

        string_to_cell_mat( aParameterList.get< std::string >( "lagrange_output_meshes" ), mOutputMeshes );

        MORIS_ERROR( mOutputMeshes( 0 ).numel() <= 1, "only one main output mesh allowed right now" );
        MORIS_ERROR( mOutputMeshes.size() <= 2,
                "Output mesh list can only have one list of main output meshes and one list of secondary output meshes" );

        if ( aParameterList.get< std::string >( "lagrange_output_mesh_names" ).empty() )
        {
            uint tOutputMeshSize = mOutputMeshes.size();

            if ( tOutputMeshSize == 1 )
            {
                mOutputMesheNames.resize( 1, "HMR_Mesh_Main" );
                mOutputNameToIndexMap[ "HMR_Mesh_Main" ] = mOutputMeshes( 0 )( 0 );
            }
            else if ( tOutputMeshSize == 2 )
            {
                uint tNumSecOutputMeshes = mOutputMeshes( 1 ).numel();
                mOutputMesheNames.resize( 1 + tNumSecOutputMeshes );
                mOutputNameToIndexMap[ "HMR_Mesh_Main" ] = mOutputMeshes( 0 )( 0 );

                mOutputMesheNames( 1 ) = "HMR_Mesh_Main";
                for ( uint Ik = 1; Ik < tNumSecOutputMeshes + 1; Ik++ )
                {
                    std::string tName              = "HMR_Sec_Mesh" + std::to_string( Ik );
                    mOutputMesheNames( Ik )        = tName;
                    mOutputNameToIndexMap[ tName ] = mOutputMeshes( 1 )( Ik - 1 );
                }
            }
        }
        else
        {
            string_to_cell( aParameterList.get< std::string >( "lagrange_output_mesh_names" ), mOutputMesheNames );

            uint tOutputMeshSize = mOutputMeshes.size();

            if ( tOutputMeshSize == 1 )
            {
                MORIS_ERROR( mOutputMesheNames.size() == 1,
                        "Number of output mesh names must be the same than number of output meshes" );

                mOutputNameToIndexMap[ mOutputMesheNames( 0 ) ] = mOutputMeshes( 0 )( 0 );
            }
            else if ( tOutputMeshSize == 2 )
            {
                uint tNumSecOutputMeshes = mOutputMeshes( 1 ).numel();
                MORIS_ERROR( mOutputMesheNames.size() == ( 1 + tNumSecOutputMeshes ),
                        "Number of output mesh names must be the same than number of output meshes" );

                mOutputNameToIndexMap[ mOutputMesheNames( 0 ) ] = mOutputMeshes( 0 )( 0 );

                for ( uint Ik = 0; Ik < mOutputMeshes( 1 ).numel(); Ik++ )
                {
                    mOutputNameToIndexMap[ mOutputMesheNames( Ik + 1 ) ] = mOutputMeshes( 1 )( Ik );
                }
            }
        }

        string_to_mat( aParameterList.get< std::string >( "lagrange_input_meshes" ), mLagrangeInputMeshes );
        string_to_mat( aParameterList.get< std::string >( "lagrange_orders" ), mLagrangeOrders );
        string_to_mat( aParameterList.get< std::string >( "lagrange_pattern" ), mLagrangePatterns );
        string_to_mat( aParameterList.get< std::string >( "bspline_orders" ), mBSplineOrders );
        string_to_mat( aParameterList.get< std::string >( "bspline_pattern" ), mBSplinePatterns );

        this->set_union_pattern( aParameterList.get< sint >( "union_pattern" ) );
        this->set_working_pattern( aParameterList.get< sint >( "working_pattern" ) );

        string_to_cell_mat( aParameterList.get< std::string >( "lagrange_to_bspline" ), mLagrangeToBSplineMesh );

        if ( aParameterList.get< sint >( "severity_level" ) != 1 )
        {
            this->set_severity_level( aParameterList.get< sint >( "severity_level" ) );
        }

        // set truncation flag
        this->set_bspline_truncation( (bool)aParameterList.get< sint >( "truncate_bsplines" ) );

        //        // set minimum initial refinement
        string_to_mat( aParameterList.get< std::string >( "initial_refinement" ), mInitialRefinementLevel );
        string_to_mat( aParameterList.get< std::string >( "initial_refinement_pattern" ), mInitialRefinementPattern );

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
        Cell< std::string > tFunctionNames = string_to_cell< std::string >( aParameterList.get< std::string >( "refinement_function_names" ) );

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

    // creates a parameter list from parameters
    ParameterList
    create_hmr_parameter_list( const Parameters* aParameters )
    {
        MORIS_ERROR( false, "create_hmr_parameter_list(), function not changed yet" );
        // create default values
        ParameterList tParameterList = prm::create_hmr_parameter_list();

        // buffer size
        tParameterList.set( "refinement_buffer", (sint)aParameters->get_refinement_buffer() );
        tParameterList.set( "staircase_buffer", (sint)aParameters->get_staircase_buffer() );

        // verbosity flag
        tParameterList.set( "severity_level", (sint)aParameters->get_severity_level() );

        // truncation flag
        tParameterList.set( "truncate_bsplines", (sint)aParameters->truncate_bsplines() );

        // initial refinement
        // tParameterList.set( "initial_refinement",     ( sint ) aParameters->get_initial_refinement() );

        tParameterList.set( "max_refinement_level", (sint)aParameters->get_max_refinement_level() );

        // side sets
        tParameterList.set( "domain_sidesets", aParameters->get_side_sets_as_string() );

        tParameterList.set( "use_multigrid", (sint)aParameters->use_multigrid() );

        tParameterList.set( "use_refinement_interrelation", (sint)aParameters->get_refinement_interrelation() );

        tParameterList.set( "renumber_lagrange_nodes", (sint)aParameters->get_renumber_lagrange_nodes() );

        tParameterList.set( "use_number_aura", (sint)aParameters->use_number_aura() );

        return tParameterList;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::copy_selected_parameters( const Parameters& aParameters )
    {
        MORIS_ERROR( false, "copy_selected_parameters(), function not changed yet" );
        // buffer size
        this->set_refinement_buffer( aParameters.get_refinement_buffer() );
        this->set_staircase_buffer( aParameters.get_staircase_buffer() );

        // verbosity flag
        this->set_severity_level( aParameters.get_severity_level() );

        // truncation flag
        this->set_bspline_truncation( aParameters.truncate_bsplines() );

        // initial refinement
        this->set_initial_refinement( aParameters.get_initial_refinement() );
        this->set_additional_lagrange_refinement( aParameters.get_additional_lagrange_refinement() );

        // side sets
        this->set_side_sets( aParameters.get_side_sets() );

        // gmsh scaling factor
        this->set_gmsh_scale( aParameters.get_gmsh_scale() );

        this->set_lagrange_orders( aParameters.get_lagrange_orders() );

        this->set_bspline_orders( aParameters.get_bspline_orders() );
        this->set_max_refinement_level( aParameters.get_max_refinement_level() );

        this->set_multigrid( aParameters.use_multigrid() );

        // set refinement interrelation parameter
        this->set_refinement_interrelation( aParameters.get_refinement_interrelation() );

        this->set_renumber_lagrange_nodes( aParameters.get_renumber_lagrange_nodes() );

        this->set_number_aura( aParameters.use_number_aura() );
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::copy_selected_parameters( ParameterList& aParameterList )
    {
        // create a temporary parameter object
        Parameters tParameters( aParameterList, nullptr );

        // copy values into myself
        this->copy_selected_parameters( tParameters );
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::error( const std::string& aMessage ) const
    {
        if ( par_rank() == 0 )
        {
            MORIS_ERROR( false, "%s", aMessage.c_str() );
        }
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::error_if_locked( const std::string& aFunctionName ) const
    {
        if ( mParametersAreLocked )
        {
            std::string tMessage = "Error: calling function Parameters->" + aFunctionName + "() is forbidden since parameters are locked.";

            this->error( tMessage );
        }
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
     * Sets processor dimensions.  Only matters for "user defined" decomp method
     * only, mDecompMethod==0.
     */
    void
    Parameters::set_processor_dimensions( const Matrix< DDUMat >& aProcessorDimensions )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_processor_dimensions" );

        mProcessorDimensions = aProcessorDimensions;
    }

    //--------------------------------------------------------------------------------

    /**
     * sets the mesh orders according to given matrix
     */
    void
    Parameters::set_lagrange_orders( const Matrix< DDUMat >& aMeshOrders )
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
    Parameters::set_bspline_orders( const Matrix< DDUMat >& aMeshOrders )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_lagrange_orders" );

        MORIS_ERROR( aMeshOrders.max() <= 3, "Polynomial degree must be between 1 and 3" );
        MORIS_ERROR( 1 <= aMeshOrders.min(), "Polynomial degree must be between 1 and 3" );

        mBSplineOrders = aMeshOrders;

        // make sure that max polynomial is up to date
        this->update_max_polynomial_and_truncated_buffer();
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
    Parameters::set_lagrange_patterns( const Matrix< DDUMat >& aPatterns )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_lagrange_patterns" );

        // test sanity of input
        MORIS_ERROR( aPatterns.length() == mLagrangeOrders.length(),
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
    Parameters::set_bspline_patterns( const Matrix< DDUMat >& aPatterns )
    {
        // test if calling this function is allowed
        this->error_if_locked( "set_bspline_patterns" );

        // test sanity of input
        MORIS_ERROR( aPatterns.length() == mBSplineOrders.length(),
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

            // get number of B-Spline meshes
            auto tNumberOfBSplineMeshes = mBSplineOrders.length();

            MORIS_ERROR( mBSplinePatterns.length() == tNumberOfBSplineMeshes,
                    "B-Spline pattern list does not match number of B-Splines" );

            // get number of Lagrange meshes
            auto tNumberOfLagrangeMeshes = mLagrangeOrders.length();

            MORIS_ERROR( mLagrangePatterns.length() == tNumberOfLagrangeMeshes,
                    "Lagrange pattern list does not match number of Lagrange meshes" );
        }
    }

    //--------------------------------------------------------------------------------

    std::string
    Parameters::get_side_sets_as_string() const
    {
        std::string aString;
        mat_to_string( mSideSets, aString );
        return aString;
    }

    //--------------------------------------------------------------------------------

    void
    Parameters::set_refinement_functions( Cell< Refinement_Function > aRefinementFunctions )
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
        const Cell< Matrix< DDUMat > >& tOutputMeshes = this->get_output_mesh();

        bool tIsOutputMesh = false;

        for ( uint k = 0; k < tOutputMeshes( 0 ).numel(); ++k )
        {
            if ( aMeshIndex == tOutputMeshes( 0 )( k ) )
            {
                tIsOutputMesh = true;
                break;
            }
        }

        return tIsOutputMesh;
    }

    //--------------------------------------------------------------------------------

}    // namespace moris::hmr
