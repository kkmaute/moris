/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Parameters.hpp
 *
 */

#pragma once

#include <string>
#include <cstdio>

#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"
#include "assert.hpp"

#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_XML_Parser.hpp"
#include "cl_Library_IO.hpp"

namespace moris::mtk
{
    class Cell;
    class Field;
}    // namespace moris::mtk
namespace moris::hmr
{
    class Element;

    // User-defined refinement function
    typedef sint ( *Refinement_Function )(
            hmr::Element*           aElement,
            const Matrix< DDRMat >& aElementLocalValues );

    // User-defined refinement function
    typedef hmr::ElementalRefienmentIndicator ( *Refinement_Function_2 )(
            mtk::Cell*                    aElement,
            std::shared_ptr< mtk::Field > aField,
            uint                          tActivationPattern,
            uint&                         aMaxLevel );

    //--------------------------------------------------------------------------------

    /**
     * \brief This struct contains user defined settings
     */
    class Parameters
    {
      public:
        static const uint mLagrangeInputPattern  = 1;
        static const uint mLagrangeOutputPattern = 3;

      private:
        //! Processor decomposition method.  1=Original MPI decomp (min processor interface). 2=Min mesh interface. 0=Manually Defined
        uint mProcDecompMethod = 1;

        //! Processor layout if mProcDecompMethod is 0 (user defined here). Product MUST = # of processors used. Can be 1, 2, or 3 dimensions.
        Vector< uint > mProcessorDimensions = { 2, 2 };

        //! number of elements per direction in overall mesh, without aura
        //! 2D or 3D is determined by length of this vector
        Vector< uint > mNumberOfElementsPerDimension = { 2, 2 };

        //! width, height and depth of domain (without aura)
        Vector< real > mDomainDimensions;

        //! coordinate of first visible node
        Matrix< DDRMat > mDomainOffset;

        //! size of staircase buffer
        uint mRefinementBuffer = 1;
        uint mStaircaseBuffer  = 1;

        //! max polynomial to be supported
        uint mMaxPolynomial = 1;

        //! scale factor for gmsh output
        real mGmshScale = 1;

        //! flag telling if truncation is used
        bool mBSplineTruncationFlag = true;

        //! tells if critical features of the settings object are locked
        bool mParametersAreLocked = false;

        //! mesh orders, by default, a linear mesh is generated
        Vector< uint > mLagrangeOrders = { 1 };

        //! mesh orders, by default, a linear mesh is generated
        Vector< uint > mBSplineOrders = { 1 };

        //! defines which Lagrange mesh is associated with which refinement pattern
        Vector< uint > mLagrangePatterns = { 0 };

        //! defines which B-Spline mesh is associated with which refinement pattern
        Vector< uint > mBSplinePatterns = { 0 };

        //! defines which B-Spline mesh is associated with which lagrange mesh
        Vector< Vector< uint > > mLagrangeToBSplineMesh = { { 0 } };

        //! default union pattern
        uint mUnionPattern = gNumberOfPatterns - 1;

        //! default pattern for iterative refinement
        uint mWorkingPattern = gNumberOfPatterns - 2;

        //! Lagrange Meshes that are used for the output meshes
        Vector< Vector< uint > > mOutputMeshes    = { { 0 } };
        Vector< std::string >    mOutputMeshNames = { "" };

        moris::map< std::string, moris_index > mOutputNameToIndexMap;

        Matrix< DDUMat > mLagrangeInputMeshes = { {} };

        Matrix< DDUMat > mBSplineInputMeshes = { {} };

        Matrix< DDUMat > mInitialRefinementLevel   = { {} };
        Matrix< DDUMat > mInitialRefinementPattern = { {} };

        uint mAdditionalLagrangeRefinementLevel = 0;
        //! defines which SideSets are to be generated
        Matrix< DDUMat > mSideSets;

        bool mUseMultigrid = false;

        bool mNumberAura = false;

        bool mRefinementForLowLevelElements = false;

        bool mAdvancedTMatrices = false;

        std::string mBackgroundMeshFileName;
        std::string mLagrangeMeshFileName;

        std::string mBasisFunctionVtkFileName;

        bool mWriteRefinementPattern = false;

        std::string mRestartFromRefinedPatternFileName;

        //! maximum level for refinement. Default value is specified
        //! by global constant
        uint mMaxRefinementLevel = gMaxNumberOfLevels - 1;

        //! Renumber Lagrange Nodes
        bool mRenumberLagrangeNodes = false;

        // HMR user defined refinement function
        Vector< Refinement_Function > mRefinementFunctions;

      public:
        /**
         * returns user defined processor decomposition method
         *
         * @return uint
         */
        auto
        get_processor_decomp_method() const
                -> decltype( mProcDecompMethod )
        {
            return mProcDecompMethod;
        }

        /**
         * sets processor decomposition method
         *
         * @param[in] aProcDecompMethod uint
         *
         * @return void
         */
        void set_processor_decomp_method( uint aProcDecompMethod );

        /**
         * returns user defined processor decomposition method
         *
         * @return uint
         */
        auto
        get_processor_dimensions() const
                -> decltype( mProcessorDimensions )
        {
            return mProcessorDimensions;
        }

        /**
         * Constructor that loads parameters from a library
         */
        Parameters( Parameter_List&                         aParameterList,
                const std::shared_ptr< moris::Library_IO >& aLibrary );

        /**
         * Trivial constructor
         */
        Parameters() = default;

        /**
         * trivial destructor
         */
        ~Parameters() = default;

        /**
         * prints user settings passed to HMR
         *
         * @return void
         */
        void print() const;

        /**
         * Sets global logger severity level
         *
         * @param[in] aSeverityLevel Logger severity level
         */
        static void
        set_severity_level( sint aSeverityLevel )
        {
            gLogger.set_severity_level( aSeverityLevel );
        }

        /**
         * Gets the logger severity level
         * @return
         */
        static sint
        get_severity_level()
        {
            return gLogger.get_severity_level();
        }

        /**
         * sets the buffer size to given value
         *
         * @param[in] aBufferSize
         *
         * @return void
         */
        void
        set_refinement_buffer( uint aBufferSize )
        {
            mRefinementBuffer = aBufferSize;
        }

        /**
         * sets the buffer size to given value
         *
         * @param[in] aBufferSize
         *
         * @return void
         */
        void
        set_staircase_buffer( uint aBufferSize )
        {
            mStaircaseBuffer = aBufferSize;
        }

        /**
         * Gets the staircase buffer
         *
         * @return staircase buffer
         */
        auto
        get_staircase_buffer() const
                -> decltype( mStaircaseBuffer )
        {
            if ( mBSplineTruncationFlag )
            {
                return std::max( std::max( mStaircaseBuffer, mMaxPolynomial ), mRefinementBuffer );
            }
            else
            {
                return std::max( mStaircaseBuffer, mRefinementBuffer );
            }
        }

        /**
         * returns the buffer size
         *
         * @return luint
         */
        auto
        get_refinement_buffer() const
                -> decltype( mRefinementBuffer )
        {
            return mRefinementBuffer;
        }

        /**
         * sets the mesh orders according to given matrix
         *
         * @param aMeshOrders Lagrange mesh orders
         */
        void set_lagrange_orders( const Vector< uint >& aMeshOrders );

        /**
         * sets the mesh orders according to given matrix
         *
         * @param aMeshOrders B-spline mesh orders
         */
        void set_bspline_orders( const Vector< uint >& aMeshOrders );

        /**
         * returns a matrix with mesh orders
         */
        auto
        get_lagrange_orders() const
                -> decltype( mLagrangeOrders )
        {
            return mLagrangeOrders;
        }

        /**
         * sets the patterns for the Lagrange Meshes
         *
         * @param[ in ] aPattern patterns set by set_mesh_orders the Lagrange meshes refer to.
         *
         */
        void set_lagrange_patterns( const Vector< uint >& aPatterns );

        /**
         * Gets the patterns for the Lagrange meshes
         *
         * @return Lagrange mesh patterns
         */
        auto
        get_lagrange_patterns() const
                -> decltype( mLagrangePatterns )
        {
            return mLagrangePatterns;
        }

        /**
         * Returns the order of a Lagrange mesh by index
         *
         * @param aIndex Lagrange mesh index
         * @return Lagrange order
         */
        auto
        get_lagrange_order( uint aIndex ) const
                -> decltype( mLagrangeOrders( aIndex ) )
        {
            return mLagrangeOrders( aIndex );
        }

        /**
         * Returns the pattern of a Lagrange mesh by index.
         *
         * @param aIndex Lagrange mesh index
         * @return Lagrange pattern
         */
        auto
        get_lagrange_pattern( uint aIndex ) const
                -> decltype( mLagrangePatterns( aIndex ) )
        {
            return mLagrangePatterns( aIndex );
        }

        /**
         * sets the patterns for the B-Spline Meshes
         *
         * @param[ in ] aPattern patterns set by set_mesh_orders the B-Spline meshes refer to.
         */
        void set_bspline_patterns( const Vector< uint >& aPatterns );

        /**
         * Returns a list of all B-spline patterns.
         *
         * @return B-spline patterns
         */
        auto
        get_bspline_patterns() const
                -> decltype( mBSplinePatterns )
        {
            return mBSplinePatterns;
        }

        /**
         * Returns the pattern of a B-spline mesh by index.
         *
         * @param aIndex B-spline mesh index
         * @return B-spline pattern
         */
        auto
        get_bspline_pattern( uint aIndex ) const
                -> decltype( mBSplinePatterns( aIndex ) )
        {
            return mBSplinePatterns( aIndex );
        }

        /**
         * Returns the order of a B-spline mesh by index.
         *
         * @param aIndex B-spline mesh index
         * @return B-spline order
         */
        auto
        get_bspline_order( uint aIndex ) const
                -> decltype( mBSplineOrders( aIndex ) )
        {
            return mBSplineOrders( aIndex );
        }

        /**
         * Gets the maximum B-spline order of any mesh.
         *
         * @return Maximum B-spline order
         */
        uint get_max_bspline_order() const;

        /**
         * Sets information about which B-spline mesh is associated with which Lagrange mesh.
         *
         * @param Lagrange to B-spline mesh relationship
         */
        void
        set_lagrange_to_bspline_mesh( const Vector< Vector< uint > >& aLagrangeToBSplineMesh )
        {
            mLagrangeToBSplineMesh = aLagrangeToBSplineMesh;
        }

        /**
         * Gets the B-spline mesh indices that are associated with a given Lagrange mesh.
         *
         * @param aLagrangeMeshIndex Lagrange mesh index
         * @return B-spline mesh indices
         */
        Vector< uint >
        get_lagrange_to_bspline_mesh( uint aLagrangeMeshIndex ) const
        {
            return mLagrangeToBSplineMesh( aLagrangeMeshIndex );
        }

        /**
         * Gets the number of B-Spline meshes
         */
        uint
        get_number_of_bspline_meshes() const
        {
            return mBSplineOrders.size();
        }

        /**
         * Gets the total number of Lagrange meshes
         */
        uint
        get_number_of_lagrange_meshes() const
        {
            return mLagrangeOrders.size();
        }

        /**
         * Gets the index of the defined Lagrange output mesh for a specified order
         */
        const Vector< Vector< uint > >&
        get_output_mesh() const
        {
            return mOutputMeshes;
        }

        /**
         * set which lagrange meshes are used for an output
         */
        void
        set_output_meshes( const Vector< Vector< uint > >& aOutputMeshes )
        {
            // test if calling this function is allowed
            this->error_if_locked( "set_output_meshes" );

            mOutputMeshes = aOutputMeshes;
        };

        /**
         * returns the index of the defined Lagrange output mesh names for a specified order
         */
        const Vector< std::string >&
        get_output_mesh_names() const
        {
            return mOutputMeshNames;
        }

        /**
         * returns mesh index for name
         */
        moris_index
        get_mesh_index_by_name( const std::string& aName ) const
        {
            MORIS_ERROR( mOutputNameToIndexMap.key_exists( aName ),
                    "Parameters::get_mesh_index_by_name() Mesh name does not exist" );

            return mOutputNameToIndexMap.find( aName );
        }

        /**
         * returns true if mesh name exists
         */
        bool
        get_mesh_name_exists( const std::string& aName ) const
        {
            return mOutputNameToIndexMap.key_exists( aName );
        }

        /**
         * Checks if the given mesh is an output mesh
         *
         * @param aMeshIndex Mesh index to check
         */
        bool is_output_mesh( uint aMeshIndex ) const;

        /**
         * returns lagrange input mesh index
         */
        const Matrix< DDUMat >&
        get_lagrange_input_mesh() const
        {
            return mLagrangeInputMeshes;
        };

        /**
         * Sets the Lagrange input meshes.
         *
         * @param aLagrangeInputMeshes Lagrange input mesh indices
         */
        void
        set_lagrange_input_mesh( const Matrix< DDUMat >& aLagrangeInputMeshes )
        {
            // test if calling this function is allowed
            this->error_if_locked( "set_lagrange_input_mesh" );

            mLagrangeInputMeshes = aLagrangeInputMeshes;
        };

        /**
         * Padding size is the bigger one of mBufferSize and mMaxPolynomial.
         * In the future, filter size will be regarded here.
         *
         * @return luint number of padding elements defining aura width
         */
        auto get_padding_size() const
                -> decltype( mStaircaseBuffer );

        /**
         * returns user defined elements per direction on domain (without aura)
         *
         * @return Matrix< DDLUMat >
         */
        auto
        get_number_of_elements_per_dimension() const
                -> decltype( mNumberOfElementsPerDimension )
        {
            return mNumberOfElementsPerDimension;
        }

        /**
         * sets elements per direction on domain (without aura) according to
         * defined value
         *
         * @param[in] aNumberOfElementsPerDimension Number of elements per dimension as a vector
         */
        void set_number_of_elements_per_dimension( const Vector< uint >& aNumberOfElementsPerDimension );


        /**
         * sets elements per direction on domain (without aura) according to
         * defined value. 2D Version.
         *
         * @param[in] aElementsX Number of elements in the x direction
         * @param[in] aElementsY Number of elements in the y direction
         */
        void set_number_of_elements_per_dimension(
                luint aElementsX,
                luint aElementsY );

        /**
         * sets elements per direction on domain (without aura) according to
         * defined value. 3D Version.
         *
         * @param[in] aElementsX Number of elements in the x direction
         * @param[in] aElementsY Number of elements in the y direction
         * @param[in] aElementsZ Number of elements in the z direction
         */
        void set_number_of_elements_per_dimension(
                luint aElementsX,
                luint aElementsY,
                luint aElementsZ );

        /**
         * determines if dimension is 1D, 2D or 3D
         *
         * return luint
         */
        luint
        get_number_of_dimensions() const
        {
            return mNumberOfElementsPerDimension.size();
        }

        /**
         * returns with, height and length of specified domain
         *
         * @param[in] Matrix< DDRMat > Mat containing length in x, y and z-direction
         *
         * @return void
         */
        void set_domain_dimensions( const Vector< real >& aDomainDimensions );

        /**
         * returns with, height and length of specified domain. 2D Version
         *
         * @param[in] real dimension in X-Direction
         * @param[in] real dimension in Y-Direction
         * @return void
         */
        void set_domain_dimensions(
                real aDomainDimensionsX,
                real aDomainDimensionsY );

        /**
         * returns with, height and length of specified domain. 3D Version
         *
         * @param[in] real dimension in X-Direction
         * @param[in] real dimension in Y-Direction
         * @param[in] real dimension in Z-Direction
         * @return void
         */
        void set_domain_dimensions(
                real aDomainDimensionsX,
                real aDomainDimensionsY,
                real aDomainDimensionsZ );

        /**
         * returns with, height and length of specified domain
         *
         * @return Matrix< DDRMat >
         */
        const Vector< real >& get_domain_dimensions() const;

        /**
         * sets the coordinate of first node of calculation domain
         *
         * @param[in]  aDomainOffset   Mat containing the coordinates
         *
         * @return void
         */
        void set_domain_offset( const Matrix< DDRMat >& aDomainOffset );

        /**
         * sets the coordinate of first node of calculation domain. 2D Version.
         *
         * @param[in] aDomainOffsetX  coordinate offset in x-direction
         * @param[in] aDomainOffsetY  coordinate offset in y-direction
         *
         * @return void
         */
        void set_domain_offset(
                real aDomainOffsetX,
                real aDomainOffsetY );

        /**
         * sets the coordinate of first node of calculation domain. 3D Version.
         *
         * @param[in] aDomainOffsetX  coordinate offset in x-direction
         * @param[in] aDomainOffsetY  coordinate offset in y-direction
         * @param[in] aDomainOffsetY  coordinate offset in y-direction
         *
         * @return void
         */
        void set_domain_offset(
                real aDomainOffsetX,
                real aDomainOffsetY,
                real aDomainOffsetZ );

        /**
         * returns coordinate of first node on calculation domain
         *
         * return Matrix< DDRMat >
         */
        auto
        get_domain_offset() const
                -> decltype( mDomainOffset )
        {
            return mDomainOffset;
        }

        /**
         * Calculates which ijk range contains the calculation domain
         * of the mesh (excludes padding elements)
         *
         * return Matrix< DDLUMat >
         */
        Matrix< DDLUMat > get_domain_ijk() const;

        /**
         * Gets the current maximum polynomial order
         *
         * @return Max order
         */
        auto
        get_max_polynomial() const
                -> decltype( mMaxPolynomial )
        {
            return mMaxPolynomial;
        }

        /**
         * Sets the mesh scale factor for saving to gmsh.
         *
         * @param aScaleFactor gmsh scale factor
         */
        void
        set_gmsh_scale( real aScaleFactor )
        {
            mGmshScale = aScaleFactor;
        }

        /**
         * Gets the mesh scale factor for saving to gmsh.
         *
         * @return gmsh scale factor
         */
        auto
        get_gmsh_scale() const
                -> decltype( mGmshScale )
        {
            return mGmshScale;
        }

        /**
         * Sets the B-spline truncation flag.
         *
         * @param aTruncateBSplines Whether or not to truncate B-splines
         */
        void
        set_bspline_truncation( bool aTruncateBSplines )
        {
            mBSplineTruncationFlag = aTruncateBSplines;
        }

        /**
         * Gets if B-spline truncation is on or not.
         *
         * @return B-spline truncation flag
         */
        bool
        truncate_bsplines() const
        {
            return mBSplineTruncationFlag;
        }

        /**
         * test if input is sane
         */
        void check_sanity() const;

        /**
         * Sets the pattern used for union meshes.
         *
         * @param aUsedUnionPattern Union pattern to set
         */
        void
        set_union_pattern( uint aUsedUnionPattern )
        {
            mUnionPattern = aUsedUnionPattern;
        }

        /**
         * returns the default pattern for union meshes
         *
         * @return Union pattern
         */
        uint
        get_union_pattern() const
        {
            return mUnionPattern;
        }

        /**
         * Gets the current working pattern.
         *
         * @return Working pattern
         */
        uint
        get_working_pattern() const
        {
            return mWorkingPattern;
        }

        /**
         * Sets the current working pattern.
         *
         * @param aWorkingPattern New working pattern
         */
        void
        set_working_pattern( uint aWorkingPattern )
        {
            mWorkingPattern = aWorkingPattern;
        }

        /**
         * lock critical parameters
         */
        void lock();

        /**
         * Sets initial refinement levels
         *
         * @param aLevel Initial refinement
         */
        void
        set_initial_refinement( const moris::Matrix< DDUMat >& aLevel )
        {
            mInitialRefinementLevel = aLevel;
        }

        /**
         * Gets initial refinement levels
         *
         * @return Initial refinement
         */
        moris::Matrix< DDUMat >
        get_initial_refinement() const
        {
            return mInitialRefinementLevel;
        }

        /**
         * Sets initial refinement patterns
         *
         * @param aPatterns Initial refinement patterns
         */
        void
        set_initial_refinement_patterns( const moris::Matrix< DDUMat >& aPatterns )
        {
            mInitialRefinementPattern = aPatterns;
        }

        /**
         * Gets initial refinement patterns
         *
         * @return Initial refinement patterns
         */
        moris::Matrix< DDUMat >
        get_initial_refinement_patterns() const
        {
            return mInitialRefinementPattern;
        }

        /**
         * Gets the initial refinement for a specific activation pattern.
         *
         * @param aActivationPattern Activation pattern to target
         * @return Initial refinement level
         */
        uint
        get_initial_refinement( uint aActivationPattern ) const
        {
            sint tInitialRefinement = 0;

            for ( uint Ik = 0; Ik < mInitialRefinementPattern.numel(); Ik++ )
            {
                if ( mInitialRefinementPattern( Ik ) == aActivationPattern )
                {
                    tInitialRefinement = mInitialRefinementLevel( Ik );
                }
            }
            return tInitialRefinement;
        }

        /**
         * Sets additional Lagrange refinement level.
         *
         * @param aLevel Additional refinement level
         */
        void
        set_additional_lagrange_refinement( uint aLevel )
        {
            mAdditionalLagrangeRefinementLevel = aLevel;
        }

        /**
         * Gets additional refinement level.
         *
         * @return Additional refinement level
         */
        uint
        get_additional_lagrange_refinement() const
        {
            return mAdditionalLagrangeRefinementLevel;
        }

        /**
         * Gets the side sets to generate.
         *
         * @return Side set indices
         */
        const Matrix< DDUMat >&
        get_side_sets() const
        {
            return mSideSets;
        }

        /**
         * Sets which side sets to generate.
         *
         * @param aSideSets Side set indices
         */
        void
        set_side_sets( const Matrix< DDUMat >& aSideSets )
        {
            mSideSets = aSideSets;
        }

        /**
         * Gets if multigrid is to be used by HMR.
         *
         * @return multigrid flag
         */
        [[nodiscard]] bool
        use_multigrid() const
        {
            return mUseMultigrid;
        }

        /**
         * Gets if the aura is to be numbered by HMR.
         *
         * @return number aura flag
         */
        [[nodiscard]] bool
        use_number_aura() const
        {
            return mNumberAura;
        }

        /**
         * Gets if refinement is to be used for low level elements
         *
         * @return Low level refinement flag
         */
        bool
        use_refinement_for_low_level_elements() const
        {
            return mRefinementForLowLevelElements;
        }

        /**
         * Gets if advanced T-matrices are to be used by HMR.
         *
         * @return Advanced T-matrix flag
         */
        bool
        use_advanced_t_matrices() const
        {
            return mAdvancedTMatrices;
        }

        /**
         * Sets if HMR is to use multigrid.
         *
         * @param aUseMultigrid Multigrid flag
         */
        void
        set_multigrid( bool aUseMultigrid )
        {
            mUseMultigrid = aUseMultigrid;
        }

        /**
         * Sets if HMR is to nubmer the aura.
         *
         * @param aNumberAura Number aura flag
         */
        void
        set_number_aura( bool aNumberAura )
        {
            mNumberAura = aNumberAura;
        }

        /**
         * Sets if refinement is to be used for low level elements
         *
         * @param aRefinementForLowLevelElements Low level element refinement flag
         */
        void
        set_refinement_for_low_level_elements( const bool aRefinementForLowLevelElements )
        {
            mRefinementForLowLevelElements = aRefinementForLowLevelElements;
        }

        /**
         * Sets if HMR is to use advanced T-matrices
         *
         * @param aAdvancedTMatrices Advanced T-matrix flag
         */
        void
        set_use_advanced_t_matrices( bool aAdvancedTMatrices )
        {
            mAdvancedTMatrices = aAdvancedTMatrices;
        }

        /**
         * Sets the file name to write the background mesh to.
         *
         * @param aWriteBackgroundMesh Background mesh file name
         */
        void
        set_background_mesh_file_name( const std::string& aWriteBackgroundMesh )
        {
            mBackgroundMeshFileName = aWriteBackgroundMesh;
        }

        /**
         * Sets the file name to write the Lagrange mesh to
         *
         * @param aWriteOutputLagrangeMesh Lagrange mesh file name
         */
        void
        set_lagrange_mesh_file_name( const std::string& aWriteOutputLagrangeMesh )
        {
            mLagrangeMeshFileName = aWriteOutputLagrangeMesh;
        }

        /**
         * Sets if HMR is to write the refinement pattern to file.
         *
         * @param aWriteRefinmentPatternFile Refinement pattern output flag
         */
        void
        set_write_refinement_pattern_file_flag( bool aWriteRefinmentPatternFile )
        {
            mWriteRefinementPattern = aWriteRefinmentPatternFile;
        }

        /**
         * Sets the HDF5 file name for HMR to read in restart refinement info from.
         *
         * @param aRestartfromRefinedPattern Restart refinement pattern output file
         */
        void
        set_restart_refinement_pattern_file( const std::string& aRestartfromRefinedPattern )
        {
            mRestartFromRefinedPatternFileName = aRestartfromRefinedPattern;
        }

        /**
         * Gets the name of the file to write the background mesh to, or an empty string for no output.
         *
         * @return Background mesh file name
         */
        const std::string&
        get_background_mesh_file_name()
        {
            return mBackgroundMeshFileName;
        }

        /**
         * Gets the name of the file to write the Lagrange mesh to, or an empty string for no output.
         *
         * @return Output Lagrange mesh file name
         */
        const std::string&
        get_lagrange_mesh_file_name()
        {
            return mLagrangeMeshFileName;
        }

        /**
         * Gets if the refinement pattern is to be written to file
         *
         * @return
         */
        [[nodiscard]] bool
        write_refinement_pattern() const
        {
            return mWriteRefinementPattern;
        }

        /**
         * Gets the HDF5 file name to read refinement restart info from.
         *
         * @return Refinement restart file name
         */
        const std::string&
        get_restart_refinement_pattern_file()
        {
            return mRestartFromRefinedPatternFileName;
        }

        /**
         * Gets the maximum refinement level.
         *
         * @return Max refinement level
         */
        uint
        get_max_refinement_level() const
        {
            return mMaxRefinementLevel;
        }

        /**
         * Sets the maximum refinement level.
         *
         * @param aLevel Max refinement level
         */
        void
        set_max_refinement_level( uint aLevel )
        {
            mMaxRefinementLevel = std::min( aLevel, gMaxNumberOfLevels - 1 );
        }

        /**
         * Sets if Lagrange nodes are to be renumbered or not.
         *
         * @param aRenumberLagrangeNodes Lagrange node renumbering flag
         */
        void
        set_renumber_lagrange_nodes( bool aRenumberLagrangeNodes )
        {
            mRenumberLagrangeNodes = aRenumberLagrangeNodes;
        }

        /**
         * Gets if Lagrange nodes are to be renumbered or not.
         *
         * @return Lagrange node renumbering flag
         */
        bool
        get_renumber_lagrange_nodes() const
        {
            return mRenumberLagrangeNodes;
        }

        /**
         * Sets the VTK file name to write basis functions to.
         *
         * @param aFileName Basis function file name
         */
        void
        set_basis_fuction_vtk_file_name( const std::string& aFileName )
        {
            mBasisFunctionVtkFileName = aFileName;
        }

        /**
         * Gets the VTK file name to write basis functions to, or an empty string if they are not to be written.
         *
         * @return Basis function file name
         */
        std::string
        get_basis_fuction_vtk_file_name() const
        {
            return mBasisFunctionVtkFileName;
        }

        /**
         * Sets the refinement functions to the parameters.
         *
         * @param aRefinementFunctions Vector of refinement functions
         */
        void set_refinement_functions( const Vector< Refinement_Function >& aRefinementFunctions );

        /**
         * Get a user-defined refinement function from the parameters
         * @param aFunctionIndex
         * @return
         */
        Refinement_Function get_refinement_function( uint aFunctionIndex );

      private:
        /**
         * calls error message only if parameters are locked
         */
        void error_if_locked( const std::string& aFunctionName ) const;

        /**
         * called from set_lagrange_orders and set_bspline_orders
         */
        void update_max_polynomial_and_truncated_buffer();

        /**
         * auto setting for dimension lengths and offset
         */
        void set_default_dimensions_and_offset();
    };
}