/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Parameters.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_PARAMETERS_HPP_
#define SRC_HMR_CL_HMR_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"
#include "assert.hpp"

#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_XML_Parser.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace mtk
    {
        class Cell;
        class Field;
    }    // namespace mtk
    namespace hmr
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

        // -----------------------------------------------------------------------------

        // fixme: to be deleted soon
        // creates a parameter list with default inputs
        void load_hmr_parameter_list_from_xml( const std::string& aFilePath,
                ParameterList&                                    aParameterList );

        //--------------------------------------------------------------------------------

        /**
         * \brief This struct contains user defined settings
         */
        class Parameters
        {
            //! Processor decomposition method.  1=Original MPI decomp (min processor interface). 2=Min mesh interface. 0=Manually Defined
            uint mProcDecompMethod = 1;

            //! Processor layout if mProcDecompMethod is 0 (user defined here). Product MUST = # of processors used. Can be 1, 2, or 3 dimensions.
            Matrix< DDUMat > mProcessorDimensions = { { 2 }, { 2 } };

            //! number of elements per direction in overall mesh, without aura
            //! 2D or 3D is determined by length of this vector
            Matrix< DDLUMat > mNumberOfElementsPerDimension = { { 2 }, { 2 } };

            //! width, height and depth of domain (without aura)
            Matrix< DDRMat > mDomainDimensions;

            //! coordinate of first visible node
            Matrix< DDRMat > mDomainOffset;

            // --- Begin changable parameters.
            //     Make sure to add them to copy_selected_parameters()

            //! size of staircase buffer
            luint mRefinementBuffer = 1;
            luint mStaircaseBuffer  = 1;

            //! max polynomial to be supported
            luint mMaxPolynomial = 1;

            //! scale factor for gmsh output
            real mGmshScale = 1;

            //! flag telling if truncation is used
            bool mBSplineTruncationFlag = true;

            // --- End changable parameters.

            //! tells if critical features of the settings object are locked
            bool mParametersAreLocked = false;

            //! mesh orders, by default, a linear mesh is generated
            Matrix< DDUMat > mLagrangeOrders = { { 1 } };

            //! mesh orders, by default, a linear mesh is generated
            Matrix< DDUMat > mBSplineOrders = { { 1 } };

            //! defines which Lagrange mesh is associated with which refinement pattern
            Matrix< DDUMat > mLagrangePatterns = { { 0 } };

            //! defines which B-Spline mesh is associated with which refinement pattern
            Matrix< DDUMat > mBSplinePatterns = { { 0 } };

            //! defines which B-Spline mesh is associated with which lagrange mesh
            Cell< Matrix< DDSMat > > mLagrangeToBSplineMesh;

            //! maps input orders with B-Splines
            //           Matrix< DDUMat> mBSplineInputMap;
            //           Matrix< DDUMat> mBSplineOutputMap;

            //! default input pattern     //FIXME delete these
            const uint mBSplineInputPattern  = 0;
            const uint mLagrangeInputPattern = 1;

            //! default output pattern
            const uint mBSplineOutputPattern  = 2;
            const uint mLagrangeOutputPattern = 3;

            //! default union pattern
            uint mUnionPattern = gNumberOfPatterns - 1;

            //! default pattern for iterative refinement
            uint mWorkingPattern = gNumberOfPatterns - 2;

            //! default pattern for output refinement
            //           const      uint mRefinedOutputPattern = 5;

            //! Map Lagrange Meshes that are used for the unity meshes
            //! position 0: first order,
            //! position 1: second order,
            //! position 2: third order
            Matrix< DDUMat > mUnionMeshes;

            //! Lagrange Meshes that are used for the output meshes
            Cell< Matrix< DDUMat > > mOutputMeshes     = { { { 0 } } };
            Cell< std::string >      mOutputMesheNames = { { { "" } } };

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

            std::string mWriteBackgroundMesh             = "";
            std::string mWriteOutputLagrangeMesh         = "";
            std::string mWriteOutputLagrangeMeshToExodus = "";

            std::string mBasisFunctionVtkFileName = "";

            bool mWriteRefinementPatternFileFlag = false;

            std::string mRestartfromRefinedPattern = "";

            //! maximum level for refinement. Default value is specified
            //! by global constant
            uint mMaxRefinementLevel = gMaxNumberOfLevels - 1;

            //! Generate information about mesh refinement level interrelation
            bool mRefinementInterrelation = false;

            //! Renumber Lagrange Nodes
            bool mRenumberLagrangeNodes = false;

            // HMR user defined refinement function
            Cell< Refinement_Function > mRefinementFunctions;

            //--------------------------------------------------------------------------------

          public:
            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------

            /**
             * sets processor decomposition method
             *
             * @param[in] aProcDecompMethod uint
             *
             * @return void
             */
            void set_processor_decomp_method( const uint& aProcDecompMethod );

            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------
            /**
             * sets processor dimensions for manually defined processor decomposition
             * method, (mProcDecompMethod==3). Matters for this case only.
             *
             * @param[in] aProcessorDimensions Matrix< DDUMat >
             *
             * @return void
             */
            void set_processor_dimensions( const Matrix< DDUMat >& aProcessorDimensions );

            //--------------------------------------------------------------------------------

            /*
             * trivial constructor
             */
            Parameters(){};

            //--------------------------------------------------------------------------------

            /*
             * parameter list constructor
             */
            Parameters( ParameterList&                   aParameterList,
                    std::shared_ptr< moris::Library_IO > aLibrary );

            //--------------------------------------------------------------------------------

            /*
             * trivial destructor
             */
            ~Parameters(){};
            //--------------------------------------------------------------------------------
            /**
             * prints user settings passed to HMR
             *
             * @return void
             */
            void print() const;

            //--------------------------------------------------------------------------------

            /**
             * sets verbosity switch
             *
             * @param[in] aSwitch    true or false
             * @return void
             */
            void
            set_severity_level( const sint aSwitch )
            {
                gLogger.set_severity_level( aSwitch );
            }

            //--------------------------------------------------------------------------------

            sint
            get_severity_level() const
            {
                return gLogger.get_severity_level();
            };

            //--------------------------------------------------------------------------------

            /**
             * sets the buffer size to given value
             *
             * @param[in] aBufferSize
             *
             * @return void
             */
            void
            set_refinement_buffer( const luint& aBufferSize )
            {
                mRefinementBuffer = aBufferSize;
            }

            //--------------------------------------------------------------------------------

            /**
             * sets the buffer size to given value
             *
             * @param[in] aBufferSize
             *
             * @return void
             */
            void
            set_staircase_buffer( const luint& aBufferSize )
            {
                mStaircaseBuffer = aBufferSize;
            }

            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------
            /**
             * sets the mesh orders according to given matrix
             */
            void set_lagrange_orders( const Matrix< DDUMat >& aMeshOrders );

            //--------------------------------------------------------------------------------

            /**
             * sets the mesh orders according to given matrix
             */
            void set_bspline_orders( const Matrix< DDUMat >& aMeshOrders );

            //--------------------------------------------------------------------------------

            /**
             * returns a matrix with mesh orders
             */
            auto
            get_lagrange_orders() const
                    -> decltype( mLagrangeOrders )
            {
                return mLagrangeOrders;
            }

            //--------------------------------------------------------------------------------

            /**
             * returns a matrix with mesh orders
             */
            auto
            get_bspline_orders() const
                    -> decltype( mBSplineOrders )
            {
                return mBSplineOrders;
            }

            //-------------------------------------------------------------------------------

            /**
             * sets the patterns for the Lagrange Meshes
             *
             * @param[ in ] aPattern patterns set by set_mesh_orders the Lagrange meshes refer to.
             *
             */
            void set_lagrange_patterns( const Matrix< DDUMat >& aPatterns );

            //-------------------------------------------------------------------------------

            /**
             * returns a moris::Mat containing the patterns the meshes are linked to
             */
            auto
            get_lagrange_patterns() const
                    -> decltype( mLagrangePatterns )
            {
                return mLagrangePatterns;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns an entry of mLagrangeOrders
             */
            auto
            get_lagrange_order( const uint& aIndex ) const
                    -> decltype( mLagrangeOrders( aIndex ) )
            {
                return mLagrangeOrders( aIndex );
            }

            //-------------------------------------------------------------------------------

            /**
             * returns an entry of mLagrangePatterns
             */
            auto
            get_lagrange_pattern( const uint& aIndex ) const
                    -> decltype( mLagrangePatterns( aIndex ) )
            {
                return mLagrangePatterns( aIndex );
            }

            //-------------------------------------------------------------------------------

            /**
             * sets the patterns for the B-Spline Meshes
             *
             * @param[ in ] aPattern patterns set by set_mesh_orders the B-Spline meshes refer to.
             */
            void set_bspline_patterns( const Matrix< DDUMat >& aPatterns );

            //-------------------------------------------------------------------------------

            /**
             * returns a moris::Mat containing the patterns the meshes are linked to
             */
            auto
            get_bspline_patterns() const
                    -> decltype( mBSplinePatterns )
            {
                return mBSplinePatterns;
            }

            //--------------------------------------------------------------------------------

            /**
             * returns an entry of mBSplinePatterns
             */
            auto
            get_bspline_pattern( const uint& aIndex ) const
                    -> decltype( mBSplinePatterns( aIndex ) )
            {
                return mBSplinePatterns( aIndex );
            }

            //--------------------------------------------------------------------------------

            /**
             * returns an entry of mBSplineOrders
             */
            auto
            get_bspline_order( const uint& aIndex ) const
                    -> decltype( mBSplineOrders( aIndex ) )
            {
                return mBSplineOrders( aIndex );
            }

            //--------------------------------------------------------------------------------

            /**
             * returns an entry of mBSplineOrders
             */
            void
            set_lagrange_to_bspline_mesh( const Cell< Matrix< DDSMat > > aLagrangeToBSplineMesh )
            {
                mLagrangeToBSplineMesh = aLagrangeToBSplineMesh;
            }

            Matrix< DDSMat >
            get_lagrange_to_bspline_mesh( const uint& aLagrangeMeshIndex ) const
            {
                return mLagrangeToBSplineMesh( aLagrangeMeshIndex );
            }

            //--------------------------------------------------------------------------------

            /**
             * returns the number of B-Spline meshes
             */
            uint
            get_number_of_bspline_meshes() const
            {
                return mBSplineOrders.length();
            }

            //--------------------------------------------------------------------------------

            /**
             * returns the number of Lagrange meshes
             */
            uint
            get_number_of_lagrange_meshes() const
            {
                return mLagrangeOrders.length();
            }
            //--------------------------------------------------------------------------------

            /**
             * returns the index of the defined Lagrange union mesh for a specified order
             */
            uint
            get_union_mesh( const uint& aOrder ) const
            {
                return mUnionMeshes( aOrder - 1 );
            }

            //--------------------------------------------------------------------------------

            /**
             * returns the index of the defined Lagrange output mesh for a specified order
             */
            const Cell< Matrix< DDUMat > >&
            get_output_mesh() const
            {
                return mOutputMeshes;
            }

            //--------------------------------------------------------------------------------

            /**
             * set which lagrange meshes are used for an output
             */
            void
            set_output_meshes( const Cell< Matrix< DDUMat > >& aOutputMeshes )
            {
                // test if calling this function is allowed
                this->error_if_locked( "set_output_meshes" );

                mOutputMeshes = aOutputMeshes;
            };

            //--------------------------------------------------------------------------------

            /**
             * set output mesh names
             */
            void
            set_output_mesh_names( const Cell< std::string >& aOutputMesheNames )
            {
                // test if calling this function is allowed
                this->error_if_locked( "set_output_meshes_names" );

                mOutputMesheNames = aOutputMesheNames;
            }

            //--------------------------------------------------------------------------------

            /**
             * returns the index of the defined Lagrange output mesh names for a specified order
             */
            const Cell< std::string >&
            get_output_mesh_names() const
            {
                return mOutputMesheNames;
            }

            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------

            /**
             * returns true if mesh name exists
             */
            bool
            get_mesh_name_exists( const std::string& aName ) const
            {
                return mOutputNameToIndexMap.key_exists( aName );
            }

            //--------------------------------------------------------------------------------

            /**
             * checks if this mesh index belongs to an output mesh
             */
            bool is_output_mesh( const uint aMeshIndex ) const;

            //--------------------------------------------------------------------------------

            /**
             * returns lagrange input mesh index
             */
            const Matrix< DDUMat >&
            get_lagrange_input_mesh() const
            {
                return mLagrangeInputMeshes;
            };

            //--------------------------------------------------------------------------------

            /**
             * set lagrange input mesh index
             */
            void
            set_lagrange_input_mesh( const Matrix< DDUMat >& aLagrangeInputMeshes )
            {
                // test if calling this function is allowed
                this->error_if_locked( "set_output_meshes" );

                mLagrangeInputMeshes = aLagrangeInputMeshes;
            };

            //--------------------------------------------------------------------------------

            /**
             * returns lagrange input mesh index
             */
            const Matrix< DDUMat >&
            get_bspline_input_mesh() const
            {
                return mBSplineInputMeshes;
            }

            //--------------------------------------------------------------------------------

            /**
             * set lagrange input mesh index
             */
            void
            set_bspline_input_mesh( const Matrix< DDUMat >& aBSplineInputMeshes )
            {
                // test if calling this function is allowed
                this->error_if_locked( "set_output_meshes" );

                mBSplineInputMeshes = aBSplineInputMeshes;
            };

            //--------------------------------------------------------------------------------

            /**
             * sets the maximum polynomial degree to given value
             *
             * @param[in] aMaxPolynomial
             *
             * @return void
             */
            // void
            // set_max_polynomial( const luint & aMaxPolynomial ) ;

            //--------------------------------------------------------------------------------

            /**
             * Padding size is the bigger one of mBufferSize and mMaxPolynomial.
             * In the future, filter size will be regarded here.
             *
             * @return luint number of padding elements defining aura width
             */
            auto get_padding_size() const
                    -> decltype( mStaircaseBuffer );

            //--------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------
            /**
             * sets elements per direction on domain (without aura) according to
             * defined value
             *
             * @param[in] aNumberOfElementsPerDimension Matrix< DDLUMat >
             *
             * @return void
             */
            void set_number_of_elements_per_dimension( const Matrix< DDLUMat >& aNumberOfElementsPerDimension );

            //--------------------------------------------------------------------------------

            /**
             * sets elements per direction on domain (without aura) according to
             * defined value. 2D Version.
             *
             * @param[in] aElementsX
             * @param[in] aElementsY
             * @return void
             */
            void set_number_of_elements_per_dimension( const luint& aElementsX,
                    const luint&                                    aElementsY );

            //--------------------------------------------------------------------------------

            /**
             * sets elements per direction on domain (without aura) according to
             * defined value. 3D Version.
             *
             * @param[in] aElementsX
             * @param[in] aElementsY
             * @param[in] aElementsZ
             * @return void
             */
            void set_number_of_elements_per_dimension( const luint& aElementsX,
                    const luint&                                    aElementsY,
                    const luint&                                    aElementsZ );

            //-------------------------------------------------------------------------------

            /**
             * determines if dimension is 1D, 2D or 3D
             *
             * return luint
             */
            luint
            get_number_of_dimensions() const
            {
                return mNumberOfElementsPerDimension.length();
            }

            //-------------------------------------------------------------------------------

            /**
             * returns with, height and length of specified domain
             *
             * @param[in] Matrix< DDRMat > Mat containing length in x, y and z-direction
             *
             * @return void
             */
            void set_domain_dimensions( const Matrix< DDRMat >& aDomainDimensions );

            //-------------------------------------------------------------------------------

            /**
             * returns with, height and length of specified domain. 2D Version
             *
             * @param[in] real dimension in X-Direction
             * @param[in] real dimension in Y-Direction
             * @return void
             */
            void set_domain_dimensions( const real& aDomainDimensionsX,
                    const real&                     aDomainDimensionsY );

            //-------------------------------------------------------------------------------

            /**
             * returns with, height and length of specified domain. 3D Version
             *
             * @param[in] real dimension in X-Direction
             * @param[in] real dimension in Y-Direction
             * @param[in] real dimension in Z-Direction
             * @return void
             */
            void set_domain_dimensions( const real& aDomainDimensionsX,
                    const real&                     aDomainDimensionsY,
                    const real&                     aDomainDimensionsZ );

            //-------------------------------------------------------------------------------

            /**
             * returns with, height and length of specified domain
             *
             * @return Matrix< DDRMat >
             */
            Matrix< DDRMat > get_domain_dimensions() const;

            //-------------------------------------------------------------------------------

            /**
             * sets the coordinate of first node of calculation domain
             *
             * @param[in]  aDomainOffset   Mat containing the coordinates
             *
             * @return void
             */
            void set_domain_offset( const Matrix< DDRMat >& aDomainOffset );

            //-------------------------------------------------------------------------------

            /**
             * sets the coordinate of first node of calculation domain. 2D Version.
             *
             * @param[in] aDomainOffsetX  coordinate offset in x-direction
             * @param[in] aDomainOffsetY  coordinate offset in y-direction
             *
             * @return void
             */
            void set_domain_offset( const real& aDomainOffsetX,
                    const real&                 aDomainOffsetY );

            //-------------------------------------------------------------------------------

            /**
             * sets the coordinate of first node of calculation domain. 3D Version.
             *
             * @param[in] aDomainOffsetX  coordinate offset in x-direction
             * @param[in] aDomainOffsetY  coordinate offset in y-direction
             * @param[in] aDomainOffsetY  coordinate offset in y-direction
             *
             * @return void
             */
            void set_domain_offset( const real& aDomainOffsetX,
                    const real&                 aDomainOffsetY,
                    const real&                 aDomainOffsetZ );

            //-------------------------------------------------------------------------------
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

            //-------------------------------------------------------------------------------

            /**
             * Calculates which ijk range contains the calculation domain
             * of the mesh (excludes padding elements)
             *
             * return Matrix< DDLUMat >
             */
            Matrix< DDLUMat > get_domain_ijk() const;

            //-------------------------------------------------------------------------------

            auto
            get_max_polynomial() const
                    -> decltype( mMaxPolynomial )
            {
                return mMaxPolynomial;
            }

            //-------------------------------------------------------------------------------

            // void
            // set_demo_knot_parameter( const real & aParam )
            //{
            //     mDemoKnotParameter = aParam;
            // }

            //-------------------------------------------------------------------------------

            // auto
            // get_demo_knot_parameter() const -> decltype ( mDemoKnotParameter )
            //{
            //     return mDemoKnotParameter;
            // }

            //-------------------------------------------------------------------------------

            void
            set_gmsh_scale( const real& aScaleFactor )
            {
                mGmshScale = aScaleFactor;
            }

            //-------------------------------------------------------------------------------
            auto
            get_gmsh_scale() const
                    -> decltype( mGmshScale )
            {
                return mGmshScale;
            }

            //-------------------------------------------------------------------------------

            void
            set_bspline_truncation( const bool aSwitch )
            {
                mBSplineTruncationFlag = aSwitch;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the flag that tells if the truncation flag is set
             */
            bool
            truncate_bsplines() const
            {
                return mBSplineTruncationFlag;
            }

            //-------------------------------------------------------------------------------

            /**
             * test if input is sane
             */
            void check_sanity() const;

            //-------------------------------------------------------------------------------

            /**
             * returns the default pattern for input meshes
             */
            uint
            get_bspline_input_pattern() const
            {
                return mBSplineInputPattern;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the default pattern for input meshes
             */
            uint
            get_lagrange_input_pattern() const
            {
                return mLagrangeInputPattern;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the default pattern for output meshes
             */
            uint
            get_bspline_output_pattern() const
            {
                return mBSplineOutputPattern;
            }

            /**
             * returns the default pattern for output meshes
             */
            uint
            get_lagrange_output_pattern() const
            {
                return mLagrangeOutputPattern;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the default pattern for union meshes
             */
            void
            set_union_pattern( const uint aUsedUnionPattern )
            {
                mUnionPattern = aUsedUnionPattern;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the default pattern for union meshes
             */
            uint
            get_union_pattern() const
            {
                return mUnionPattern;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns the working pattern
             */
            uint
            get_working_pattern() const
            {
                return mWorkingPattern;
            }

            void
            set_working_pattern( const uint aWorkingPattern )
            {
                mWorkingPattern = aWorkingPattern;
            }

            //-------------------------------------------------------------------------------
            /**
             * Copy selected parameters from other parameter object
             * Note that not all parameters can be copied
             */
            void copy_selected_parameters( const Parameters& aParameters );

            //-------------------------------------------------------------------------------

            /**
             * Copy selected parameters from other parameter list
             * Note that not all parameters can be copied
             */
            void copy_selected_parameters( ParameterList& aParameterList );

            //-------------------------------------------------------------------------------

            /**
             * lock critical parameters
             */
            void lock();

            //-------------------------------------------------------------------------------

            void
            set_initial_refinement( const moris::Matrix< DDUMat >& aLevel )
            {
                mInitialRefinementLevel = aLevel;
            }

            //-------------------------------------------------------------------------------

            moris::Matrix< DDUMat >
            get_initial_refinement() const
            {
                return mInitialRefinementLevel;
            }

            //-------------------------------------------------------------------------------

            void
            set_initial_refinement_patterns( const moris::Matrix< DDUMat >& aPatterns )
            {
                mInitialRefinementPattern = aPatterns;
            }

            //-------------------------------------------------------------------------------

            moris::Matrix< DDUMat >
            get_initial_refinement_patterns() const
            {
                return mInitialRefinementPattern;
            }

            //-------------------------------------------------------------------------------

            uint
            get_initial_refinement( uint tActivationPattern ) const
            {
                sint tInitialRefinement = 0;

                for ( uint Ik = 0; Ik < mInitialRefinementPattern.numel(); Ik++ )
                {
                    if ( mInitialRefinementPattern( Ik ) == tActivationPattern )
                    {
                        tInitialRefinement = mInitialRefinementLevel( Ik );
                    }
                }
                return tInitialRefinement;
            }

            //-------------------------------------------------------------------------------

            void
            set_additional_lagrange_refinement( const uint& aLevel )
            {
                mAdditionalLagrangeRefinementLevel = aLevel;
            }

            //-------------------------------------------------------------------------------

            uint
            get_additional_lagrange_refinement() const
            {
                return mAdditionalLagrangeRefinementLevel;
            }

            //-------------------------------------------------------------------------------

            //           Matrix< DDUMat> get_bspline_input_map() const
            //           {
            //               return mBSplineInputMap;
            //           }

            //-------------------------------------------------------------------------------

            //           Matrix< DDUMat> get_bspline_output_map() const
            //           {
            //               return mBSplineOutputMap;
            //           }

            //-------------------------------------------------------------------------------

            //           void set_bspline_input_map( const Matrix< DDUMat> & aBSplineInputMap )
            //           {
            //               mBSplineInputMap = aBSplineInputMap;
            //           }

            //-------------------------------------------------------------------------------

            //           void set_bspline_output_map( const Matrix< DDUMat> & aBSplineOutputMap )
            //           {
            //               mBSplineOutputMap = aBSplineOutputMap;
            //           }

            //-------------------------------------------------------------------------------

            const Matrix< DDUMat >&
            get_side_sets() const
            {
                return mSideSets;
            }

            //-------------------------------------------------------------------------------

            void
            set_side_sets( const Matrix< DDUMat >& aSideSets )
            {
                mSideSets = aSideSets;
            }

            //-------------------------------------------------------------------------------

            /**
             * returns a string with the specified side set ordinals
             */
            std::string get_side_sets_as_string() const;

            //-------------------------------------------------------------------------------

            bool
            use_multigrid() const
            {
                return mUseMultigrid;
            }

            //-------------------------------------------------------------------------------

            bool
            use_number_aura() const
            {
                return mNumberAura;
            }

            //-------------------------------------------------------------------------------

            bool
            use_refinement_for_low_level_elements() const
            {
                return mRefinementForLowLevelElements;
            }

            //-------------------------------------------------------------------------------

            bool
            use_advanced_t_matrices() const
            {
                return mAdvancedTMatrices;
            }

            //-------------------------------------------------------------------------------

            void
            set_multigrid( const bool aSwitch )
            {
                mUseMultigrid = aSwitch;
            }

            //-------------------------------------------------------------------------------

            void
            set_number_aura( const bool aSwitch )
            {
                mNumberAura = aSwitch;
            }

            //-------------------------------------------------------------------------------

            void
            set_refinement_for_low_level_elements( const bool aSwitch )
            {
                mRefinementForLowLevelElements = aSwitch;
            }

            //-------------------------------------------------------------------------------

            void
            set_use_advanced_t_matrices( const bool aSwitch )
            {
                mAdvancedTMatrices = aSwitch;
            }

            //-------------------------------------------------------------------------------

            void
            set_write_background_mesh( const std::string& aWriteBackgroundMesh )
            {
                mWriteBackgroundMesh = aWriteBackgroundMesh;
            }

            //-------------------------------------------------------------------------------

            void
            set_write_output_lagrange_mesh( const std::string& aWriteOutputLagrangeMesh )
            {
                mWriteOutputLagrangeMesh = aWriteOutputLagrangeMesh;
            }

            //-------------------------------------------------------------------------------

            void
            set_write_output_lagrange_mesh_to_exodus( const std::string& aWriteOutputLagrangeMesh )
            {
                mWriteOutputLagrangeMeshToExodus = aWriteOutputLagrangeMesh;
            }

            //-------------------------------------------------------------------------------

            void
            set_write_refinement_pattern_file_flag( bool aWriteRefinmentPatternFile )
            {
                mWriteRefinementPatternFileFlag = aWriteRefinmentPatternFile;
            }

            //-------------------------------------------------------------------------------

            void
            set_restart_refinement_pattern_file( const std::string& aRestartfromRefinedPattern )
            {
                mRestartfromRefinedPattern = aRestartfromRefinedPattern;
            }

            //-------------------------------------------------------------------------------

            const std::string&
            get_write_background_mesh()
            {
                return mWriteBackgroundMesh;
            }

            //-------------------------------------------------------------------------------

            const std::string&
            get_write_output_lagrange_mesh()
            {
                return mWriteOutputLagrangeMesh;
            }

            //-------------------------------------------------------------------------------

            const std::string&
            get_write_output_lagrange_mesh_to_exodus()
            {
                return mWriteOutputLagrangeMeshToExodus;
            }

            //-------------------------------------------------------------------------------

            bool
            get_write_refinement_pattern_file_flag()
            {
                return mWriteRefinementPatternFileFlag;
            }

            //-------------------------------------------------------------------------------

            const std::string&
            get_restart_refinement_pattern_file()
            {
                return mRestartfromRefinedPattern;
            }

            //-------------------------------------------------------------------------------

            uint
            get_max_refinement_level() const
            {
                return mMaxRefinementLevel;
            }

            //-------------------------------------------------------------------------------

            void
            set_max_refinement_level( const uint aLevel )
            {
                mMaxRefinementLevel = std::min( aLevel, gMaxNumberOfLevels - 1 );
            }
            //-------------------------------------------------------------------------------

            bool
            get_refinement_interrelation() const
            {
                return mRefinementInterrelation;
            }

            //-------------------------------------------------------------------------------

            void
            set_renumber_lagrange_nodes( const bool aSwitch )
            {
                mRenumberLagrangeNodes = aSwitch;
            }

            //-------------------------------------------------------------------------------

            bool
            get_renumber_lagrange_nodes() const
            {
                return mRenumberLagrangeNodes;
            }

            //-------------------------------------------------------------------------------

            void
            set_basis_fuction_vtk_file_name( const std::string aFileName )
            {
                mBasisFunctionVtkFileName = aFileName;
            }

            //-------------------------------------------------------------------------------

            std::string
            get_basis_fuction_vtk_file_name() const
            {
                return mBasisFunctionVtkFileName;
            }

            //-------------------------------------------------------------------------------

            void
            set_refinement_interrelation( const bool aSwitch )
            {
                mRefinementInterrelation = aSwitch;
            }

            /**
             *
             * @param aRefinementFunctions
             */
            void set_refinement_functions( Cell< Refinement_Function > aRefinementFunctions );

            /**
             * Get a user-defined refinement function from the parameters
             * @param aFunctionIndex
             * @return
             */
            Refinement_Function get_refinement_function( uint aFunctionIndex );

            //-------------------------------------------------------------------------------

          private:
            //-------------------------------------------------------------------------------

            /**
             * returns an error message for an invalid parameter
             */
            void error( const std::string& aMessage ) const;

            //-------------------------------------------------------------------------------

            /**
             * calls error message only if parameters are locked
             */
            void error_if_locked( const std::string& aFunctionName ) const;

            //-------------------------------------------------------------------------------

            /**
             * called from set_lagrange_orders and set_bspline_orders
             */
            void update_max_polynomial_and_truncated_buffer();

            //-------------------------------------------------------------------------------

            /**
             * auto setting for dimension lengths and offset
             */
            void set_default_dimensions_and_offset();

            //-------------------------------------------------------------------------------

            //           void set_mesh_orders( const Matrix< DDUMat > & aBSplineOrders,
            //                                 const Matrix< DDUMat > & aLagrangeOrders );

            //-------------------------------------------------------------------------------
        }; /* Parameters */

        // -----------------------------------------------------------------------------

        // creates a parameter list with default options
        //        ParameterList create_hmr_parameter_list();

        // -----------------------------------------------------------------------------

        /**
         * creates a parameter list from a parameter object
         */
        ParameterList create_hmr_parameter_list( const Parameters* aParameters );

        // -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_PARAMETERS_HPP_ */

