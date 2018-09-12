/*
 * cl_HMR_Parameters.hpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_PARAMETERS_HPP_
#define SRC_HMR_CL_HMR_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"

#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_XML_Parser.hpp"       //CON/src
#include "cl_Param_List.hpp"       //CON/src

namespace moris
{
    namespace hmr
    {


       typedef Param_List< boost::variant< bool, sint, real, std::string > > ParameterList;

// -----------------------------------------------------------------------------

        // creates a parameter list with default inputs
        //ParameterList
        //create_parameter_list();

// -----------------------------------------------------------------------------

        // creates a parameter list with default inputs
        ParameterList
        load_parameter_list_from_xml( const std::string & aFilePath );

//--------------------------------------------------------------------------------

        /**
         * \brief This struct contains user defined settings
         */
        class Parameters
        {
           //! number of elements per direction in overall mesh, without aura
           //! 2D or 3D is determined by length of this vector
           Mat< luint > mNumberOfElementsPerDimension ;

           //! width, height and depth of domain (without aura)
           Mat< real >  mDomainDimensions;

           //! coordinate of first visible node
           Mat< real >  mDomainOffset;

           // --- Begin changable parameters.
           //     Make sure to add them to copy_selected_parameters()

           //! size of staircase buffer
           luint        mBufferSize              = 1;

           //! max polynomial to be supported
           luint        mMaxPolynomial           = 2;

           //! tells if debug flags are to be printed
           bool         mVerbose                 = true ;

           //! max surface level for refinement
           uint         mMaxSurfaceLevel = 3;

           //! max level for refinement
           uint         mMaxVolumeLevel = 2;

           //! for demo mode
           //real         mDemoKnotParameter = 1;

           //! scale factor for gmsh output
           real         mGmshScale = 1;

           //! flag telling if truncation is used
           bool         mBSplineTruncationFlag = true;

           // --- End changable parameters.

           //! tells if critical features of the settings object are locked
           bool         mParametersAreLocked = false;

           //! mesh orders, by default, a linear mesh is generated
           Mat< uint >  mLagrangeOrders = { { 1 } };

           //! mesh orders, by default, a linear mesh is generated
           Mat< uint >  mBSplineOrders = { { 1 } };

           //! defines which Lagrange mesh is associated with which refinement pattern
           Mat< uint > mLagrangePatterns = { { 0 } };

           //! defines which B-Spline mesh is associated with which refinement pattern
           Mat< uint > mBSplinePatterns = { { 0 } };

           //! Links the Lagrange mesh to a B-Spline Mesh
           Mat< uint > mLagrangeToBSpline = { { 0 } };

           //! default input pattern
           const      uint mInputPattern = 0;

           //! default output pattern
           const      uint mOutputPattern = 1;

           //! default union pattern
           const      uint mUnionPattern = 2;

           //! default pattern for output refinement
           const      uint mRefinedOutputPattern = 3;

           //! Lagrange Meshes that are used for the unity meshes
           Mat< uint >     mUnionMeshes;

           //! Lagrange Meshes that are used for the output meshes
           Mat< uint >     mOutputMeshes;

           //! Lagrange Mesh that is used for the refined output
           uint             mRefinedOutputMesh = 3;

//--------------------------------------------------------------------------------
        public:
//--------------------------------------------------------------------------------

          /*
           * trivial constructor
           */
          Parameters(){};

//--------------------------------------------------------------------------------

          /*
           * parameter list constructor
           */
          Parameters( ParameterList & aParameterList );

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
           void
           print() const;

//--------------------------------------------------------------------------------

           /**
            * returns verbosity switch
            *
            * @return bool
            */
           auto
           is_verbose() const
               -> decltype ( mVerbose )
           {
               return mVerbose;
           }

//--------------------------------------------------------------------------------

           /**
            * sets verbosity switch
            *
            * @param[in] aSwitch    true or false
            * @return void
            */
           void
           set_verbose( const bool aSwitch )
           {
               mVerbose = aSwitch;
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
           set_buffer_size( const luint & aBufferSize )
           {
               mBufferSize = aBufferSize;
           }

//--------------------------------------------------------------------------------

           /**
            * returns the buffer size
            *
            * @return luint
            */
           auto
           get_buffer_size()
               const -> decltype( mBufferSize )
           {
               return mBufferSize;
           }

//--------------------------------------------------------------------------------

           /**
            * a function which sets orders for lagrange and B-Spline meshes
            * in the most simple way
            */
           void
           set_mesh_orders_simple( const uint & aMaxOrder );

//--------------------------------------------------------------------------------
           /**
            * sets the mesh orders according to given matrix
            */
           void
           set_lagrange_orders( const Mat< uint > & aMeshOrders );

//--------------------------------------------------------------------------------


           /**
            * sets the mesh orders according to given matrix
            */
           void
           set_bspline_orders( const Mat< uint > & aMeshOrders );

//--------------------------------------------------------------------------------

           /**
            * returns a matrix with mesh orders
            */
           auto
           get_lagrange_orders() const -> decltype( mLagrangeOrders )
           {
               return mLagrangeOrders;
           }

//--------------------------------------------------------------------------------

           /**
            * returns a matrix with mesh orders
            */
           auto
           get_bspline_orders() const -> decltype( mBSplineOrders )
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
           void
           set_lagrange_patterns( const Mat< uint > & aPatterns );

//-------------------------------------------------------------------------------

           /**
            * returns a moris::Mat containing the patterns the meshes are linked to
            */
           auto
           get_lagrange_patterns() const -> decltype ( mLagrangePatterns )
           {
               return mLagrangePatterns;
           }

//-------------------------------------------------------------------------------

           /**
            * returns an entry of mLagrangeOrders
            */
           auto
           get_lagrange_order( const uint & aIndex ) const
               -> decltype( mLagrangeOrders( aIndex ) )
           {
               return mLagrangeOrders( aIndex );
           }

//-------------------------------------------------------------------------------

           /**
            * returns an entry of mLagrangePatterns
            */
           auto
           get_lagrange_pattern( const uint & aIndex ) const
               -> decltype( mLagrangePatterns( aIndex ) )
           {
               return mLagrangePatterns( aIndex );
           }


//--------------------------------------------------------------------------------

           /**
            * define which Lagrange mesh is linked to which B-Spline mesh
            */
           void
           set_lagrange_to_bspline(  const Mat< uint > & aBSplineMeshIndices );

//-------------------------------------------------------------------------------

           /**
            * returns the matrix telling which Lagrange mesh is linked with which
            * B-Spline mesh
            */
           auto
           get_lagrange_to_bspline() const -> decltype( mLagrangeToBSpline )
           {
               return mLagrangeToBSpline;
           }

//-------------------------------------------------------------------------------

           /**
            * returns an individual entry of Lagrange to B-Spline
            */
           auto
           get_lagrange_to_bspline( const uint & aIndex ) const
               -> decltype( mLagrangeToBSpline ( aIndex ) )
           {
               return mLagrangeToBSpline ( aIndex );
           }

//-------------------------------------------------------------------------------

           /**
            * sets the patterns for the B-Spline Meshes
            *
           * @param[ in ] aPattern patterns set by set_mesh_orders the B-Spline meshes refer to.
            */
           void
           set_bspline_patterns( const Mat< uint > & aPatterns );

//-------------------------------------------------------------------------------

           /**
            * returns a moris::Mat containing the patterns the meshes are linked to
            */
           auto
           get_bspline_patterns() const -> decltype ( mBSplinePatterns )
           {
               return mBSplinePatterns;
           }

//--------------------------------------------------------------------------------

           /**
            * returns an entry of mBSplinePatterns
            */
           auto
           get_bspline_pattern( const uint & aIndex ) const
               -> decltype( mBSplinePatterns( aIndex ) )
           {
               return mBSplinePatterns( aIndex );
           }

//--------------------------------------------------------------------------------

           /**
            * returns an entry of mBSplineOrders
            */
           auto
           get_bspline_order( const uint & aIndex ) const
               -> decltype( mBSplineOrders( aIndex ) )
           {
               return mBSplineOrders( aIndex );
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
           get_union_mesh( const uint & aOrder ) const
           {
               return mUnionMeshes( aOrder-1 );
           }

//--------------------------------------------------------------------------------

           /**
            * returns the index of the defined Lagrange output mesh for a specified order
            */
           uint
           get_output_mesh( const uint & aOrder ) const
           {
               return mOutputMeshes( aOrder-1 );
           }

//--------------------------------------------------------------------------------

           /**
            * returns the mesh for refined output
            */
           auto
           get_refined_output_mesh() const -> decltype( mRefinedOutputMesh )
           {
               return mRefinedOutputMesh;
           }

//--------------------------------------------------------------------------------


           /**
            * sets the maximum polynomial degree to given value
            *
            * @param[in] aMaxPolynomial
            *
            * @return void
            */
           //void
           //set_max_polynomial( const luint & aMaxPolynomial ) ;

//--------------------------------------------------------------------------------

           /**
            * Padding size is the bigger one of mBufferSize and mMaxPolynomial.
            * In the future, filter size will be regarded here.
            *
            * @return luint number of padding elements defining aura width
            */
           auto
           get_padding_size() const
               -> decltype ( mBufferSize ) ;

//--------------------------------------------------------------------------------

           /**
            * returns user defined elements per direction on domain (without aura)
            *
            * @return Mat<luint>
            */
           auto
           get_number_of_elements_per_dimension() const
               -> decltype ( mNumberOfElementsPerDimension )
           {
               return mNumberOfElementsPerDimension ;
           }

//--------------------------------------------------------------------------------
           /**
            * sets elements per direction on domain (without aura) according to
            * defined value
            *
            * @param[in] aNumberOfElementsPerDimension Mat<luint>
            *
            * @return void
            */
           void
           set_number_of_elements_per_dimension(
                   const Mat<luint> & aNumberOfElementsPerDimension );

//--------------------------------------------------------------------------------

           /**
            * sets elements per direction on domain (without aura) according to
            * defined value. 2D Version.
            *
            * @param[in] aElementsX
            * @param[in] aElementsY
            * @return void
            */
           void
           set_number_of_elements_per_dimension(
                           const luint & aElementsX,
                           const luint & aElementsY );

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
           void
           set_number_of_elements_per_dimension(
                   const luint & aElementsX,
                   const luint & aElementsY,
                   const luint & aElementsZ );

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
            * @param[in] Mat<real> Mat containing length in x, y and z-direction
            *
            * @return void
            */
           void
           set_domain_dimensions( const Mat<real> & aDomainDimensions );

//-------------------------------------------------------------------------------

           /**
            * returns with, height and length of specified domain. 2D Version
            *
            * @param[in] real dimension in X-Direction
            * @param[in] real dimension in Y-Direction
            * @return void
            */
           void
           set_domain_dimensions(
                   const real & aDomainDimensionsX,
                   const real & aDomainDimensionsY );

//-------------------------------------------------------------------------------

           /**
            * returns with, height and length of specified domain. 3D Version
            *
            * @param[in] real dimension in X-Direction
            * @param[in] real dimension in Y-Direction
            * @param[in] real dimension in Z-Direction
            * @return void
            */
           void
           set_domain_dimensions(
                   const real & aDomainDimensionsX,
                   const real & aDomainDimensionsY,
                   const real & aDomainDimensionsZ);

//-------------------------------------------------------------------------------

           /**
            * returns with, height and length of specified domain
            *
            * @return Mat<real>
            */
           Mat< real >
           get_domain_dimensions() const ;

//-------------------------------------------------------------------------------

           /**
            * sets the coordinate of first node of calculation domain
            *
            * @param[in]  aDomainOffset   Mat containing the coordinates
            *
            * @return void
            */
           void
           set_domain_offset( const Mat<real> & aDomainOffset );

//-------------------------------------------------------------------------------

           /**
            * sets the coordinate of first node of calculation domain. 2D Version.
            *
            * @param[in] aDomainOffsetX  coordinate offset in x-direction
            * @param[in] aDomainOffsetY  coordinate offset in y-direction
            *
            * @return void
            */
           void
           set_domain_offset(
                   const real & aDomainOffsetX,
                   const real & aDomainOffsetY );

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
           void
           set_domain_offset(
                   const real & aDomainOffsetX,
                   const real & aDomainOffsetY,
                   const real & aDomainOffsetZ );

//-------------------------------------------------------------------------------
           /**
            * returns coordinate of first node on calculation domain
            *
            * return Mat<real>
            */
           auto
           get_domain_offset() const -> decltype( mDomainOffset )
           {
               return mDomainOffset;
           }

//-------------------------------------------------------------------------------

           /**
            * Calculates which ijk range contains the calculation domain
            * of the mesh (excludes padding elements)
            *
            * return Mat<luint>
            */
           Mat<luint>
           get_domain_ijk() const;

//-------------------------------------------------------------------------------

           void
           set_max_surface_level( const uint & aLevel )
           {
               mMaxSurfaceLevel = aLevel;
           }

//-------------------------------------------------------------------------------

           auto
           get_max_surface_level() const -> decltype ( mMaxSurfaceLevel )
           {
               return mMaxSurfaceLevel;
           }

//-------------------------------------------------------------------------------

           void
           set_max_volume_level( const uint & aLevel )
           {
               mMaxVolumeLevel = aLevel;
           }

//-------------------------------------------------------------------------------

           auto
           get_max_volume_level() const -> decltype ( mMaxVolumeLevel )
           {
               return mMaxVolumeLevel;
           }

//-------------------------------------------------------------------------------

           auto
           get_max_polynomial() const -> decltype ( mMaxPolynomial )
           {
               return mMaxPolynomial;
           }

//-------------------------------------------------------------------------------

           //void
           //set_demo_knot_parameter( const real & aParam )
           //{
           //    mDemoKnotParameter = aParam;
           //}

//-------------------------------------------------------------------------------

           //auto
           //get_demo_knot_parameter() const -> decltype ( mDemoKnotParameter )
           //{
           //    return mDemoKnotParameter;
           //}

//-------------------------------------------------------------------------------

           void
           set_gmsh_scale( const real& aScaleFactor )
           {
               mGmshScale = aScaleFactor;
           }

//-------------------------------------------------------------------------------
           auto
           get_gmsh_scale() const -> decltype( mGmshScale )
           {
               return  mGmshScale;
           }

//-------------------------------------------------------------------------------

           void
           set_bspline_truncation( const bool aSwitch )
           {
               mBSplineTruncationFlag = aSwitch;

               if ( aSwitch )
               {
                   mBufferSize = mMaxPolynomial;
               }
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
            * sets the values for  mLagrangePatterns and mBSplinePatterns
            * to default values
            */
           void
           set_mesh_order( const uint & aInterpolationOrder );

//-------------------------------------------------------------------------------

           /**
            * test if input is sane
            */
           void
           check_sanity() const;

//-------------------------------------------------------------------------------

           /**
            * returns the default pattern for input meshes
            */
           auto
           get_input_pattern() const -> decltype( mInputPattern )
           {
               return mInputPattern;
           }

//-------------------------------------------------------------------------------

           /**
            * returns the default pattern for output meshes
            */
           auto
           get_output_pattern() const -> decltype( mOutputPattern )
           {
               return mOutputPattern;
           }

//-------------------------------------------------------------------------------

           /**
            * returns the default pattern for union meshes
            */
           auto
           get_union_pattern() const -> decltype( mUnionPattern )
           {
               return mUnionPattern;
           }

//-------------------------------------------------------------------------------

           /**
            * returns the default pattern for union meshes
            */
           auto
           get_refined_output_pattern() const -> decltype( mRefinedOutputPattern )
           {
               return mRefinedOutputPattern;
           }

//-------------------------------------------------------------------------------

           /**
            * Copy selected parameters from other parameter object
            * Note that not all parameters can be copied
            */
           void
           copy_selected_parameters( const Parameters & aParameters );

//-------------------------------------------------------------------------------

           /**
            * Copy selected parameters from other parameter list
            * Note that not all parameters can be copied
            */
           void
           copy_selected_parameters( ParameterList & aParameterList );

//-------------------------------------------------------------------------------

           /**
            * lock critical parameters
            */
           void
           lock();

//-------------------------------------------------------------------------------
        private:
//-------------------------------------------------------------------------------

           /**
            * returns an error message for an invalid parameter
            */
           void
           error( const std::string & aMessage ) const;

//-------------------------------------------------------------------------------

           /**
            * calls error message only if parameters are locked
            */
           void
           error_if_locked( const std::string & aFunctionName  ) const;

//-------------------------------------------------------------------------------

           /**
            * called from set_lagrange_orders and set_bspline_orders
            */
           void
           update_max_polynomial_and_truncated_buffer();

//-------------------------------------------------------------------------------

           /**
            * auto setting for dimension lengths and offset
            */
           void
           set_default_dimensions_and_offset();

//-------------------------------------------------------------------------------

           /**
            * converts a string to a real matrix
            */
           void
           string_to_mat( const std::string & aString, Mat< real > & aMat ) const;

//-------------------------------------------------------------------------------

           /**
            * converts a string to an luint matrix
            */
           void
           string_to_mat( const std::string & aString, Mat< luint > & aMat ) const;

//-------------------------------------------------------------------------------

        }; /* Parameters */
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_PARAMETERS_HPP_ */
