/*
 * cl_HMR_Parameters.hpp
 *
 *  Created on: May 5, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_PARAMETERS_HPP_
#define SRC_HMR_CL_HMR_PARAMETERS_HPP_

#include <cstdio>

#include "assert.hpp"
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------

        /**
         * \brief This struct contains user defined settings
         */
        class Parameters
        {

           //! number of elements per direction in overall mesh, without aura
           //! 2D or 3D is determined by length of this vector
           Mat< luint > mNumberOfElementsPerDimension ;

           //! size of staircase buffer
           luint        mBufferSize              = 1;

           //! max polynomial to be supported
           luint        mMaxPolynomial           = 2;

           //! tells if debug flags are to be printed
           bool         mVerbose                 = true ;

           //! width, height and depth of domain (without aura)
           Mat< real >  mDomainDimensions;

           //! coordinate of first visible node
           Mat< real >  mDomainOffset            = { { 0 }, { 0 }, { 0 } };

           //! max surface level for refinement
           uint         mMaxSurfaceLevel = 3;

           //! max level for refinement
           uint         mMaxVolumeLevel = 2;

           //! for demo mode
           real         mDemoKnotParameter = 1;

           //! scale factor for gmsh output
           real         mGmshScale = 1;

           //! flag telling if truncation is used
           bool         mBSplineTruncationFlag = true;

           //! tells if critical features of the settings object are locked
           bool         mParametersAreLocked = false;
//--------------------------------------------------------------------------------
        public:
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
            * sets the maximum polynomial degree to given value
            *
            * @param[in] aMaxPolynomial
            *
            * @return void
            */
           void
           set_max_polynomial( const luint & aMaxPolynomial ) ;

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
           set_domain_dimensions(
                   const Mat<real> & aDomainDimensions )
           {
                   mDomainDimensions = aDomainDimensions;
           }

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
           set_domain_offset( const Mat<real> & aDomainOffset )
           {
               mDomainOffset = aDomainOffset;
           }

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

           void
           set_demo_knot_parameter( const real & aParam )
           {
               mDemoKnotParameter = aParam;
           }

//-------------------------------------------------------------------------------

           auto
           get_demo_knot_parameter() const -> decltype ( mDemoKnotParameter )
           {
               return mDemoKnotParameter;
           }

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
           }


//-------------------------------------------------------------------------------

           bool
           truncate_bsplines() const
           {
               return mBSplineTruncationFlag;
           }

//-------------------------------------------------------------------------------
        }; /* Parameters */
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_PARAMETERS_HPP_ */
