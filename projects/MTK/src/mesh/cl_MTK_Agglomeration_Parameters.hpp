/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.hpp
 *
 */
#pragma once

#include "cl_MTK_Enums.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_Library_IO.hpp"

namespace moris::mtk
{
    using Agglomeration_Reference_Function = real ( * )( const Matrix< DDRMat >&, uint );    // Pointer to shape diameter reference function that takes facet coords and index and returns the reference value

    //--------------------------------------------------------------------------------------------------------------

    class Agglomeration_Parameters
    {
        real mExp;      // Exponent controls sharpness of agglomeration
        real mShift;    // Shifts the value by this much

      public:
        Agglomeration_Parameters(
                real aAgglomerationExponent = 2.0,
                real aAgglomerationShift    = 0.0 );

        virtual real get_ref( const Matrix< DDRMat >& aCoords, uint aIndex = 0 ) const = 0;

        /**
         * Gets the exponent for the agglomeration function
         */
        real get_exp() const;

        /**
         * Gets the shift for the agglomeration function
         */
        real get_shift() const;
    };

    //--------------------------------------------------------------------------------------------------------------

    class Agglomeration_Parameters_Constant : public Agglomeration_Parameters
    {
      private:
        real mRef;

      public:
        Agglomeration_Parameters_Constant(
                real aAgglomerationExponent  = 2.0,
                real aAgglomerationShift     = 0.0,
                real aAgglomerationReference = 1.0 );

        virtual real get_ref( const Matrix< DDRMat >& aCoords, uint aIndex = 0 ) const override;
    };

    //--------------------------------------------------------------------------------------------------------------

    class Agglomeration_Parameters_Variable : public Agglomeration_Parameters
    {
      private:
        Agglomeration_Reference_Function mRefFunc;

      public:
        Agglomeration_Parameters_Variable(
                real                             aAgglomerationExponent  = 2.0,
                real                             aAgglomerationShift     = 0.0,
                Agglomeration_Reference_Function aAgglomerationReference = nullptr );
        /**
         * Evaluates the reference function for the target for agglomeration at the given coordinates and index
         */
        virtual real get_ref( const Matrix< DDRMat >& aCoords, uint aIndex = 0 ) const override;
    };

    /**
     * Helper function to create agglomeration parameters struct from parameter list.
     * Checks the parameter list for the type of agglomeration parameters and creates the appropriate struct.
     *
     */
    std::unique_ptr< Agglomeration_Parameters > create_agglomeration_parameters( const Parameter_List& aParameterList, const std::shared_ptr< Library_IO > aLibrary );

}    // namespace moris::mtk
