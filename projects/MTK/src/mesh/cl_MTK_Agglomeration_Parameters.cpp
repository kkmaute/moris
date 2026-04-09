/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.cpp
 *
 */

#include "cl_MTK_Agglomeration_Parameters.hpp"

namespace moris::mtk
{

    Agglomeration_Parameters::Agglomeration_Parameters(
            real aAgglomerationExponent,
            real aAgglomerationShift )
            : mExp( aAgglomerationExponent )
            , mShift( aAgglomerationShift )
    {
        MORIS_ASSERT( (uint)mExp % 2 == 0, "Agglomeration_Parameters - Exponent must be even to ensure violation value is positive" );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Agglomeration_Parameters::get_exp() const
    {
        return mExp;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Agglomeration_Parameters::get_shift() const
    {
        return mShift;
    }

    //--------------------------------------------------------------------------------------------------------------

    Agglomeration_Parameters_Constant::Agglomeration_Parameters_Constant(
            real aAgglomerationExponent,
            real aAgglomerationShift,
            real aAgglomerationReference )
            : Agglomeration_Parameters( aAgglomerationExponent, aAgglomerationShift )
            , mRef( aAgglomerationReference ) {};

    //--------------------------------------------------------------------------------------------------------------

    Agglomeration_Parameters_Variable::Agglomeration_Parameters_Variable(
            real                             aAgglomerationExponent,
            real                             aAgglomerationShift,
            Agglomeration_Reference_Function aAgglomerationReference )
            : Agglomeration_Parameters( aAgglomerationExponent, aAgglomerationShift )
            , mRefFunc( aAgglomerationReference ) {};

    //--------------------------------------------------------------------------------------------------------------

    real Agglomeration_Parameters_Constant::get_ref( const Matrix< DDRMat >& aCoords, uint aIndex ) const
    {
        return mRef;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Agglomeration_Parameters_Variable::get_ref( const Matrix< DDRMat >& aCoords, uint aIndex ) const
    {
        return mRefFunc( aCoords, aIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::unique_ptr< Agglomeration_Parameters > create_agglomeration_parameters( const Parameter_List& aParameterList, const std::shared_ptr< Library_IO > aLibrary )
    {
        std::string tReferenceFunctionName = aParameterList.get< std::string >( "agglomeration_reference_function_name" );

        if ( tReferenceFunctionName == "" )
        {
            // No reference function specified, use constant reference value
            return std::make_unique< Agglomeration_Parameters_Constant >(
                    aParameterList.get< real >( "agglomeration_exponent" ),
                    aParameterList.get< real >( "agglomeration_shift" ),
                    aParameterList.get< real >( "agglomeration_reference" ) );
        }
        else
        {
            // Reference function specified, look it up and use it
            Agglomeration_Reference_Function tRefFunc = aLibrary->load_function< Agglomeration_Reference_Function >( tReferenceFunctionName );

            return std::make_unique< Agglomeration_Parameters_Variable >(
                    aParameterList.get< real >( "agglomeration_exponent" ),
                    aParameterList.get< real >( "agglomeration_shift" ),
                    tRefFunc );
        }
    }
}    // namespace moris::mtk