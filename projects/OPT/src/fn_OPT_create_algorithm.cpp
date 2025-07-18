/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_create_algorithm.cpp
 *
 */

#include "fn_OPT_create_algorithm.hpp"
#include "cl_OPT_Algorithm_GCMMA.hpp"
#include "cl_OPT_Algorithm_SQP.hpp"
#include "cl_OPT_Algorithm_LBFGS.hpp"
#include "cl_OPT_Algorithm_Sweep.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::opt
{
    std::shared_ptr< Algorithm >
    create_algorithm( const Parameter_List& aAlgorithmParameterList )
    {
        std::string tAlgorithmName = aAlgorithmParameterList.get< std::string >( "algorithm" );
        if ( tAlgorithmName == "gcmma" )
        {
            return std::make_shared< OptAlgGCMMA >( aAlgorithmParameterList );
        }
        else if ( tAlgorithmName == "sqp" )
        {
            return std::make_shared< Algorithm_SQP >( aAlgorithmParameterList );
        }
        else if ( tAlgorithmName == "lbfgs" )
        {
            return std::make_shared< Algorithm_LBFGS >( aAlgorithmParameterList );
        }
        else if ( tAlgorithmName == "sweep" )
        {
            return std::make_shared< Algorithm_Sweep >( aAlgorithmParameterList );
        }
        else
        {
            MORIS_ERROR( false, "%s is not recognized as a valid Algorithm type in fn_OPT_create_algorithm.", tAlgorithmName.c_str() );
            return nullptr;
        }
    }
}    // namespace moris::opt
