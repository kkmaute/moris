

#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_HMR.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_WRK_Performer_Manager.hpp"


#include "fn_Exec_load_user_library.hpp"

namespace moris
{
namespace wrk
{
//------------------------------------------------------------------------------
Performer_Manager::Performer_Manager( std::shared_ptr< Library_IO > aLibrary ) :mLibrary(aLibrary)
{

}

//------------------------------------------------------------------------------

void Performer_Manager::initialize_performers()
{
    mHMRPerformer.resize( 1 );
    mGENPerformer.resize( 1 );
    mXTKPerformer.resize( 1 );
    mMTKPerformer.resize( 2 );
    mMDLPerformer.resize( 1 );

    // load the HMR parameter list
    std::string tHMRString = "HMRParameterList";
    MORIS_PARAMETER_FUNCTION tHMRParameterListFunc = mLibrary->load_parameter_file( tHMRString );
    moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
    tHMRParameterListFunc( tHMRParameterList );

    std::string tGENString = "GENParameterList";
    MORIS_PARAMETER_FUNCTION tGENParameterListFunc = mLibrary->load_parameter_file( tGENString );
    moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
    tGENParameterListFunc( tGENParameterList );

    std::string tXTKString = "XTKParameterList";
    MORIS_PARAMETER_FUNCTION tXTKParameterListFunc = mLibrary->load_parameter_file( tXTKString );
    moris::Cell< moris::Cell< ParameterList > > tXTKParameterList;
    tXTKParameterListFunc( tXTKParameterList );

    mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList( 0 )( 0 ) );

    mMTKPerformer( 0 ) =std::make_shared< mtk::Mesh_Manager >();

    // Create GE with parameter list
    mGENPerformer( 0 ) = std::make_shared< ge::GEN_Geometry_Engine >( tGENParameterList(0)(0) );

    mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

    mXTKPerformer( 0 ) = std::make_shared< xtk::Model >( tXTKParameterList( 0 )( 0 ) );

}

void Performer_Manager::set_performer_cooperations()
{
	mHMRPerformer( 0 )->set_performer( mMTKPerformer( 0 ) );

    mGENPerformer( 0 )->set_performer( mHMRPerformer( 0 ) );
    mGENPerformer( 0 )->set_library( mLibrary );

    mXTKPerformer( 0 )->set_geometry_engine( mGENPerformer( 0 ).get() );
    mXTKPerformer( 0 )->set_performer( mMTKPerformer( 0 ) );

}

//------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
