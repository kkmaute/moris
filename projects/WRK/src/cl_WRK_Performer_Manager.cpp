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

#include "cl_Library_IO.hpp"

namespace moris
{
    namespace wrk
    {
        //------------------------------------------------------------------------------

        // Parameter function
        typedef void ( *Parameter_Function ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );

        //------------------------------------------------------------------------------

        Performer_Manager::Performer_Manager( std::shared_ptr< Library_IO > aLibrary )
        : mLibrary(aLibrary)
        {
        }

        //------------------------------------------------------------------------------

        Performer_Manager::~Performer_Manager()
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
            Parameter_Function tHMRParameterListFunc = mLibrary->load_function<Parameter_Function>( tHMRString );
            moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
            tHMRParameterListFunc( tHMRParameterList );

            std::string tGENString = "GENParameterList";
            Parameter_Function tGENParameterListFunc = mLibrary->load_function<Parameter_Function>( tGENString );
            moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
            tGENParameterListFunc( tGENParameterList );

            // create HMR performer
            mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList( 0 )( 0 ), mLibrary );

            // create MTK performer - will be used for HMR mesh
            mMTKPerformer( 0 ) =std::make_shared< mtk::Mesh_Manager >();

            // Create GE performer
            mGENPerformer( 0 ) = std::make_shared< ge::Geometry_Engine >( tGENParameterList, mLibrary );

            // create MTK performer - will be used for XTK mesh
            mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

            // create MDL performer
            mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mLibrary, 0 );
        }

        //------------------------------------------------------------------------------

        void Performer_Manager::set_performer_cooperations()
        {
            // Set performer to HMR
            mHMRPerformer( 0 )->set_performer( mMTKPerformer( 0 ) );

            // Set performer to MDL
            mMDLPerformer( 0 )->set_performer( mMTKPerformer( 1 ) );
        }

        //------------------------------------------------------------------------------

        void Performer_Manager::create_xtk()
        {
            // Read parameter list from shared object
            Parameter_Function tXTKParameterListFunc = mLibrary->load_function<Parameter_Function>( "XTKParameterList" );
            moris::Cell< moris::Cell< ParameterList > > tXTKParameterList;
            tXTKParameterListFunc( tXTKParameterList );

            // Create XTK
            mXTKPerformer( 0 ) = std::make_shared< xtk::Model >( tXTKParameterList( 0 )( 0 ) );

            // Reset output MTK performer
            mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();
            mMDLPerformer( 0 )->set_performer( mMTKPerformer( 1 ) );

            // Set performers
            mXTKPerformer( 0 )->set_geometry_engine( mGENPerformer( 0 ).get() );
            mXTKPerformer( 0 )->set_input_performer( mMTKPerformer( 0 ) );
            mXTKPerformer( 0 )->set_output_performer( mMTKPerformer( 1 ) );
        }

        //------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
