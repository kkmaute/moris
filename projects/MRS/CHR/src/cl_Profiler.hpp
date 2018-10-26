/*
 * cl_Profiler.hpp
 *
 *  Created on: Oct 21, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_
#define PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_

#include <string>

#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

namespace moris
{
    class Profiler
    {
        const std::string mLogfile;
//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------

        Profiler( const std::string & aLogfile ) : mLogfile( aLogfile )
        {
#ifdef WITHGPERFTOOLS
            ProfilerStart( aLogfile.c_str() );
#endif
        }

//------------------------------------------------------------------------------

        Profiler() : Profiler( "/tmp/gprofmoris.log" )
        {

        }

//------------------------------------------------------------------------------

        stop()
        {
#ifdef WITHGPERFTOOLS
    ProfilerStop();
#endif
        }

//------------------------------------------------------------------------------
    };

}



#endif /* PROJECTS_MRS_CHR_SRC_CL_PROFILER_HPP_ */
