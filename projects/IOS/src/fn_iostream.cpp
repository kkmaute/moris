#include "fn_iostream.hpp" // IOS/src

// C++ header files.
#include <iostream>

// Third-party header files.
#include <boost/iostreams/stream.hpp>

// Namespaces.
namespace io = boost::iostreams;

namespace moris
{
    // Stream buffer sinks.
    io::stream_buffer<moris::ios::clog_Sink> clog_buf((moris::ios::clog_Sink()));
    io::stream_buffer<moris::ios::cout_Sink> cout_buf((moris::ios::cout_Sink()));
    io::stream_buffer<moris::ios::cerr_Sink> cerr_buf((moris::ios::cerr_Sink()));

    // Output streams.
    std::ostream clog(&clog_buf);
    std::ostream cout(&cout_buf);
    std::ostream cerr(&cerr_buf);
}
