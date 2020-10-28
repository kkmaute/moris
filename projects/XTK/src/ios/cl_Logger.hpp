/*
 * cl_Logger.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: ktdoble
 */

#ifndef SRC_IOS_CL_LOGGER_HPP_
#define SRC_IOS_CL_LOGGER_HPP_

#include <mpi.h>

//#include <boost/log/expressions.hpp>
//#include <boost/log/sources/global_logger_storage.hpp>
//#include <boost/log/support/date_time.hpp>
//#include <boost/log/trivial.hpp>
//#include <boost/log/utility/setup.hpp>

#include "fn_num_digits.hpp"

#include <iostream>
#define XTK_LOG_FUNCTION() std::cout
#define XTK_INFO  std::cout
#define XTK_WARN  std::cerr
#define XTK_ERROR std::cerr

//#define XTK_LOG_FUNCTION() BOOST_LOG_FUNCTION()
//#define XTK_INFO  BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::info)
//#define XTK_WARN  BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::warning)
//#define std::cout BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::error)

////Narrow-char thread-safe logger.
//typedef boost::log::sources::severity_logger_mt<boost::log::trivial::severity_level> logger_t;
//
////declares a global logger with a custom initialization
//BOOST_LOG_GLOBAL_LOGGER(my_logger, logger_t)
//
//namespace attrs   = boost::log::attributes;
//namespace expr    = boost::log::expressions;
//namespace logging = boost::log;
//
////Defines a global logger initialization routine
//BOOST_LOG_GLOBAL_LOGGER_INIT(my_logger, logger_t)
//{
//    logger_t lg;
//
//    int tProcRank = 0;
//    int tProcSize = 0;
//    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
//    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);
//
//    size_t tNumDigits = xtk::num_digits(tProcSize-1);
//
//    logging::add_common_attributes();
//
//    std::string time_stamp_suffix = ".%Y-%m-%d_%H:%M:%S";
//
//    std::string tFileName = "log/XTK_LOG"+time_stamp_suffix+"."+std::to_string(tProcRank)+".log";
//
//    logging::add_file_log(
//            boost::log::keywords::file_name = tFileName,
//            boost::log::keywords::format = (
//                    expr::stream << expr::format_date_time<     boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
//                    << " [" << expr::attr<     boost::log::trivial::severity_level >("Severity") << "]: "
//                    << "[Rank: " << std::dec << std::setw(tNumDigits) << std::setfill('0') << tProcRank << "]"
//                    << expr::smessage
//            )
//    );
//
//    logging::add_console_log(
//            std::cout,
//            boost::log::keywords::format = (
//                    expr::stream << expr::format_date_time<     boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
//                    << " [" << expr::attr<     boost::log::trivial::severity_level >("Severity") << "]: "
//                    << "[Rank: " << std::dec << std::setw(tNumDigits) << std::setfill('0') << tProcRank << "] "
//                    << expr::smessage
//            )
//    );
//
//    logging::core::get()->set_filter
//    (
//        logging::trivial::severity >= logging::trivial::info
//    );
//
//    return lg;
//}



#endif /* SRC_IOS_CL_LOGGER_HPP_ */
