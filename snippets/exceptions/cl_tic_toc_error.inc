moris::exceptions::tic_toc_error error;
std::ostringstream error_msg;

// Throw the exception.
try
{
	throw error;
}
catch(moris::exceptions::tic_toc_error & e)
{
	error_msg << e.what();
}
