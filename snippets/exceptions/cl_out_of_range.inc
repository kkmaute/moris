// Index 5 is out-of-range for some container.
moris::exceptions::out_of_range error(5);
std::ostringstream error_msg;

// Throw the exception.
try
{
	throw error;
}
catch(moris::exceptions::out_of_range const & e)
{
	error_msg << e.what();
}
