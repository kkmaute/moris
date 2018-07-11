// 5 == 4.
// Size 5 is invalid size, should be 4.
moris::exceptions::invalid_size error(5);
std::ostringstream error_msg;

// Throw the exception.
try
{
	throw error;
}
catch(moris::exceptions::invalid_size & e)
{
	error_msg << e.what();
}
