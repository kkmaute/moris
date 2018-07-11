// Read mode.
std::string read_mode = moris::ios::to_stdio(std::ios::in);

// Write mode.
std::string write_mode = moris::ios::to_stdio(std::ios::out | std::ios::trunc);

// Append mode.
std::string append_mode = moris::ios::to_stdio(std::ios::app | std::ios::out);

// Read/write mode.
std::string readwrite_mode = moris::ios::to_stdio(std::ios::in | std::ios::out);

// Read/truncate mode.
std::string trunc_mode = moris::ios::to_stdio(std::ios::in | std::ios::out | std::ios::trunc);

// Read/update mode.
std::string update_mode = moris::ios::to_stdio(std::ios::app | std::ios::in | std::ios::out);