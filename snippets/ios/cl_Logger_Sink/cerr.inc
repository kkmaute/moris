moris::ios::cerr_Sink cerr_sink;

// Writes 25 characters: "An error severity message"
moris::size_t cerr_msg1 = cerr_sink.write("An error severity message", 25);

// Writes 8 characters: "An error"
moris::size_t cerr_msg2 = cerr_sink.write("An error severity message", 8);