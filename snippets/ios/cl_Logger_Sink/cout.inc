moris::ios::cout_Sink cout_sink;

// Writes 33 characters: "An informational severity message"
moris::size_t cout_msg1 = cout_sink.write("An informational severity message", 33);

// Writes 16 characters: "An informational"
moris::size_t cout_msg2 = cout_sink.write("An informational severity message", 16);