moris::ios::clog_Sink clog_sink;

// Writes 24 characters: "A trace severity message"
moris::size_t clog_msg1 = clog_sink.write("A trace severity message", 24);

// Writes 7 characters: "A trace"
moris::size_t clog_msg2 = clog_sink.write("A trace severity message", 7);