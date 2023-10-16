// Sleep for tTestingInterval in ms.
std::this_thread::sleep_for( moris::chronos::milliseconds( tTestingInterval ) );

// Stop timer.
t.stop();

// Sleep for tTestingInterval in ms.
std::this_thread::sleep_for( moris::chronos::milliseconds( tTestingInterval ) );

// Get wall time. Time > tTestingInterval but < 2 * tTestingInterval because timer was stopped.
moris::size_t wall_time_hours        = t.toc< moris::chronos::hours >().wall;
moris::size_t wall_time_minutes      = t.toc< moris::chronos::minutes >().wall;
moris::size_t wall_time_seconds      = t.toc< moris::chronos::seconds >().wall;
moris::size_t wall_time_milliseconds = t.toc< moris::chronos::milliseconds >().wall;
moris::size_t wall_time_microseconds = t.toc< moris::chronos::microseconds >().wall;
moris::size_t wall_time_nanoseconds  = t.toc< moris::chronos::nanoseconds >().wall;
