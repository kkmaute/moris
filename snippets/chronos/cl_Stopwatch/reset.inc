// Reset timer.
t.reset();

// Sleep for 1 second anew.
std::this_thread::sleep_for(moris::chronos::seconds(1));

// Get wall time. Time > 1s but < 2s.
moris::real wall_time_hours = t.toc<moris::chronos::hours>().wall;
moris::real wall_time_minutes = t.toc<moris::chronos::minutes>().wall;
moris::real wall_time_seconds = t.toc<moris::chronos::seconds>().wall;
moris::real wall_time_milliseconds = t.toc<moris::chronos::milliseconds>().wall;
moris::real wall_time_microseconds = t.toc<moris::chronos::microseconds>().wall;
moris::real wall_time_nanoseconds = t.toc<moris::chronos::nanoseconds>().wall;
