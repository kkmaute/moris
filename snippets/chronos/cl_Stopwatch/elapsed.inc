// Sleep for 1ms.
std::this_thread::sleep_for(moris::chronos::milliseconds(1));

// Get wall time in different formats.
// 1. 0.001447640s wall, 0.010000000s user + 0.000000000s system = 0.010000000s CPU (690.8%)
// 2. 0.001526597s wall
// 3. 0.010000000s user
// 4. 0.000000000s system
// 5. 0.010000000s CPU
// 6. (651.4%)
std::string wall_time_1 = t.elapsed();
std::string wall_time_2 = t.elapsed("%ws wall");
std::string wall_time_3 = t.elapsed("%us user");
std::string wall_time_4 = t.elapsed("%ss system");
std::string wall_time_5 = t.elapsed("%ts CPU");
std::string wall_time_6 = t.elapsed("(%p%)");
