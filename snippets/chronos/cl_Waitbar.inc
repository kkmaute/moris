// Initialize timer.
moris::tic t;

// Lengthy computation on container.
std::vector<moris::real> big_container(10000,0.0);

// Display indication of progress for a lengthy computation.
moris::chronos::Waitbar show_progress(big_container.size());

for (moris::uint it=0; it<big_container.size(); ++it)
{
    // Do the computation...
    std::this_thread::sleep_for(moris::chronos::microseconds(1));

    // Display progress.
    ++show_progress;
}

// Measure elapsed time.
moris::size_t wall_time_microseconds = t.toc<moris::chronos::microseconds>().wall;
