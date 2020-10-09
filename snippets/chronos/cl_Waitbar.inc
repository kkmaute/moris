// Initialize timer.
moris::tic t;

// Lengthy computation on container.
std::vector<moris::real> big_container(10000,0.0);

// Display indication of progress for a lengthy computation.
moris::chronos::Waitbar show_progress(big_container.size());

for (auto it: big_container)
{
    // Do the computation...
    std::this_thread::sleep_for(moris::chronos::microseconds(1));

    // Do some dummy operation on variable "it" to avoid it being flagged as unused
    it += 1.0;

    // Display progress.
    ++show_progress;
}

// Measure elapsed time.
moris::size_t wall_time_microseconds = t.toc<moris::chronos::microseconds>().wall;
