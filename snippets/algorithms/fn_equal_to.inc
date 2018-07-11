moris::real pi = 3.14159265358979323846264338327950;

moris::bool_t pi_equal_16digits = moris::equal_to(pi, 3.1415926535897932);
moris::bool_t pi_equal_12digits = moris::equal_to(pi, 3.141592653589);
moris::bool_t pi_equal_8digits  = moris::equal_to(pi, 3.14159265, 1.0e+08);
moris::bool_t pi_equal_4digits  = moris::equal_to(pi, 3.1415, 1.0e+12);
moris::bool_t pi_equal_2digits  = moris::equal_to(pi, 3.14, 1.0e+16);

moris::bool_t pi_not_equal_8digits  = moris::equal_to(pi, 3.14159265);
moris::bool_t pi_not_equal_4digits  = moris::equal_to(pi, 3.1415);
moris::bool_t pi_not_equal_2digits  = moris::equal_to(pi, 3.14);
