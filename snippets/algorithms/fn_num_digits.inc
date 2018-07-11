moris::size_t a = 1;
moris::size_t b = 10;
moris::size_t c = 100;
moris::size_t d = 1000;

moris::size_t num_digits_is_1 = moris::num_digits(a);
moris::size_t num_digits_is_2 = moris::num_digits(b);
moris::size_t num_digits_is_3 = moris::num_digits(c);
moris::size_t num_digits_is_4 = moris::num_digits(d);

// Does not work with floating point numbers.
moris::real e = 1;
moris::real f = 10;
moris::real g = 100;
moris::real h = 1000;

moris::size_t num_digits_is_324 = moris::num_digits(e);
moris::size_t num_digits_is_325 = moris::num_digits(f);
moris::size_t num_digits_is_326 = moris::num_digits(g);
moris::size_t num_digits_is_327 = moris::num_digits(h);

// Negative sign counts as a digit.
moris::sint i = -10;

moris::size_t neg_num_digits_is_3 = moris::num_digits(i);
