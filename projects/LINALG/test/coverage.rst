LINALG unit-test coverage
=========================

``coverage.json`` maps every public header in ``../src/`` to registered unit-test
sources. Backend implementation headers are exercised through the public MORIS
interface. The tutorials and performance timings do not count as assertion
coverage. Run the inventory check after adding or removing an API or test::

    python3 share/cmake/tests/test_linalg_coverage.py

This is an API/test inventory, not a claim of 100% instrumented line or branch
coverage. It does not enumerate every possible scalar type, shape, backend error,
or third-party implementation branch. Successful compilation alone does not
validate a numerical result; the mapped Catch cases contain result assertions.

Coverage added by the audit
--------------------------

* ``UT_LINALG_Matrix_API.cpp`` checks iterators, column-major pointer access,
  capacity, vector construction, resizing, member reshape with padding and
  truncation, complex member transpose, const spans/views, output row/column
  overloads, deep copies, scalar assignment, indexed real/complex extrema,
  const matrix comparisons, stream output, and the explicit move helper.
* ``UT_LINALG_Function_Edges.cpp`` checks reindexing with repeated indices,
  singular/rank-deficient and rectangular matrices, complex/rectangular trace,
  empty reductions, vector reduction orientation, concatenation with aliasing
  and empty operands, sort-index permutations, exact nonzeros and duplicate
  indices, intersection index correspondence, unsorted interpolation samples
  and aliased queries, rectangular identities, linspace endpoints, negative
  diagonal offsets, rectangular solves, and wide QR/LU factorizations.
* ``UT_LINALG_Print.cpp`` asserts the printed content of matrix, row-vector,
  indexed, initializer-list, real-matrix, matrix-container, and pointer-container
  printers, including the convenience macros. Stream/file redirection is
  restored before assertions and on exceptions.
* ``UT_LINALG_Complex_Functions.cpp`` checks upper Cholesky, inverse, determinant,
  absolute values, nonfinite values, rectangular SVD with real singular values,
  and a rectangular Sylvester equation against an independently known solution.
* Binary I/O tests cover empty and complex matrices, inaccessible paths, exact
  4-MiB blocks, and partial tails. New filenames include MPI size and rank to
  avoid collisions between process-count variants.

Existing assertion coverage
---------------------------

The per-function ``fn_*.cpp`` and per-operator ``op_*.cpp`` files cover the basic
arithmetic, comparisons, predicates, decompositions, solve, interpolation-space
utilities, and I/O. ``fn_diag.cpp`` covers ``diag_vec`` as well as ``diag_mat``.
``UT_LINALG_Backend_Parity.cpp`` and ``UT_LINALG_Complex_Parity.cpp`` provide
cross-backend numerical references and view/storage checks.
``UT_LINALG_Expressions.cpp`` checks native expression types, deferred evaluation,
compound operations, destination storage reuse, product/transpose/reshape aliases,
complex nonconjugating free transpose, and native-operand operator lookup.
``op_move.cpp`` checks moves, copies, vector relocation, and borrowed-view moves.

Configuration coverage
----------------------

Run the Catch suite with both ``MORIS_USE_ARMA`` and ``MORIS_USE_EIGEN``, enabling
only one backend per build. Check both runtime-check configurations:
``MORIS_HAVE_DEBUG=ON`` and ``OFF``. Bounds-error expectations must be guarded by
``MORIS_HAVE_DEBUG``, not by an Armadillo-only preprocessor definition.
Numerical tests remain active in both configurations. Symbolic information and
runtime checks are separate CMake options.

For a configured build::

    cmake --build <build> --target LINALG-test --parallel 3
    ctest --test-dir <build> -R '^LINALG-test-(1|2|4)-procs$' --output-on-failure

The inventory check does not replace these backend/configuration runs.

Audit validation
----------------

The initial audit covered 113 Catch cases per backend (28 added by that
audit). GCC 13.3.1 builds with Armadillo 14.4.1 and Eigen 3.4.0 passed:

==================  =======================  =============
Backend/build       MORIS_HAVE_DEBUG         MPI ranks
==================  =======================  =============
Armadillo Release   OFF                      1, 2, 4
Armadillo Debug     ON                       1
Eigen Release       OFF                      1, 2, 4
Eigen Debug         ON                       1, 2, 4
==================  =======================  =============

Subsequent application-integration checks expanded the suite to 118 cases.
All 118 passed on one MPI rank with Armadillo Release (runtime checks OFF)
and Eigen Debug (runtime checks ON). The additional regressions cover native
row/column expressions and overlap, fixed/dynamic value conversions, scalar
arithmetic, vector extrema, and transpose/vectorization composition.

The complete enabled Eigen Release build, including unit-test executables
and examples, also passes with runtime checks OFF. Eleven diffusion and
max-damage FEM cases pass with that library; the ten diffusion cases also
pass with Armadillo Debug. Additional checks cover third-order geometry
interpolation, SDF transformations and 2D/3D vertex rotations, and QUAD4
subdivision. Mixed fixed/dynamic addition and subtraction retain native
expression evaluation in the LINALG regressions.
Subsequent complete enabled Armadillo and Eigen Release builds each passed all
147 registered CTest entries with runtime checks OFF. These runs include LINALG,
MTK, INT, PETSc/BDDC, and application examples on their configured MPI counts.
Math libraries used one thread per rank, inherited PETSc options were cleared,
and MPI CPU binding was disabled for concurrent test launches. The bedding and
adaptive-boxbeam variants use isolated working directories and callback-library
paths. The application fixes cover default-empty horizontal concatenation,
stable cached operands during lazy FEM expressions, explicit MORIS scalar
reductions of native subviews, and PETSc 3.23.2 retained Dirichlet coarse-basis
initialization. Exhaustive instrumented branch coverage remains separate
validation work.

The guarded Eigen LU path adds square-system regressions for extreme uniform
scaling, independently scaled RHS columns, native expressions, input/output
aliases, complex coefficients, singular inconsistent systems, and ill-conditioned
systems. An Eigen-specific Wilkinson pivot-growth case checks recovery when an
acceptable condition estimate accompanies an inaccurate LU solution; another
case preserves COD's rank truncation near machine precision. All 122 Eigen and
121 Armadillo LINALG cases passed with runtime checks OFF, and the complete Eigen
build and all 147 CTest entries passed after introducing solver selection.
Single-rank, single-thread measurements reduced the repeated-solve phase from
about 18 seconds to about 10 seconds, including conditioning and per-RHS
backward-error checks. Timings remain diagnostics rather than test requirements.

Reduction ownership regressions expand the LINALG suite to 124 Eigen cases and
122 Armadillo cases. Evaluated product operands are shared immutably when Eigen
copies reduction expressions; coefficient-wise operands retain deferred reads.
The tests cover copied/moved expression lifetime, composition, and destination
aliases. All 124 Eigen cases also pass with runtime debug checks ON.

The corresponding FEM checks use conservative Frobenius-error/spectral-lower-bound
ratios instead of SVD-based assertion norms. Their tolerance is unchanged, and
normalization keeps the bounds valid at scales of 1e-200 and 1e200. L2 residual
and Jacobian lumping are compared against independently assembled dense full
space-time references for zero, mixed, and full lumping with zero/nonunit
weights. The formulation's coefficient order and consistent time coupling are
retained. Both complete builds and 147-entry CTest suites pass; the final INT
suites contain 274 passing cases per backend. Controlled single-rank runs reduce
Eigen INT time from about 84.3 to 67.1 seconds, with Armadillo taking about
65.5 seconds. The five L2 cases fall from 4.44 to 0.70 seconds with Eigen.

The added tests exposed and now protect against default Armadillo sort-index
conversion failures, deprecated indexed extrema calls, integer vector-size
constructor selection, complex member-transpose conjugation differences, and
row-vector intersection shape differences. Release Eigen compilation also
checks that address-only alias detection accepts uninitialized destination
storage without GCC's ``-Wmaybe-uninitialized`` warning. Existing bounds tests
now use the backend-independent MORIS runtime-check option.
