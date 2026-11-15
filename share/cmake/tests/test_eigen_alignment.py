#!/usr/bin/env python3
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
"""Exercise Eigen's public ABI with SSE and AVX objects in the same executable.

Run: python3 share/cmake/tests/test_eigen_alignment.py /path/to/include/eigen3
"""
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap


def main():
    if sys.platform != "linux" or "avx2" not in Path("/proc/cpuinfo").read_text():
        print("Mixed SSE/AVX alignment check requires an AVX2 Linux host; skipped.")
        return
    eigen_include = Path(sys.argv[1]).resolve()
    if not (eigen_include / "Eigen/Dense").is_file():
        raise SystemExit("Provide the directory containing Eigen/Dense.")
    cmake_dir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moris-eigen-alignment-") as temporary:
        root = Path(temporary)
        (root / "CMakeLists.txt").write_text(textwrap.dedent('''
            cmake_minimum_required(VERSION 3.17)
            project(EigenAlignment LANGUAGES CXX)
            include("${MORIS_CMAKE_DIR}/moris_options.cmake")
            set(MORIS_USE_ARMA OFF)
            set(MORIS_USE_EIGEN ON)
            set(MORIS_USE_EIGEN_BLAS OFF)
            set(MORIS_USE_OPENBLAS ON)
            include("${MORIS_CMAKE_DIR}/moris_features.cmake")
            list(TRANSFORM MORIS_DEFINITIONS REPLACE "^-D" "")
            add_compile_definitions(${MORIS_DEFINITIONS})
            include_directories("${EIGEN_INCLUDE}")
            set(CMAKE_CXX_STANDARD 17)
            add_library(provider OBJECT provider.cpp)
            target_compile_options(provider PRIVATE -O0 -msse2 -mno-avx)
            add_library(driver OBJECT driver.cpp)
            target_compile_options(driver PRIVATE -O0 -mavx2)
            # Place the baseline allocator first, reproducing ELF template-symbol
            # interposition between differently compiled modules and plugins.
            add_executable(check $<TARGET_OBJECTS:provider> $<TARGET_OBJECTS:driver>)
        '''))
        (root / "provider.cpp").write_text(textwrap.dedent('''
            #include <Eigen/Core>
            #if EIGEN_MAX_ALIGN_BYTES != 32 || EIGEN_MAX_STATIC_ALIGN_BYTES != 32
            #error Missing public Eigen ABI contract
            #endif
            extern "C" void* allocate_values(std::size_t n)
            { return Eigen::internal::aligned_malloc(n); }
            extern "C" void free_values(void* p)
            { Eigen::internal::aligned_free(p); }
        '''))
        (root / "driver.cpp").write_text(textwrap.dedent('''
            #include <Eigen/Dense>
            #if EIGEN_MAX_ALIGN_BYTES != 32 || EIGEN_MAX_STATIC_ALIGN_BYTES != 32
            #error Missing public Eigen ABI contract
            #endif
            extern "C" void* allocate_values(std::size_t);
            extern "C" void free_values(void*);
            int main()
            {
                double* buffers[64];
                for (auto& p : buffers)
                {
                    p = static_cast<double*>(allocate_values(4 * sizeof(double)));
                }
                for (auto p : buffers)
                {
                    Eigen::Map<Eigen::MatrixXd> a(p, 2, 2);
                    a << -1, 1, 1, 1;
                    Eigen::MatrixXd b = a.eval();
                    Eigen::MatrixXd inverse = a.inverse();
                    if ((inverse * b - Eigen::MatrixXd::Identity(2, 2)).norm() > 1e-12)
                    {
                        return 1;
                    }
                }
                for (auto p : buffers)
                {
                    free_values(p);
                }
                return 0;
            }
        '''))
        build = root / "build"
        subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                        "-DMORIS_CMAKE_DIR=" + str(cmake_dir),
                        "-DEIGEN_INCLUDE=" + str(eigen_include)], check=True)
        subprocess.run(["cmake", "--build", str(build), "--parallel", "2"], check=True)
        subprocess.run([str(build / "check")], check=True)
    print("Mixed SSE/AVX allocation, copy, inverse, and deallocation passed.")


if __name__ == "__main__":
    main()
