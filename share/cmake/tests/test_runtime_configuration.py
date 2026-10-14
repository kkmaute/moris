#!/usr/bin/env python3
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
"""Check runtime prefixes, cache updates, Intel parity, and target isolation."""
import os
from pathlib import Path
import subprocess
import tempfile

cmake_dir = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="moris-runtime-") as temporary:
    root = Path(temporary)
    for prefix, directory in (("first", "lib"), ("second", "lib64")):
        folder = root / prefix / directory
        folder.mkdir(parents=True)
        for library in ("z", "ssl", "crypto", "gfortran", "dl", "ifcore", "ifcoremt", "rt"):
            (folder / ("lib" + library + ".so")).touch()
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(Runtime LANGUAGES CXX)
set(CONFIG_INSTALL_DIR share)
set(MORIS moris)
include("${MORIS_CMAKE_DIR}/moris_runtime.cmake")
get_directory_property(global_links LINK_LIBRARIES)
if(global_links)
    message(FATAL_ERROR "Runtime libraries leaked into directory requirements")
endif()
get_target_property(runtime_links moris_runtime INTERFACE_LINK_LIBRARIES)
foreach(library IN ITEMS MORIS_LIB_Z MORIS_LIB_Crypto MORIS_LIB_Ssl MORIS_LIB_GFortran)
    if(NOT "${${library}}" IN_LIST runtime_links)
        message(FATAL_ERROR "Missing runtime requirement: ${library}")
    endif()
endforeach()
''')
    for intel in ("OFF", "ON"):
        build = root / ("build-" + intel)
        for prefix, directory in (("first", "lib"), ("second", "lib64"), (None, "lib64"), ("missing", "lib")):
            env = dict(os.environ)
            for var in ("ZLIB_DIR", "OPENSSL_DIR", "DL_DIR", "RT_DIR", "INTEL_RUNTIME_DIR"):
                env.pop(var, None)
                if prefix is not None:
                    env[var] = str(root / prefix)
            options = ["-DMORIS_LIB_GFortran=" + str(root / "first/lib/libgfortran.so")]
            if prefix == "first":
                # Migration must preserve old overrides when no prefix overrides them.
                options += ["-DMORIS_Lib_Crypto=" + str(root / "first/lib/libcrypto.so")]
            result = subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                                     "-DMORIS_CMAKE_DIR=" + str(cmake_dir),
                                     "-DMORIS_USE_INTEL=" + intel, *options],
                                    env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    universal_newlines=True)
            if prefix == "missing":
                assert result.returncode and "Cannot find" in result.stdout, result.stdout
            else:
                assert result.returncode == 0, result.stdout
                cache = (build / "CMakeCache.txt").read_text()
                expected = root / (prefix or "second") / directory / "libz.so"
                assert "MORIS_LIB_Z:FILEPATH=" + str(expected) in cache, cache
                assert "MORIS_Lib_Crypto:" not in cache, cache
print("Runtime prefix/cache handling, Intel parity, and target isolation passed.")
