"""Integration tests: python3 -m unittest discover -s share/scripts/tests -v."""

import concurrent.futures
import ctypes
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[3]
SCRIPT = REPO / "share/scripts/create_shared_object.sh"


class InputBuildTest(unittest.TestCase):
    """Check isolated input builds against a read-only MORIS installation."""

    def setUp(self):
        # Spaces exercise shell quoting throughout both compilation scripts.
        self.temp = tempfile.TemporaryDirectory(prefix="moris input test ")
        self.root = Path(self.temp.name)
        self.moris = self.root / "read only moris"
        self.source = self.moris / "source"
        self.build = self.moris / "build"
        self.work = self.root / "work dir"
        self.scratch = self.root / "scratch dir"
        for path in (self.source, self.work, self.scratch):
            path.mkdir(parents=True)
        # The wrapper locates its delegated builder under MORISROOT. Include
        # that script in the fixture so concurrent tests exercise the public API.
        scripts = self.moris / "share/scripts"
        scripts.mkdir(parents=True)
        shutil.copyfile(str(SCRIPT), str(scripts / SCRIPT.name))
        # An old lock in its original location must neither block new builds
        # nor require deletion from the read-only installation.
        legacy_lock = self.moris / "projects/mains/input_file.locked"
        legacy_lock.parent.mkdir(parents=True)
        legacy_lock.write_text("stale lock from the previous input builder\n")
        (self.source / "placeholder.cpp").write_text("// original placeholder\n")
        (self.source / "fixture.hpp").write_text("constexpr int fixture_value = 40;\n")
        # Exercise the same target exporter used by projects/mains, without
        # requiring the full MORIS dependency stack.
        (self.source / "CMakeLists.txt").write_text(f'''
cmake_minimum_required(VERSION 3.17)
project(InputFixture LANGUAGES CXX)
set(MORIS_CMAKE_DIR "{REPO}/share/cmake")
include("${{MORIS_CMAKE_DIR}}/utilities/dynamic_link_input.cmake")
set(CMAKE_CXX_FLAGS "${{CMAKE_CXX_FLAGS}} -DFLAG_VALUE=1")
add_compile_definitions(DIRECTORY_VALUE=1)
include_directories("${{CMAKE_CURRENT_SOURCE_DIR}}")
add_library(shared_object_file SHARED placeholder.cpp)
target_compile_features(shared_object_file PRIVATE cxx_std_17)
target_compile_options(shared_object_file PRIVATE -Wall -Werror)
target_compile_definitions(shared_object_file PRIVATE TARGET_VALUE=1)
generate_input_build_project(shared_object_file)
''')
        self.run_ok(["cmake", "-S", str(self.source), "-B", str(self.build),
                     "-DCMAKE_BUILD_TYPE=Release"])
        self.env = dict(os.environ, MORISROOT=str(self.moris), TMPDIR=str(self.scratch))
        self.env.pop("MORIS_INPUT_BUILD_CONFIG", None)
        self.original = self.snapshot()
        # Prevent runtime writes, and also compare contents/mtimes afterward
        # so the test detects modifications even when run with elevated access.
        for path in sorted(self.moris.rglob("*"), reverse=True):
            path.chmod(0o555 if path.is_dir() or path.suffix == ".sh" else 0o444)
        self.moris.chmod(0o555)
        self.write_input("input file")

    def tearDown(self):
        # Restore directory permissions so TemporaryDirectory can remove it.
        self.moris.chmod(0o755)
        for path in self.moris.rglob("*"):
            if path.is_dir():
                path.chmod(0o755)
        self.temp.cleanup()

    def snapshot(self):
        """Detect both replaced files and otherwise invisible source touches."""
        return {str(p.relative_to(self.moris)): (p.read_bytes(), p.stat().st_mtime_ns)
                for p in self.moris.rglob("*") if p.is_file()}

    def run_ok(self, args, **kwargs):
        # A lock regression must fail with a timeout instead of hanging CI.
        result = subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=30, **kwargs)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result

    def write_input(self, name):
        # Compilation and the exported answer jointly check inherited flags,
        # definitions, include paths, configuration, and language features.
        (self.work / f"{name}.cpp").write_text('''
#include "fixture.hpp"
#ifndef NDEBUG
#error Release configuration was not preserved
#endif
extern "C" int answer() {
    if constexpr (TARGET_VALUE == DIRECTORY_VALUE) {
        return fixture_value + FLAG_VALUE + TARGET_VALUE;
    }
    return 0;
}
''')

    def command(self, name="input file", *extra):
        """Use a workdir distinct from the caller's current directory."""
        return ["sh", str(SCRIPT), str(self.work), "build", name, *extra]

    def assert_unchanged(self):
        """The source, original build, and legacy lock must remain untouched."""
        self.assertEqual(self.original, self.snapshot())

    def test_read_only_build_and_load(self):
        """A read-only source produces a loadable library outside MORIS."""
        (self.work / "input file.cpp").chmod(0o444)
        self.run_ok(self.command(), env=self.env, cwd=self.root)
        lib = ctypes.CDLL(str(self.work / "input file.so"))
        self.assertEqual(lib.answer(), 42)
        self.assertEqual(list(self.scratch.iterdir()), [])
        self.assert_unchanged()

    def test_concurrent_builds(self):
        """Public wrapper calls finish concurrently despite a stale lock."""
        self.write_input("second")
        wrapper = REPO / "share/scripts/moris_create_shared_object.sh"
        env = dict(self.env, MORISBUILDOPT="build", MORISBUILDDBG="build")
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(self.run_ok, ["sh", str(wrapper), "opt", name],
                                   env=env, cwd=str(self.work))
                       for name in ("input file", "second")]
            for future in futures:
                future.result()
        for name in ("input file", "second"):
            self.assertTrue((self.work / f"{name}.so").is_file())
        self.assertEqual(list(self.scratch.iterdir()), [])
        self.assert_unchanged()

    def test_compiler_failure_cleans_scratch(self):
        """Failure propagates without publishing a library or leaving scratch."""
        (self.work / "broken.cpp").write_text("#error intentional failure\n")
        result = subprocess.run(self.command("broken"), env=self.env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.work / "broken.so").exists())
        self.assertEqual(list(self.scratch.iterdir()), [])
        self.assert_unchanged()

    def test_debug_option_retains_private_build(self):
        """The legacy dbg flag retains private artifacts for inspection."""
        self.run_ok(self.command("input file", "dbg"), env=self.env)
        retained = list(self.scratch.iterdir())
        self.assertEqual(len(retained), 1)
        self.assertTrue((retained[0] / "build/CMakeCache.txt").is_file())
        self.assert_unchanged()

    def test_missing_input_fails(self):
        """Reject invalid input before creating a temporary build."""
        result = subprocess.run(self.command("missing"), env=self.env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(list(self.scratch.iterdir()), [])
        self.assert_unchanged()

    def test_empty_build_configuration(self):
        """A parent with no build type still exports a usable input project."""
        build = self.root / "unconfigured build"
        self.run_ok(["cmake", "-S", str(self.source), "-B", str(build)])
        (self.work / "simple.cpp").write_text('extern "C" int answer() { return 42; }\n')
        self.run_ok(["sh", str(SCRIPT), str(self.work), str(build), "simple"],
                    env=self.env)
        self.assertTrue((self.work / "simple.so").exists())
        self.assert_unchanged()

    @unittest.skipUnless(shutil.which("ninja"), "Ninja is unavailable")
    def test_multi_configuration_debug(self):
        """A selected multi-configuration export preserves Debug flags."""
        build = self.root / "multi config build"
        self.run_ok(["cmake", "-S", str(self.source), "-B", str(build),
                     "-G", "Ninja Multi-Config"])
        (self.work / "simple.cpp").write_text('''
#ifdef NDEBUG
#error Expected Debug configuration
#endif
extern "C" int answer() { return 42; }
''')
        self.run_ok(["sh", str(SCRIPT), str(self.work), str(build), "simple"],
                    env=dict(self.env, MORIS_INPUT_BUILD_CONFIG="Debug"))
        self.assertTrue((self.work / "simple.so").exists())
        self.assert_unchanged()


if __name__ == "__main__":
    unittest.main()
