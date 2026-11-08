"""Portable setup tests use temporary folders and mocks; never install packages."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import setup_environment as setup


class SetupEnvironmentTests(unittest.TestCase):
    def setUp(self):
        project_root = Path(__file__).resolve().parent
        self.temporary = tempfile.TemporaryDirectory(prefix="trip setup test ", dir=project_root)
        self.root = Path(self.temporary.name).resolve()
        self.assertEqual(self.root.parent, project_root)
        self.addCleanup(self.temporary.cleanup)
        (self.root / "requirements.txt").write_text("numpy>=1.24\n", encoding="utf-8")
        self.environment = self.root / ".venv"
        self.python = setup.environment_python(self.environment)

    def create_fake_environment(self, directory=None):
        self.python.parent.mkdir(parents=True, exist_ok=True)
        self.python.write_text("mock interpreter", encoding="utf-8")

    def environment_check_result(self, *, version=(3, 12, 0), prefix=None, base_prefix=None):
        return subprocess.CompletedProcess([], 0, stdout=json.dumps({
            "version": version,
            "prefix": str(prefix or self.environment),
            "base_prefix": str(base_prefix or self.root / "base-python"),
        }), stderr="")

    def run_setup(self):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return setup.setup_environment(self.root)

    def test_creates_project_local_environment_and_installs_requirements(self):
        with patch.object(setup.venv, "EnvBuilder") as builder, \
                patch.object(setup.subprocess, "run") as run:
            builder.return_value.create.side_effect = self.create_fake_environment
            run.side_effect = [self.environment_check_result(), subprocess.CompletedProcess([], 0)]
            self.assertEqual(self.run_setup(), 0)
            builder.assert_called_once_with(with_pip=True)
            builder.return_value.create.assert_called_once_with(self.environment)
            install = run.call_args_list[1]
            self.assertEqual(install.args[0][0], str(self.python))
            self.assertEqual(install.args[0][-2:], ["-r", str(self.root / "requirements.txt")])
            self.assertEqual(install.kwargs["cwd"], self.root)

    def test_reuses_existing_environment_without_recreating_it(self):
        self.create_fake_environment()
        sentinel = self.environment / "keep.txt"
        sentinel.write_text("keep", encoding="utf-8")
        with patch.object(setup.venv, "EnvBuilder") as builder, \
                patch.object(setup.subprocess, "run", side_effect=[
                    self.environment_check_result(), subprocess.CompletedProcess([], 0)]):
            self.assertEqual(self.run_setup(), 0)
            builder.assert_not_called()
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep")

    def test_missing_requirements_does_not_create_environment(self):
        (self.root / "requirements.txt").unlink()
        with patch.object(setup.venv, "EnvBuilder") as builder:
            self.assertEqual(self.run_setup(), 2)
            builder.assert_not_called()
        self.assertFalse(self.environment.exists())

    def test_incomplete_existing_directory_is_preserved(self):
        self.environment.mkdir()
        sentinel = self.environment / "keep.txt"
        sentinel.write_text("keep", encoding="utf-8")
        with patch.object(setup.venv, "EnvBuilder") as builder, \
                patch.object(setup.subprocess, "run") as run:
            self.assertEqual(self.run_setup(), 2)
            builder.assert_not_called()
            run.assert_not_called()
        self.assertTrue(sentinel.is_file())

    def test_old_environment_python_is_rejected_before_install(self):
        self.create_fake_environment()
        with patch.object(setup.subprocess, "run", return_value=self.environment_check_result(
                version=(3, 9, 9))) as run:
            self.assertEqual(self.run_setup(), 2)
            self.assertEqual(run.call_count, 1)

    def test_system_interpreter_cannot_receive_dependencies(self):
        self.create_fake_environment()
        with patch.object(setup.subprocess, "run", return_value=self.environment_check_result(
                prefix=self.root / "base-python")) as run:
            self.assertEqual(self.run_setup(), 2)
            self.assertEqual(run.call_count, 1)

    def test_failed_install_keeps_environment_for_retry(self):
        self.create_fake_environment()
        with patch.object(setup.subprocess, "run", side_effect=[
                self.environment_check_result(), subprocess.CompletedProcess([], 1)]):
            self.assertEqual(self.run_setup(), 1)
        self.assertTrue(self.python.is_file())

    def test_failed_creation_does_not_delete_partial_files(self):
        def fail_creation(directory):
            self.create_fake_environment(directory)
            raise OSError("mock creation failure")
        with patch.object(setup.venv, "EnvBuilder") as builder, \
                patch.object(setup.subprocess, "run") as run:
            builder.return_value.create.side_effect = fail_creation
            self.assertEqual(self.run_setup(), 1)
            run.assert_not_called()
        self.assertTrue(self.python.is_file())


if __name__ == "__main__":
    unittest.main()
