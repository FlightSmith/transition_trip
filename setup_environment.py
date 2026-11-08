"""Create/reuse this project's virtual environment and install requirements.

Run with Python 3.10 or newer. The environment is always located beside this
script, so setup also works when invoked from another working directory.
Existing environments are never deleted or recreated automatically.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
import venv


MINIMUM_PYTHON = (3, 10)


def environment_python(environment_dir: Path) -> Path:
    """Return the interpreter location appropriate to the current platform."""
    return environment_dir / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")


def inspect_environment(python_path: Path, environment_dir: Path) -> None:
    """Verify that pip will run in the intended project environment."""
    code = (
        "import json,sys; print(json.dumps({'version':list(sys.version_info[:3]),"
        "'prefix':sys.prefix,'base_prefix':sys.base_prefix}))"
    )
    try:
        result = subprocess.run([str(python_path), "-I", "-c", code],
                                capture_output=True, text=True, check=False)
    except OSError as exc:
        raise ValueError(f"The existing environment's Python could not start: {exc}") from exc
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise ValueError(f"The environment's Python did not start successfully. {detail}")
    try:
        info = json.loads(result.stdout)
        version = tuple(int(part) for part in info["version"][:2])
        prefix = Path(info["prefix"]).resolve()
        base_prefix = Path(info["base_prefix"]).resolve()
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError("Could not verify the environment's Python version and location.") from exc
    if version < MINIMUM_PYTHON:
        raise ValueError("The existing environment requires Python 3.10 or newer.")
    if prefix != environment_dir.resolve() or prefix == base_prefix:
        raise ValueError("The interpreter is not using this project's .venv environment.")


def command_text(arguments: list[str]) -> str:
    return subprocess.list2cmdline(arguments) if sys.platform == "win32" else shlex.join(arguments)


def setup_environment(project_root: Path) -> int:
    """Install the declared dependencies; return 0 on success, 1/2 on failure."""
    if sys.version_info < MINIMUM_PYTHON:
        print("Setup needs Python 3.10 or newer. Install a supported Python and try again.",
              file=sys.stderr)
        return 2
    project_root = Path(project_root).resolve()
    requirements = project_root / "requirements.txt"
    environment_dir = project_root / ".venv"
    python_path = environment_python(environment_dir)
    if not requirements.is_file():
        print(f"Cannot find {requirements}. Keep the project files together.", file=sys.stderr)
        return 2
    if environment_dir.is_symlink():
        print("The project .venv is a symbolic link. Use a normal project-local environment; "
              "setup has not changed it.", file=sys.stderr)
        return 2
    if environment_dir.exists():
        if not environment_dir.is_dir() or not python_path.is_file():
            print(f"{environment_dir} already exists but is not a complete Python environment. "
                  "Inspect it and rename it if appropriate, then run setup again. "
                  "Nothing has been deleted.", file=sys.stderr)
            return 2
        print(f"Reusing {environment_dir}", flush=True)
    else:
        print(f"Creating {environment_dir}", flush=True)
        try:
            venv.EnvBuilder(with_pip=True).create(environment_dir)
        except (OSError, subprocess.SubprocessError) as exc:
            print(f"Could not create the environment: {exc}\n"
                  "Check that this Python installation includes venv and pip. "
                  "If a partial .venv was created, inspect and rename it before retrying. "
                  "No existing environment has been deleted.", file=sys.stderr)
            return 1
    try:
        inspect_environment(python_path, environment_dir)
    except ValueError as exc:
        print(f"Environment check failed: {exc}\n"
              "Inspect or rename .venv before setting up a replacement. "
              "No existing environment has been deleted.", file=sys.stderr)
        return 2
    print("Installing the project dependencies from requirements.txt...", flush=True)
    command = [str(python_path), "-I", "-m", "pip", "install",
               "--disable-pip-version-check", "-r", str(requirements)]
    try:
        result = subprocess.run(command, cwd=project_root, check=False)
    except OSError as exc:
        print(f"Could not start dependency installation: {exc}", file=sys.stderr)
        return 1
    if result.returncode:
        print("Dependency installation did not complete. Check the connection or the pip "
              "error above, then run setup again. The environment has been kept.", file=sys.stderr)
        return 1
    print("\nSetup complete. Start the tool with:")
    print(command_text([str(python_path), str(project_root / "run_trip.py")]))
    if sys.platform == "win32":
        print("You can also double-click run_windows.cmd.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(argv)
    return setup_environment(Path(__file__).resolve().parent)


if __name__ == "__main__":
    raise SystemExit(main())
