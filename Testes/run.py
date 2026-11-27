"""Run offline regression tests or compare the native Fortran and Python programs."""
import argparse
from pathlib import Path
import subprocess
import sys


TESTS = Path(__file__).resolve().parent
TOOLS = TESTS.parent / "chatgpt-python-unfold" / "tools"


def offline(output):
    commands = [
        ["-m", "unittest", "discover", "-v", "-s", str(TESTS)],
        [str(TOOLS / "check_graphene.py"), "--output", str(output)],
        [str(TOOLS / "check_projection.py"), str(output / "unfolded.csv"),
         "--output", str(output / "projection.json")],
    ]
    for command in commands:
        status = subprocess.run([sys.executable, *command]).returncode
        if status:
            return status
    return 0


def native(output, dftb, unfold, mio):
    failures = []
    for case in ("graphene", "diamond"):
        command = [
            sys.executable, str(TOOLS / "check_native.py"), case,
            "--dftb", str(dftb), "--unfold", str(unfold),
            "--mio", str(mio), "--output", str(output),
        ]
        if subprocess.run(command).returncode:
            failures.append(case)
    if failures:
        print("Failed native comparisons: " + ", ".join(failures), file=sys.stderr)
    return int(bool(failures))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("offline", "native"), nargs="?", default="offline")
    parser.add_argument("--output", type=Path, help="directory for CSV and validation reports")
    group = parser.add_argument_group("native comparison")
    group.add_argument("--dftb", type=Path, help="path to the DFTB+ 26.1 executable")
    group.add_argument("--unfold", type=Path, help="path to the compiled Fortran unfold executable")
    group.add_argument("--mio", type=Path, help="mio-1-1 directory with skfiles/ and extras/")
    args = parser.parse_args()
    if args.mode == "offline":
        if any(path is not None for path in (args.dftb, args.unfold, args.mio)):
            parser.error("--dftb, --unfold and --mio apply to native mode")
        return offline((args.output or Path("results/testes")).resolve())

    if not all((args.dftb, args.unfold, args.mio)):
        parser.error("native mode requires --dftb, --unfold and --mio")
    for executable in (args.dftb, args.unfold):
        if not executable.is_file():
            parser.error(f"executable file not found: {executable}")
    if not args.mio.is_dir():
        parser.error(f"mio-1-1 directory not found: {args.mio}")
    return native((args.output or Path("results/native")).resolve(),
                  args.dftb.resolve(), args.unfold.resolve(), args.mio.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
