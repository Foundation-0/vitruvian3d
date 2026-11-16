"""Regression test for the bare-assert verifier fault.

A corrupted certificate must fail bare AND under `python3 -O`, and the
certificate file must be left unchanged. A valid certificate must still
verify under both modes. Run from anywhere:

    python3 research/tests/test_assert_fix.py
    python3 -O research/tests/test_assert_fix.py

Standard library only. This file itself contains no bare `assert`.
"""

import ast
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESEARCH = ROOT / "research"


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def run(args, cwd):
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=1200)


def test_no_bare_asserts():
    files = sorted(RESEARCH.rglob("*.py"))
    check(len(files) > 0, "no .py files found under research/")
    remaining = {}
    for path in files:
        tree = ast.parse(path.read_text())
        count = sum(1 for node in ast.walk(tree) if isinstance(node, ast.Assert))
        if count:
            remaining[str(path.relative_to(ROOT))] = count
    check(not remaining, "bare asserts remain: %s" % (remaining,))
    print("PASS: no bare asserts in %d research .py files" % len(files))


def _copy_passage(tmp):
    shutil.copy2(RESEARCH / "passage_verify.py", tmp / "passage_verify.py")
    shutil.copy2(RESEARCH / "passage_certificate.json", tmp / "passage_certificate.json")


def _corrupt_interval_digit(cert_path):
    text = cert_path.read_text()
    old = '"363969/1000000"'
    check(old in text, "expected interval string not in certificate")
    cert_path.write_text(text.replace(old, '"463969/1000000"', 1))


def test_passage_verify_corrupted():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        _copy_passage(tmp)
        _corrupt_interval_digit(tmp / "passage_certificate.json")
        before = (tmp / "passage_certificate.json").read_bytes()
        for mode in ([], ["-O"]):
            result = run([sys.executable, *mode, "passage_verify.py"], cwd=tmp)
            label = "python3 -O" if mode else "python3"
            check(result.returncode != 0,
                  "corrupted certificate passed under %s" % label)
            check((tmp / "passage_certificate.json").read_bytes() == before,
                  "verifier rewrote the certificate under %s" % label)
            print("PASS: corrupted certificate rejected under %s (exit %d)"
                  % (label, result.returncode))


def test_passage_verify_valid():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        _copy_passage(tmp)
        expected = json.loads((RESEARCH / "passage_verify.json").read_text())
        for mode in ([], ["-O"]):
            result = run([sys.executable, *mode, "passage_verify.py"], cwd=tmp)
            label = "python3 -O" if mode else "python3"
            check(result.returncode == 0,
                  "valid certificate failed under %s: %s" % (label, result.stderr[-500:]))
            check(json.loads(result.stdout) == expected,
                  "valid certificate output differs from committed receipt under %s" % label)
            print("PASS: valid certificate verifies under %s" % label)


def test_passage_verify_missing_and_empty():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        shutil.copy2(RESEARCH / "passage_verify.py", tmp / "passage_verify.py")
        result = run([sys.executable, "passage_verify.py"], cwd=tmp)
        check(result.returncode != 0, "missing certificate passed")
        print("PASS: missing certificate rejected (exit %d)" % result.returncode)
        (tmp / "passage_certificate.json").write_text("")
        for mode in ([], ["-O"]):
            result = run([sys.executable, *mode, "passage_verify.py"], cwd=tmp)
            label = "python3 -O" if mode else "python3"
            check(result.returncode != 0, "empty certificate passed under %s" % label)
            print("PASS: empty certificate rejected under %s (exit %d)"
                  % (label, result.returncode))


def test_passage_window_verify_corrupted():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        shutil.copytree(RESEARCH, tmp / "research")
        cert = tmp / "research" / "passage_certificate.json"
        _corrupt_interval_digit(cert)
        for mode in ([], ["-O"]):
            result = run([sys.executable, *mode, "passage_window_verify.py"],
                         cwd=tmp / "research")
            label = "python3 -O" if mode else "python3"
            check(result.returncode != 0,
                  "corrupted certificate passed passage_window_verify under %s" % label)
            print("PASS: passage_window_verify rejects corrupted certificate under %s (exit %d)"
                  % (label, result.returncode))


def test_transport_check_valid():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        shutil.copytree(RESEARCH, tmp / "research")
        check_py = tmp / "research" / "transport-conductivity-20261112" / "check.py"
        result = run([sys.executable, str(check_py),
                      "--source-root", str(tmp),
                      "--output", str(tmp / "out.json")], cwd=tmp)
        check(result.returncode == 0,
              "transport check failed on valid source: %s" % result.stderr[-500:])
        produced = json.loads((tmp / "out.json").read_text())
        committed = json.loads((RESEARCH / "transport-conductivity-20261112"
                                / "certificate.json").read_text())
        produced.pop("python", None)
        committed.pop("python", None)
        check(produced == committed,
              "transport check output differs from committed certificate.json")
        print("PASS: transport check passes on valid source")


def test_transport_check_corrupted():
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        shutil.copytree(RESEARCH, tmp / "research")
        cert = tmp / "research" / "passage_certificate.json"
        _corrupt_interval_digit(cert)
        check_py = tmp / "research" / "transport-conductivity-20261112" / "check.py"
        for mode in ([], ["-O"]):
            result = run([sys.executable, *mode, str(check_py),
                          "--source-root", str(tmp),
                          "--output", str(tmp / "out.json")], cwd=tmp)
            label = "python3 -O" if mode else "python3"
            check(result.returncode != 0,
                  "corrupted certificate passed transport check under %s" % label)
            print("PASS: transport check rejects corrupted certificate under %s (exit %d)"
                  % (label, result.returncode))


def main():
    test_no_bare_asserts()
    test_passage_verify_corrupted()
    test_passage_verify_valid()
    test_passage_verify_missing_and_empty()
    test_passage_window_verify_corrupted()
    test_transport_check_valid()
    test_transport_check_corrupted()
    print("PASS: all assert-fix regression tests")


if __name__ == "__main__":
    main()
