"""Check that passage_verify.py rejects a corrupted certificate in both run modes.

A copy of the verifier and certificate is made in a temporary folder. One half-angle
interval digit is changed in the copy. The copy must fail with a non-zero exit under
plain python and under python -O. The unchanged copy must pass in both modes.
The repository files and receipts are not written.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
GOOD = '"90993/250000"'
BAD = '"99993/250000"'


def run(folder, optimize):
    command = [sys.executable] + (['-O'] if optimize else []) + [str(folder / 'passage_verify.py')]
    return subprocess.run(command, cwd=folder, capture_output=True, text=True)


def main():
    certificate = (HERE / 'passage_certificate.json').read_text()
    if certificate.count(GOOD) != 1:
        raise SystemExit('Unexpected certificate: the digit to corrupt was not found once.')
    with tempfile.TemporaryDirectory() as tmp:
        good, bad = Path(tmp) / 'good', Path(tmp) / 'bad'
        for folder, text in ((good, certificate), (bad, certificate.replace(GOOD, BAD))):
            folder.mkdir()
            shutil.copy(HERE / 'passage_verify.py', folder)
            (folder / 'passage_certificate.json').write_text(text)
        for optimize in (False, True):
            mode = 'python -O' if optimize else 'python'
            accepted = run(good, optimize)
            if accepted.returncode != 0:
                raise SystemExit(f'{mode}: the unchanged certificate was rejected.\n{accepted.stderr}')
            rejected = run(bad, optimize)
            if rejected.returncode == 0:
                raise SystemExit(f'{mode}: the corrupted certificate was accepted.')
            if 'VerificationError' not in rejected.stderr:
                raise SystemExit(f'{mode}: the corrupted certificate failed for another reason.\n{rejected.stderr}')
            print(f'{mode}: unchanged certificate accepted, corrupted certificate rejected')
    print('PASS: corrupted passage certificate rejected with and without -O.')


if __name__ == '__main__':
    main()
