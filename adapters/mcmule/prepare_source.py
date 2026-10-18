#!/usr/bin/env python3
"""Apply one explicit-interface dispatch workaround; never edit amplitudes."""
from pathlib import Path
import hashlib
import subprocess

COMMIT = 'de32e762ca44ab82345d4a511f7b64ab61ba9bd4'
URL = 'https://gitlab.com/mule-tools/mcmule.git'
NEEDLE = '''    case(4)
      call gen(x, vecs(:,1), masses(1), vecs(:,2), masses(2), vecs(:,3), masses(3), vecs(:,4), masses(4), weight)'''
REPLACEMENT = '''    case(4)
      ! Explicit interface avoids a gfortran 16 / arm64 indirect-call ABI failure.
      if (associated(gen, psx2)) then
        call psx2(x, vecs(:,1), masses(1), vecs(:,2), masses(2), vecs(:,3), masses(3), vecs(:,4), masses(4), weight)
      else
        call gen(x, vecs(:,1), masses(1), vecs(:,2), masses(2), vecs(:,3), masses(3), vecs(:,4), masses(4), weight)
      endif'''


def prepare(source):
    source = Path(source)
    if subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip() != COMMIT:
        raise ValueError('Need the pinned McMule v0.8.0 commit')
    path = 'pre/phase_space.f95'
    original = subprocess.check_output(['git', 'show', 'HEAD:' + path], cwd=source)
    if original.decode().count(NEEDLE) != 1:
        raise ValueError('Unexpected upstream phase-space dispatch')
    prepared = original.decode().replace(NEEDLE, REPLACEMENT).encode()
    if (source / path).read_bytes() not in (original, prepared):
        raise ValueError('Unrecognized local phase-space changes')
    changes = subprocess.check_output(['git', 'diff', '--name-only', 'HEAD'], cwd=source, text=True).splitlines()
    if set(changes) - {path}:
        raise ValueError('Other upstream sources have changed')
    (source / path).write_bytes(prepared)
    return dict(path=path, upstream_sha256=hashlib.sha256(original).hexdigest(),
                prepared_sha256=hashlib.sha256(prepared).hexdigest(),
                purpose='Explicit PSX2 dispatch; no amplitude, phase-space formula or subtraction changes')


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    print(json.dumps(prepare(p.parse_args().source), indent=2))
