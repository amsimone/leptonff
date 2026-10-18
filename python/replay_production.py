#!/usr/bin/env python3
"""Replay a completed production manifest after verifying its archive checksum."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import configparser


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument('manifest', type=Path)
    p.add_argument('--engine', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a, options = p.parse_known_args()
    manifest = a.manifest.resolve()
    meta = json.loads(manifest.read_text())
    archive = meta.get('production_archive', {})
    if meta['status'] != 'completed' or not archive.get('sha256'):
        raise RuntimeError('A completed, checksummed production archive is required')
    if not archive.get('full_born_angles') or archive.get('contains_qed'):
        raise RuntimeError('Archive must cover full Born angles and contain no QED FSR')
    records = (manifest.parent/archive['file']).resolve()
    if digest(records) != archive['sha256']:
        raise RuntimeError('Production archive checksum mismatch')
    if digest(manifest.parent/'input.ini') != meta['config_sha256']:
        raise RuntimeError('Production input configuration has changed')
    config=configparser.ConfigParser(interpolation=None)
    config.read_string('[version]\n'+(manifest.parent/'input.ini').read_text())
    if config.get('cuts','makecuts').strip().lower()!='.false.':
        raise RuntimeError('Born-lepton cuts were enabled in the production input')
    if meta['process']=='z':
        requested={'--mass-min':60.,'--mass-max':120.}
        for i,option in enumerate(options[:-1]):
            if option in requested:requested[option]=float(options[i+1])
        if (requested['--mass-min']<config.getfloat('masscuts','m34min') or
                requested['--mass-max']>config.getfloat('masscuts','m34max')):
            raise RuntimeError('Requested observed mass window exceeds the generated parent-mass coverage')
    engine = a.engine
    if engine is None:
        candidates = [Path(__file__).resolve().with_name('leptonff_replay'),
                      Path(__file__).resolve().parents[1]/'build-host/leptonff_replay']
        engine = next((path for path in candidates if path.is_file()), None)
        if engine is None and shutil.which('leptonff_replay'):
            engine = Path(shutil.which('leptonff_replay'))
    if engine is None:
        raise RuntimeError('Install leptonff_replay or pass --engine')
    engine_hash=digest(engine)
    manifest_hash=digest(manifest)
    record_stat=records.stat()
    command = [str(engine.resolve()), str(records), '--output', str(a.output.resolve()), *options]
    subprocess.run(command, check=True)
    after=records.stat()
    if (digest(engine)!=engine_hash or digest(manifest)!=manifest_hash or
            (after.st_ino,after.st_size,after.st_mtime_ns)!=(record_stat.st_ino,record_stat.st_size,record_stat.st_mtime_ns)):
        raise RuntimeError('Replay executable or production input changed during evaluation')
    result = json.loads(a.output.read_text())
    if (result['physics_id'] != archive['physics_id'] or result['run_id'] != archive['run_id']
            or result['random_seed'] != archive['random_seed'] or result['process'] != meta['process']):
        raise RuntimeError('Replayed production metadata differs from the manifest')
    result['provenance'] = dict(manifest=str(manifest), manifest_sha256=manifest_hash,
                                archive_sha256=archive['sha256'], engine_sha256=engine_hash,
                                command=command, host_integration_warnings=meta.get('integration_warnings', []))
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
