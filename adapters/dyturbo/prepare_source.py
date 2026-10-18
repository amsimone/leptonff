#!/usr/bin/env python3
"""Prepare the requested DYTurbo 1.4.2 archive in an isolated source directory."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile

VERSION = '1.4.2'
URL = 'https://dyturbo.hepforge.org/downloads/dyturbo-1.4.2.tar.gz'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def prepare(archive, directory):
    archive, directory = Path(archive).resolve(), Path(directory).resolve()
    if not archive.is_file():
        raise FileNotFoundError(f'Download {URL} to {archive} first')
    if directory.exists():
        raise RuntimeError('Refusing to overwrite an existing source directory')
    rootname = 'dyturbo-'+VERSION
    with tarfile.open(archive, 'r:gz') as file:
        members = file.getmembers()
        if len(members)>50000 or sum(m.size for m in members)>2*1024**3:
            raise RuntimeError('Unexpectedly large source archive')
        for member in members:
            path = PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0]!=rootname:
                raise RuntimeError('Unsafe or unexpected archive path: '+member.name)
            if not (member.isdir() or member.isfile() or member.issym() or member.islnk()):
                raise RuntimeError('Unsupported special archive member: '+member.name)
        directory.mkdir(parents=True)
        source = directory/rootname
        try:
            for member in members:
                destination = directory/member.name
                if not destination.resolve().is_relative_to(source):
                    raise RuntimeError('Archive destination leaves the source directory')
                if member.issym() or member.islnk():
                    if PurePosixPath(member.linkname).is_absolute():
                        raise RuntimeError('Absolute archive link')
                    target = destination.parent/member.linkname if member.issym() else directory/member.linkname
                    if not target.resolve().is_relative_to(source):
                        raise RuntimeError('Archive link leaves the source directory')
                file.extract(member, directory)
            version_files = [p for p in [source/'configure.ac',source/'configure'] if p.is_file()]
            if not version_files or not any(re.search(r'(?<!\d)1\.4\.2(?!\d)', p.read_text(errors='replace')) for p in version_files):
                raise RuntimeError('DYTurbo 1.4.2 version not confirmed in configure sources')
            candidates = []
            signature = re.compile(r'(theta[01]|costh|fillhisto|resint|ctint|namespace\s+DYTurbo|class\s+Kinematics)', re.I)
            for folder in ['src','histo','resum','counterterm','integrands','phasespace','born','user']:
                for path in sorted((source/folder).rglob('*')):
                    if path.suffix.lower() not in ['.h','.hpp','.c','.cc','.cpp','.cxx','.f','.f90'] or not path.is_file():
                        continue
                    matches = [(i,line.strip()) for i,line in enumerate(path.read_text(errors='replace').splitlines(),1) if signature.search(line)]
                    if matches:
                        candidates.append(dict(file=str(path.relative_to(source)),sha256=digest(path),matches=matches[:12]))
            report = dict(version=VERSION,upstream=URL,archive_sha256=digest(archive),source_directory=str(source),
                          status='prepared; native coupling still requires source-level implementation and validation',
                          interface_candidates=candidates)
            (directory/'source_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
            return report
        except Exception:
            shutil.rmtree(directory)
            raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    p.add_argument('--output',type=Path,default=Path('external/dyturbo-1.4.2-work'))
    a = p.parse_args()
    print(json.dumps(prepare(a.archive,a.output),indent=2))


if __name__=='__main__':
    main()
