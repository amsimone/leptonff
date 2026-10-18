"""Check isolated source preparation; these are synthetic archive fixtures."""
import io
from pathlib import Path
import sys
import tarfile
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'adapters/dyturbo'))
from prepare_source import prepare

with tempfile.TemporaryDirectory() as directory:
    base=Path(directory)
    def archive(name,extra=None):
        path=base/name
        with tarfile.open(path,'w:gz') as file:
            value=b'AC_INIT([DYTurbo],[1.4.2])\n'
            entry=tarfile.TarInfo('dyturbo-1.4.2/configure.ac');entry.size=len(value);file.addfile(entry,io.BytesIO(value))
            if extra:file.addfile(extra,io.BytesIO(b''))
        return path
    good=archive('good.tar.gz')
    result=prepare(good,base/'prepared')
    assert result['version']=='1.4.2' and (base/'prepared/source_manifest.json').is_file()
    for index,member in enumerate([tarfile.TarInfo('../outside'),tarfile.TarInfo('/absolute')]):
        bad=archive(f'bad{index}.tar.gz',member)
        try:prepare(bad,base/f'bad{index}')
        except RuntimeError:pass
        else:raise AssertionError('Unsafe path accepted')
    link=tarfile.TarInfo('dyturbo-1.4.2/link');link.type=tarfile.SYMTYPE;link.linkname='../../outside'
    bad=archive('link.tar.gz',link)
    try:prepare(bad,base/'linked')
    except RuntimeError:pass
    else:raise AssertionError('External link accepted')
    assert not (base/'linked').exists()
print('Versioned DYTurbo preparation and extraction checks passed on synthetic fixtures')
