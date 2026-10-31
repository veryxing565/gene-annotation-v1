#!/usr/bin/env python3
"""Build deterministic, allowlisted source archives without network access."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
TOP = {'README.md','LICENSE','THIRD_PARTY_NOTICES.md','VERSION','CHANGELOG.md',
       'CITATION.cff','CONTRIBUTING.md','SECURITY.md','.gitignore','.gitattributes',
       'requirements.txt','environment.yml','run_gene_annotation.py'}
FOLDERS = {'scripts','pipeline','workflow','tests','config','examples','docs','provenance','.github'}


def digest(data): return hashlib.sha256(data).hexdigest()


def sources():
    result = {}
    for p in sorted(ROOT.rglob('*')):
        rel = p.relative_to(ROOT)
        if '__pycache__' in rel.parts: continue
        if p.is_symlink(): raise ValueError('Symlink in source tree: '+str(rel))
        if not p.is_file(): continue
        if len(rel.parts)==1 and rel.name not in TOP: continue
        if len(rel.parts)>1 and rel.parts[0] not in FOLDERS: continue
        if p.suffix in {'.pyc','.pyo'} or '.local.' in p.name: continue
        data = p.read_bytes()
        if len(data)>2_000_000 or b'\0' in data: raise ValueError('Unexpected binary/large file: '+str(rel))
        data.decode('utf-8')
        if b'\r' in data: raise ValueError('Release sources must use LF: '+str(rel))
        result[rel.as_posix()] = data
    if not TOP.issubset(result): raise ValueError('Missing sources: '+repr(TOP-set(result)))
    return result


def executable(name): return name.endswith(('.py','.pl','.sh'))


def build(output):
    output=output.resolve()
    if output==ROOT or ROOT in output.parents: raise ValueError('Build outside the source tree')
    version=(ROOT/'VERSION').read_text().strip()
    if not version or any(c not in '0123456789.abcdefghijklmnopqrstuvwxyz-' for c in version): raise ValueError('Unsafe version')
    stem='gene-annotation-v1-'+version
    names=[stem+'.zip',stem+'.tar.gz',stem+'.files.json','SHA256SUMS']
    output.mkdir(parents=True,exist_ok=True)
    if any((output/n).exists() for n in names): raise FileExistsError('Use a new artifact directory')
    files=sources()
    record=dict(version=version,source_only=True,files=[dict(path=n,bytes=len(d),sha256=digest(d)) for n,d in files.items()])
    files['MANIFEST.sha256']=''.join(digest(d)+'  '+n+'\n' for n,d in files.items()).encode()
    with zipfile.ZipFile(output/names[0],'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name,data in files.items():
            item=zipfile.ZipInfo(stem+'/'+name,date_time=(2026,10,31,0,0,0)); item.create_system=3
            item.external_attr=(0o100755 if executable(name) else 0o100644)<<16
            item.compress_type=zipfile.ZIP_DEFLATED; archive.writestr(item,data)
    with (output/names[1]).open('xb') as raw:
        with gzip.GzipFile(fileobj=raw,filename='',mode='wb',mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w') as archive:
                for name,data in files.items():
                    item=tarfile.TarInfo(stem+'/'+name); item.size=len(data); item.mtime=0
                    item.mode=0o755 if executable(name) else 0o644; item.uid=item.gid=0
                    item.uname=item.gname=''; archive.addfile(item,io.BytesIO(data))
    (output/names[2]).write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    (output/names[3]).write_text(''.join(digest((output/n).read_bytes())+'  '+n+'\n' for n in names[:3]),encoding='ascii')
    print(json.dumps(dict(version=version,source_files=len(record['files']),artifacts=[str(output/n) for n in names]),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
