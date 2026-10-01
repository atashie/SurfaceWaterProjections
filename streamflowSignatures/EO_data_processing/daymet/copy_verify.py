"""Copy the run folder's deliverables to the exFAT drive and prove every byte arrived.

The drive has silently truncated files before (mtime preserved), so each file is
md5-hashed locally, copied, flushed, re-read FROM THE DRIVE and hashed again; any
mismatch or size difference is fatal. MANIFEST_md5.txt on the drive holds one line per
file ("md5  size  relpath"); re-copying a file REPLACES its line (it used to append a
second one). The manifest is rewritten through a temporary file and an fsync.
Usage: python copy_verify.py <src_root> <dst_root> <relpath> [<relpath> ...]
"""
import fcntl
import hashlib
import os
import shutil
import sys


def md5(path, block=1 << 24, nocache=False):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        if nocache:                      # read from the device, not the page cache (macOS)
            fcntl.fcntl(fh.fileno(), fcntl.F_NOCACHE, 1)
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def main():
    src, dst, rels = sys.argv[1], sys.argv[2], sys.argv[3:]
    files = []
    for r in rels:
        p = os.path.join(src, r)
        if os.path.isdir(p):
            for root, _, names in os.walk(p):
                files += [os.path.relpath(os.path.join(root, n), src) for n in sorted(names)
                          if not n.startswith(".")]
        else:
            files.append(r)
    lines, bad, total = [], 0, 0
    for rel in files:
        s, d = os.path.join(src, rel), os.path.join(dst, rel)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        hs = md5(s)
        shutil.copyfile(s, d)
        with open(d, "rb+") as fh:
            os.fsync(fh.fileno())
        hd, ss, sd = md5(d, nocache=True), os.path.getsize(s), os.path.getsize(d)
        ok = hs == hd and ss == sd
        bad += not ok
        total += ss
        lines.append(f"{hs}  {ss:>14}  {rel}")
        if not ok:
            print(f"MISMATCH {rel}: src {hs} {ss} B, drive {hd} {sd} B")
    man = os.path.join(dst, "MANIFEST_md5.txt")
    entries = {}                         # relpath -> line; insertion order kept, re-copies replace
    if os.path.exists(man):
        for ln in open(man):
            f = ln.rstrip("\n").split(None, 2)
            if len(f) == 3:
                entries[f[2]] = ln.rstrip("\n")
    for ln in lines:
        entries[ln.split(None, 2)[2]] = ln
    with open(man + ".tmp", "w") as fh:
        fh.write("\n".join(entries.values()) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(man + ".tmp", man)
    print(f"copied {len(files)} files, {total / 1e9:.2f} GB; md5 + size verified on the drive: "
          f"{len(files) - bad} OK, {bad} mismatches; manifest {man}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
