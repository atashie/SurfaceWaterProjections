"""Copy the run folder's deliverables to the exFAT drive and prove every byte arrived.

The drive has silently truncated files before (mtime preserved). So each file is:
  * md5-hashed locally;
  * copied to `<name>.new` on the drive and flushed (F_FULLFSYNC);
  * re-read FROM THE DRIVE and hashed again.
Only then does the copy replace the drive's file and its manifest line. A mismatch or size
difference leaves the earlier drive copy and its line untouched, and the script exits 1.

MANIFEST_md5.txt on the drive holds one line per file ("md5  size  relpath"). It is
rewritten after each verified file, through a temporary file and an fsync. A file whose
manifest line records OTHER content is refused unless --replace is given, so a truncated
or wrong source cannot quietly overwrite a known-good copy.

The re-read rate is printed. A rate far above the drive's own (the exFAT thumb drive reads
~30-50 MB/s) means the pages came from memory, not the device: F_NOCACHE does not evict
pages that are already resident. Such a file is noted; re-hash it after a remount to be sure.

Usage: python copy_verify.py [--replace] <src_root> <dst_root> <relpath> [<relpath> ...]
Hidden files and folders are skipped.
"""
import argparse
import fcntl
import hashlib
import os
import shutil
import sys
import time


def md5(path, block=1 << 24, nocache=False):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        if nocache:                      # read from the device, not the page cache (macOS)
            fcntl.fcntl(fh.fileno(), fcntl.F_NOCACHE, 1)
        for b in iter(lambda: fh.read(block), b""):
            h.update(b)
    return h.hexdigest()


def full_fsync(fh):
    os.fsync(fh.fileno())
    if hasattr(fcntl, "F_FULLFSYNC"):    # macOS: also flush the device's own write cache
        fcntl.fcntl(fh.fileno(), fcntl.F_FULLFSYNC)


def read_manifest(man):
    entries = {}                         # relpath -> line; insertion order kept, re-copies replace
    if os.path.exists(man):
        for ln in open(man):
            f = ln.rstrip("\n").split(None, 2)
            if len(f) == 3:
                entries[os.path.normpath(f[2])] = ln.rstrip("\n")
    return entries


def write_manifest(man, entries):
    with open(man + ".tmp", "w") as fh:
        fh.write("\n".join(entries.values()) + "\n")
        full_fsync(fh)
    os.replace(man + ".tmp", man)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replace", action="store_true",
                    help="allow replacing a drive file whose manifest line records other content")
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("rels", nargs="+")
    a = ap.parse_args()
    files = []
    for r in a.rels:
        p = os.path.join(a.src, r)
        if os.path.isdir(p):
            for root, dirs, names in os.walk(p):
                dirs[:] = sorted(d for d in dirs if not d.startswith("."))
                files += [os.path.normpath(os.path.relpath(os.path.join(root, n), a.src)) for n in sorted(names)
                          if not n.startswith(".")]
        else:
            files.append(os.path.normpath(r))
    man = os.path.join(a.dst, "MANIFEST_md5.txt")
    entries = read_manifest(man)
    sums = {rel: (md5(os.path.join(a.src, rel)), os.path.getsize(os.path.join(a.src, rel))) for rel in files}
    clash = [rel for rel in files if rel in entries and entries[rel].split(None, 2)[:2] != [sums[rel][0], str(sums[rel][1])]]
    if clash and not a.replace:
        sys.exit("the manifest records other content for these files; nothing copied (pass --replace to overwrite):\n  "
                 + "\n  ".join(f"{rel}: drive {entries[rel].split(None, 2)[0]}, source {sums[rel][0]}" for rel in clash))
    bad, total = 0, 0
    for rel in files:
        s, d = os.path.join(a.src, rel), os.path.join(a.dst, rel)
        hs, ss = sums[rel]
        os.makedirs(os.path.dirname(d) or ".", exist_ok=True)
        shutil.copyfile(s, d + ".new")
        with open(d + ".new", "rb+") as fh:
            full_fsync(fh)
        t0 = time.time()
        hd, sd = md5(d + ".new", nocache=True), os.path.getsize(d + ".new")
        rate = sd / max(time.time() - t0, 1e-6) / 1e6
        if hs != hd or ss != sd:
            os.remove(d + ".new")
            bad += 1
            print(f"MISMATCH {rel}: src {hs} {ss} B, drive {hd} {sd} B; the earlier drive copy is kept")
            continue
        os.replace(d + ".new", d)
        entries[rel] = f"{hs}  {ss:>14}  {rel}"
        write_manifest(man, entries)
        total += ss
        note = ("  (faster than the drive reads: the re-read may have come from memory; re-hash after a remount "
                "to be sure)") if rate > 500 and sd > 1e6 else ""
        print(f"OK {rel}: {ss:,} B, re-read at {rate:,.0f} MB/s{note}", flush=True)
    print(f"copied {len(files) - bad} of {len(files)} files, {total / 1e9:.2f} GB; md5 + size verified on the drive: "
          f"{len(files) - bad} OK, {bad} mismatches; manifest {man}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
