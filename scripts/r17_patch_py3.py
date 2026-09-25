# -*- coding: utf-8 -*-
"""Minimal Python 2 -> 3 compatibility patcher for bulik/ldsc v1.0.1.

Applied transforms (mechanical only; no algorithmic change):
  1. xrange(  -> range(
  2. .ix[      -> .loc[
  3. print >>fh, expr  ->  fh.write(str(expr) + "\n")
  4. bare print statements -> print(...) function call
  5. iteritems() -> items()
Every patched line is logged to patch_log.txt in the same directory.
"""
import io
import os
import re
import sys

FILES = [
    "ldsc.py",
    "munge_sumstats.py",
    os.path.join("ldscore", "__init__.py"),
    os.path.join("ldscore", "parse.py"),
    os.path.join("ldscore", "sumstats.py"),
    os.path.join("ldscore", "regressions.py"),
    os.path.join("ldscore", "jackknife.py"),
    os.path.join("ldscore", "irwls.py"),
    os.path.join("ldscore", "ldscore.py"),
    os.path.join("ldscore", "score.py"),
]

log_lines = []


def patch_source(src):
    out = []
    for lineno, line in enumerate(src.splitlines(keepends=True), 1):
        orig = line
        stripped = line.lstrip()
        indent = line[: len(line) - len(stripped)]

        # 3. print >>fh, expr
        m = re.match(r"^(\s*)print\s*>>(\w+),\s*(.*)$", line.rstrip("\n"))
        if m:
            ind, fh, expr = m.groups()
            line = '%s%s.write(str(%s) + "\\n")\n' % (ind, fh, expr)
        elif re.match(r"^\s*print\s+[^=(]", line) and not stripped.startswith("print("):
            # 4. bare print statement -> print function
            expr = stripped[len("print"):].strip()
            line = "%sprint(%s)\n" % (indent, expr)
        else:
            # 1/2/5 plain textual replacements
            line = line.replace("xrange(", "range(")
            line = line.replace(".ix[", ".loc[")
            line = line.replace(".iteritems()", ".items()")

        if line != orig:
            log_lines.append("L%d: %r -> %r" % (lineno, orig.rstrip("\n"), line.rstrip("\n")))
        out.append(line)
    return "".join(out)


def main(root):
    for rel in FILES:
        path = os.path.join(root, rel)
        if not os.path.exists(path):
            continue
        with io.open(path, "r", encoding="utf-8") as fh:
            src = fh.read()
        patched = patch_source(src)
        with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(patched)
    with io.open(os.path.join(root, "patch_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log_lines))
    print("patched %d lines across %d files" % (len(log_lines), len(FILES)))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
