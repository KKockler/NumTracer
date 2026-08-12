#!/usr/bin/env python3
"""Helper for noinline_sweep.sh: extract a flow's map-TU compile line and rewrite it for probing.

Kept as a separate file because the compile command is full of shell metacharacters and passing it
back through a command substitution mangles it.

Writes into <tmp>:
  cwd.txt    the directory the compile must run in
  build.sh   the compile, with -o/-c stripped, the arch forced, `-Xptxas -v` added, and the object
             path taken as $1
"""

from __future__ import annotations

import json
import pathlib
import shlex
import sys


def main() -> int:
    app, flow, tmp, arch = sys.argv[1:5]
    db = json.load(open(pathlib.Path(app) / "build" / "compile_commands.json"))
    hit = [e for e in db if e["file"].endswith(f"{flow}/src/CT_map_1.cc")]
    if not hit:
        print(f"no compile entry for {flow}/src/CT_map_1.cc", file=sys.stderr)
        return 1

    toks = shlex.split(hit[0]["command"])
    out, drop_next = [], False
    for t in toks:
        if drop_next:
            drop_next = False
            continue
        if t == "-o":
            drop_next = True
            continue
        if t == "-c" or t.startswith("-arch="):
            continue
        out.append(t)
    out += [f"-arch={arch}", "-Xptxas", "-v", "-c"]

    # compile_commands.json records the *host* compiler (/usr/bin/c++) because Kokkos interposes
    # nvcc_wrapper through CMake's compile launcher (Kokkos_COMPILE_LAUNCHER), which does not appear
    # in the database. Replaying the recorded line verbatim therefore hands nvcc's flags (-Xcudafe,
    # -Xptxas, -arch) to plain g++ and everything fails. Reproduce the real chain.
    launcher = pathlib.Path(app) / "build" / "CMakeCache.txt"
    wrapper = None
    for line in launcher.read_text().splitlines():
        if line.startswith("Kokkos_NVCC_WRAPPER:"):
            wrapper = line.split("=", 1)[1].strip()
    if wrapper and pathlib.Path(wrapper).exists():
        out[0] = wrapper
        out.insert(1, "-ccbin")
        out.insert(2, __import__("os").environ.get("CCBIN", "/usr/bin/g++-14"))

    tmpdir = pathlib.Path(tmp)
    (tmpdir / "cwd.txt").write_text(hit[0]["directory"])
    (tmpdir / "build.sh").write_text("#!/usr/bin/env bash\nset -e\nexec " + " ".join(shlex.quote(x) for x in out) + ' -o "$1"\n')
    return 0


if __name__ == "__main__":
    sys.exit(main())
