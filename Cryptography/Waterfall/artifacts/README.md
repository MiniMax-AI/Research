# Verification Artifacts — Practical Collision Cryptanalysis of Waterfall

This directory contains the self-contained C collision solvers, differential trail headers, SAT lower-bound verifier, pre-computed colliding message pairs, and unmodified official NIST SHA-3 Round 1 reference implementation accompanying the paper [Practical Collision Cryptanalysis of Waterfall](../paper/main.pdf).

All artifacts verify practical full-state collisions ($`93`$ internal $`32`$-bit state words = $`2{,}976`$ bits) across all four NIST security levels (`Waterfall-224`, `Waterfall-256`, `Waterfall-384`, and `Waterfall-512`) against both the unmodified official NIST C reference (`Waterfall_ref.c`) and an independent streaming evaluator (`verify_reference.c`).

## Environment and Dependencies

- **Supported Platforms & Architecture**: Linux and macOS (`x86_64` and `arm64`, $`64`$-bit little-endian POSIX environments with $`32`$-bit `unsigned int`, as assumed by the official NIST reference implementation `Waterfall_ref.c`).
- **C Compiler**: Any standard C99/C11 compiler (`gcc` or `clang` with `-O2`/`-O3` and POSIX `clock_gettime`; override via `CC=clang ./run.sh`). Zero external C libraries are required. On `x86_64` CPUs supporting BMI2, hardware `_pext_u32`/`_pdep_u32` instructions are used automatically; on `arm64` (e.g., Apple Silicon) and non-BMI2 `x86_64` targets, portable C99 inline equivalents are used automatically with bit-identical output.
- **Python Environment** (`python3` $`\ge 3.8`$):
  - **`--demo` mode**: 100% Python 3 standard library (zero third-party packages required).
  - **`--full` mode**: Runs end-to-end with 100% Python 3 standard library — if `pycryptosat` is not installed, `gen_patterns.py` automatically uses the pre-generated `pattern_{7,8,9}blk.h` headers while `run.sh --full` still executes the live $`20`$-seed multi-collision sweep. To optionally re-solve the CNF+XOR SAT instances of Theorem 2 ($`B \le 6`$ unsatisfiable and $`B \in \{7, 8, 9\}`$ satisfiable) live from scratch, install `pycryptosat`:
    ```bash
    pip install pycryptosat
    ```

## Directory Layout

| Path | Description |
| :--- | :--- |
| `run.sh` | Top-level convenience wrapper forwarding all arguments (`--demo` or `--full`) to `unified_collision/run.sh`. |
| `fetch_third_party.sh` | Downloads and SHA-256-verifies the official NIST SHA-3 Round 1 Waterfall C reference (`Waterfall_ref.c`, `hash.h`) into `unified_collision/reference/` when missing. |
| `unified_collision/` | Unified $`7`$-, $`8`$-, and $`9`$-block C collision solvers (`unified_collide.c`), differential pattern headers (`pattern_{7,8,9}blk.h`), algebraic $`\mathsf{Xbox}`$ differential model (`xdiff6.h`, `model.h`), Theorem 2 CNF+XOR SAT verifier (`gen_patterns.py`), original $`44`$-block dual-excursion (`solve5.c`) and $`9`$-block quiet-counter coset-scan (`attack7.c`) ablation solvers (`ablation_solvers/`), pre-computed $`448`$-byte and $`512`$-byte collision witnesses (`witnesses/`), and NIST/independent C verification drivers (`reference/`). |

## Quick Start / Usage

Run directly from `artifacts/` (or after `cd unified_collision`). On the first run, `run.sh` automatically invokes `./fetch_third_party.sh --ensure` to fetch and SHA-256-verify `Waterfall_ref.c` and `hash.h` (printing a terminal notice); subsequent runs are completely offline:

```bash
./fetch_third_party.sh   # Optional: explicitly pre-fetch / verify NIST Waterfall reference files
./run.sh --demo   # Verifies 448 B & 512 B witnesses, dual-excursion (solve5.c, 2,816 B) & quiet-counter coset-scan (attack7.c, 576 B) solvers, and live 9/8/7-block collisions across Waterfall-224/256/384/512 (~1–4 s)
./run.sh --full   # Additionally runs Theorem 2 SAT lower-bound check (gen_patterns.py) and a 20-seed live 512-byte multi-collision sweep (~2–5 s)
```

Equivalently, from inside `unified_collision/`:

```bash
cd unified_collision
./run.sh --demo
./run.sh --full
```

All temporary binaries and intermediate outputs are built inside an isolated `mktemp -d` directory and cleaned up automatically upon exit.

## External Reference Sources (Fetched on Demand)

The third-party NIST SHA-3 Round 1 Waterfall reference implementation (`Waterfall_ref.c`, `hash.h` by Bob Hattersley, 2008) is not redistributed in this repository; `./fetch_third_party.sh` (called automatically by `./run.sh` when missing) downloads and SHA-256-verifies `Waterfall.zip` from upstream (`ethz-spylab/cryptanalysis-benchmark` @ `d243a0a0` / NIST SHA-3 archive mirror) into gitignored paths under `unified_collision/reference/`.
