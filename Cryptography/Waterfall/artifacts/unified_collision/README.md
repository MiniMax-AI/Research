# Unified Minimal-Length Full-State Collision Artifact (`unified_collision`)

This directory contains the unified C collision solver (`unified_collide.c`), the single-excursion differential pattern headers (`pattern_7blk.h`, `pattern_8blk.h`, `pattern_9blk.h`), the CNF+XOR SAT lower-bound verifier (`gen_patterns.py`), the two original component-ablation C solvers ($`44`$-block dual-excursion solver `ablation_solvers/route_a/solve5.c` and $`9`$-block quiet-counter coset-scan solver `ablation_solvers/route_b/attack7.c`), pre-computed colliding message pairs, and the unmodified official NIST SHA-3 Round 1 reference implementation (`Waterfall_ref.c`) for Section 3 and Section 5 of the [paper](../../paper/main.pdf).

## Environment and Dependencies

- **Supported Platforms & Architecture**: Linux and macOS (`x86_64` and `arm64`, $`64`$-bit little-endian POSIX environments with $`32`$-bit `unsigned int`, matching the assumptions of `reference/Waterfall_ref.c`).
- **C Compiler**: Standard C99/C11 compiler (`gcc` or `clang` with `-O2`/`-O3` and POSIX `clock_gettime`; configurable via `CC=clang ./run.sh`). Zero external C libraries are required. `unified_collide.c` uses hardware BMI2 `_pext_u32`/`_pdep_u32` instructions when compiled on BMI2-capable `x86_64` hosts and includes a built-in portable C99 fallback for `arm64` and non-BMI2 targets (`run.sh` also supplies an automatic temporary compatibility header when building `ablation_solvers/route_b/attack7.c` on non-BMI2 platforms).
- **Python Environment** (`python3` $`\ge 3.8`$):
  - **`--demo` mode**: 100% Python 3 standard library (no third-party packages needed).
  - **`--full` mode**: Runs out-of-the-box with 100% Python 3 standard library — when `pycryptosat` is absent, `gen_patterns.py` prints a notice and uses the pre-generated `pattern_{7,8,9}blk.h` headers while `run.sh --full` still executes the full $`20`$-seed multi-collision verification sweep. Installing `pycryptosat` (`pip install pycryptosat`) optionally enables `gen_patterns.py` to solve the CNF+XOR SAT instances from scratch and re-verify Theorem 2 ($`B \le 6`$ unsatisfiable; $`B \in \{7, 8, 9\}`$ satisfiable). The offline table generator `ablation_solvers/route_a/gen_data3.py` (whose output `fq_data2.h` is already pre-generated) also uses only the Python 3 standard library.

## Key Files

| File / Directory | Description |
| :--- | :--- |
| `run.sh` | Automated build and verification script (`--demo` and `--full`); compiles into a temporary directory and leaves zero artifacts behind. |
| `unified_collide.c` | Unified C solver implementing Algorithm 1 ($`32`$-element byte-cascade projection, $`\mathsf{Xbox}^{-1}`$ immediate block-counter steering, and step $`t - K_i`$ ring-buffer pre-steering). |
| `pattern_7blk.h` | Minimal steerable $`7`$-block ($`119`$-step, $`448`$-byte) single-excursion differential pattern. |
| `pattern_8blk.h` | Fast $`8`$-block ($`136`$-step, $`512`$-byte) single-excursion differential pattern ($`3.5\text{ ms}`$ solve time). |
| `pattern_9blk.h` | Quiet-counter $`9`$-block ($`153`$-step, $`576`$-byte) single-excursion differential pattern. |
| `gen_patterns.py` | CNF+XOR SAT verifier for Theorem 2 (unsatisfiability of $`B \le 6`$ blocks and satisfiability of $`B \in \{7, 8, 9\}`$ blocks among $`3`$-step look-back steerable trails). |
| `xdiff6.h`, `model.h` | Closed-form $`4`$-round $`\mathsf{Xbox}`$ differential algebra and $`93`$-word Waterfall state model. |
| `ablation_solvers/route_a/` | Original dual-excursion quiet-counter solver (`solve5.c`, $`\Delta s = 544`$ shift + $`32`$-element byte-cascade projection, $`44`$ blocks = $`2{,}816\text{ B}`$, $`45.0\text{ ms}`$). |
| `ablation_solvers/route_b/` | Original quiet-counter coset-scan single-excursion solver (`attack7.c`, modulo-$`32`$ pool parity cancellation + full $`2^{13}`$-coset scan, $`9`$ blocks = $`576\text{ B}`$, $`171.0\text{ ms}`$). |
| `witnesses/` | Pre-computed $`448`$-byte ($`7`$-block) and $`512`$-byte ($`8`$-block) colliding message pairs (`collide_448B_m{1,2}.bin`, `collide_512B_m{1,2}.bin`). |
| `reference/` | Audit driver (`nist_driver.c`) and independent streaming verifier (`verify_reference.c`); `Waterfall_ref.c` and `hash.h` (official NIST SHA-3 Round 1 C reference) are fetched and SHA-256-verified automatically on first run via `../fetch_third_party.sh --ensure`. |

## Usage

Automated end-to-end verification:

```bash
./run.sh --demo
./run.sh --full
```

Manual standalone compilation of the $`7`$-, $`8`$-, or $`9`$-block unified solver (outputs the two colliding messages as hex lines on `stdout` for a given `<seed>` and `<digest_bits>` $`\in \{224, 256, 384, 512\}`$):

```bash
TMP=$(mktemp -d)
cc -O3 -DPATTERN_FILE="\"pattern_8blk.h\"" -o "$TMP/ucol8" unified_collide.c
"$TMP/ucol8" 0 256
rm -rf "$TMP"
```
