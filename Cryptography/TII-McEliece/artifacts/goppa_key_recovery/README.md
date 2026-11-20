# Artifact: TII McEliece (Binary Goppa Code) Key-Recovery Suite

This suite implements and verifies the dual-engine structural key-recovery attack against the TII McEliece Challenges (`Track 2` = practical public-key recovery in the official challenge repository `ElenaKirshanova/tii_decoding_challenge`, fetched and SHA-256-verified automatically into `tii_repo/` via `../fetch_third_party.sh`) presented in [../../paper/main.pdf](../../paper/main.pdf) (Sections 2–5 and Appendix A): Route A (`route_a/`, Primal Shortening + Multi-Target Information-Set Decoding + Minimum-Weight Locator Derivative Identity $`\sigma'_S(X) = \bar{g}(X)^2`$ modulo $`T_{\mathrm{pub}}`$) and Route B (`route_b/`, Shortened-Dual $`p_{\mathsf{HOV}} = 3`$ Higher-Order Vanishing + Single-Direction Frobenius-Reduced $`\mathsf{HOVER}`$), recovering $`94`$ official TII Track-2 secret keys (`TII-83` through `TII-254`, including $`55`$ first-blood solutions and $`39`$ accelerated instances).

## Key Files

| File / Directory | Description |
| :--- | :--- |
| `run.sh` | Unified entry point (`--demo` and `--full`); compiles native C/C++ helpers into an auto-cleaned `mktemp -d` directory and runs live recovery (when `sage` is available) plus full witness verification. |
| `verify_witnesses_pure.py` | Zero-dependency `Python 3.8+` standard-library $`\mathbb{F}_2 / \mathbb{F}_{2^m}[X]`$ verifier testing the native `minwords_mt` / `rankmask` engines, Theorem 1–2 identities, SHA-256 digests, and exact `check_solution_track2.sage` parity reconstruction ($`\mathrm{RREF}(H_{\mathrm{raw}}(\hat{L}, \hat{g})) = H_{\mathrm{pub}}`$). |
| `tii_runner.py` | End-to-end SageMath driver for live Route A and Route B key recovery and official `check_solution_track2.sage` execution. |
| `route_a/minwords_mt.c` | Lock-free OpenMP bit-parallel Lee–Brickell ($`p = 3`$), Ball-Collision ($`p = 5, 6`$), and L1d-bitmap Stern–Dumer ($`p = 6`$) multi-target $`\mathsf{ISD}`$ sampler. |
| `route_a/rankmask.cpp` | Bit-sliced exhaustive rank-census helper for low-rank subspace isolation over $`\mathbb{F}_2`$. |
| `route_a/goppa_attack.py` | Primal-shortened quadratic relation kernel, $`T_{\mathrm{pub}}`$ polarization quotient, systematic $`\mathbb{F}_2`$-linear off-diagonal lifting, and $`O(n)`$ closed-form support extension. |
| `route_b/` | Shortened-dual $`p_{\mathsf{HOV}} = 3`$ $`\mathsf{HOV}`$ and single-direction Frobenius-reduced $`\mathsf{HOVER}`$ solver (`hovlib.py`, `hover_direction.py`, `hover_fibers.py`, `hover_recover.py`, `labels.py`, `selector.py`, `selincr.py`, `bitsolve.py`). |
| `tii_repo/` *(fetched on demand)* | Official TII challenge generator and verifier (`McElieceKeyGen.sage`, `check_solution_track2.sage`, `key_rec_instances`, and $`94`$ public keys in `public_keyRec/pk_McEliece_<sec>.txt`), downloaded and SHA-256-verified automatically on first run via `../fetch_third_party.sh --ensure`. |
| `witnesses/` | All $`94`$ recovered secret keys (`sk_McEliece_<sec>.txt`), SHA-256 and SageMath `10.9` verification certificates (`results.json`, `pass.txt`, `run.log`), and per-challenge timing manifests (`verified_tii_<sec>.json`, `verified_tii_94_summary.json`). |
| `run_demo.log` | Reference log captured on the `64`-core Xeon benchmark server. |

## Environment and Dependencies

- **Platforms**: Linux and macOS (`x86_64` and `arm64`). Reference benchmarks in [../../paper/main.pdf](../../paper/main.pdf) were recorded on a `64`-physical-core Dual Intel Xeon Gold 6338 server running Linux `x86_64` ($`92`$ single-node instances in $`177.82`$ physical core-hours), plus multi-node AMD EPYC CPU clusters for `TII-228` ($`mt = 126`$) and `TII-250` ($`mt = 120`$).
- **Required**: `bash`, a C11/C++17 compiler (`gcc`/`g++` or `clang`/`clang++`), and `Python 3.8+` (standard library only: `hashlib`, `json`, `os`, `pathlib`, `re`, `shutil`, `struct`, `subprocess`, `sys`, `tempfile`).
  - **OpenMP (`-fopenmp`)**: Built into GNU `gcc`/`g++` on Linux; on macOS, install GNU GCC via `brew install gcc` for multi-core $`\mathsf{ISD}`$ sampling. When `$CC`/`$CXX` are unset, `run.sh`, `verify_witnesses_pure.py`, and `tii_runner.py` automatically detect `gcc`/`g++` or Homebrew `gcc-*`/`g++-*` with `-fopenmp`, and fall back cleanly to single-threaded C11/C++17 compilation (`#ifdef _OPENMP`) if OpenMP is unavailable.
- **Optional**: SageMath 10.0+ (`sage` on `$PATH`, bundling `numpy`; reference runs used SageMath `10.9`). When `sage` is absent, `run.sh` runs `verify_witnesses_pure.py` to compile and exercise the native C/C++ helpers and verify all shipped witnesses against `tii_repo/public_keyRec/` using pure-Python $`\mathbb{F}_2 / \mathbb{F}_{2^m}[X]`$ arithmetic.

## Usage

```bash
./run.sh --demo   # ~0.5 s pure-Python + C/C++ fallback; ~8 s with live SageMath Route A recovery
./run.sh --full   # ~2.0 s pure-Python + C/C++ fallback (94/94 keys); ~2–5 min with live SageMath Route A + B
```

- **`--demo` (default)**: compiles `route_a/minwords_mt.c` and `route_a/rankmask.cpp` into a temporary directory; if `sage` is on `$PATH`, performs live Route A key recovery and official `check_solution_track2.sage` verification on `TII-133, 83, 91, 253` (override OpenMP thread count via `TII_THREADS=<N>`); then runs `verify_witnesses_pure.py --demo` to verify the Theorem 1 identity $`\sigma'_S(X) = \bar{g}(X)^2`$ and Theorem 2 pair-count parity $`\binom{2t+1}{2} \equiv t \pmod{2}`$ on $`64`$ live-sampled weight-$`15`$ codewords of `TII-133`, audit all $`94`$ SHA-256 digests and structural key parameters, and execute the exact Track-2 acceptance test (`check_key == 1`) on $`13`$ representative Route A and Route B challenges (`83, 133, 91, 125, 253, 254, 228, 250, 172, 187, 240, 246, 252`).
- **`--full`**: runs live SageMath Route A (`133, 83, 91, 125, 253, 254`) and Route B (`172`) recovery when `sage` is available, and executes the full RREF and generator orthogonality acceptance test ($`\mathrm{RREF}(H_{\mathrm{raw}}(\hat{L}, \hat{g})) = H_{\mathrm{pub}}`$) across all $`94 / 94`$ solved TII Track-2 secret keys.

### Standalone Pure-Python Verification & Live Single-Instance Recovery

```bash
# Zero-dependency pure-Python + C/C++ verification (no SageMath needed):
python3 verify_witnesses_pure.py --demo
python3 verify_witnesses_pure.py --full

# Live single-instance key recovery via SageMath into a scratch directory:
SAGE_PY="$(sage -c 'import sys; print(sys.executable)')"
OUT_DIR="$(mktemp -d)"
TII_OUTPUT_DIR="$OUT_DIR" "$SAGE_PY" tii_runner.py A 253 8
TII_OUTPUT_DIR="$OUT_DIR" "$SAGE_PY" tii_runner.py B 172
rm -rf "$OUT_DIR"
```
