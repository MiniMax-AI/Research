# Artifacts: Primal-Shortened ISD Key Recovery on Binary Goppa Codes

This directory contains the self-contained verification suite, live key-recovery solvers, official public keys, and all $`94`$ recovered private-key witnesses accompanying [../paper/main.pdf](../paper/main.pdf).

The artifact suite (`goppa_key_recovery/`) implements and verifies our dual-engine structural key-recovery framework on Track 2 (practical public-key recovery) of the TII McEliece Challenges:

- **Route A (`goppa_key_recovery/route_a/`)**: Primal code shortening to a compact $`n'`$-column window ($`n' \ll n`$), lock-free multi-threaded Information-Set Decoding ($`\mathsf{ISD}`$) harvesting minimum odd-weight codewords ($`w = 2t + 1`$), quadratic interpolation of the characteristic-$`2`$ locator derivative-square identity $`\sigma'_S(X) = \bar{g}(X)^2`$ modulo the public polarization form $`T_{\mathrm{pub}}`$, systematic $`\mathbb{F}_2`$-linear off-diagonal lifting, and $`O(n)`$ closed-form support extension — solving $`89`$ Track-2 challenges (`TII-83` through `TII-254`).
- **Route B (`goppa_key_recovery/route_b/`)**: Shortened-dual degree-$`3`$ ($`p_{\mathsf{HOV}} = 3`$) higher-order vanishing ($`\mathsf{HOV}`$) and single-direction Frobenius-reduced $`\mathsf{HOVER}`$ solver with incremental endomorphism stabilization and truncated Segre–Veronese fiber labeling — solving $`5`$ high-rate Track-2 challenges in the $`(m, t) = (10, 11)`$ family (`TII-172`, `TII-187`, `TII-240`, `TII-246`, `TII-252`).
- **Complete Witness & Verification Bundle (`goppa_key_recovery/witnesses/`, `tii_repo/`, `verify_witnesses_pure.py`)**: All $`94`$ recovered Track-2 secret keys ($`55`$ first-blood solutions and $`39`$ speedups), the unmodified official TII public keys and SageMath checker (`tii_repo/check_solution_track2.sage`), and a zero-dependency pure-Python $`\mathbb{F}_2 / \mathbb{F}_{2^m}[X]`$ verifier checking SHA-256 digests, Theorem 1–2 identities, generator orthogonality, and exact reduced row-echelon parity reconstruction $`\mathrm{RREF}(H_{\mathrm{raw}}(\hat{L}, \hat{g})) = H_{\mathrm{pub}}`$ across all $`94`$ keys.

## Environment and Dependencies

- **Supported Platforms & Architecture**:
  - Linux and macOS on both `x86_64` and `arm64` (Apple Silicon / AArch64) are supported out of the box for witness verification (`--demo` and `--full`) as well as single- and multi-core live key recovery.
  - **Reference Benchmark Hardware**: The single-node wall-clock and compute figures reported in [../paper/main.pdf](../paper/main.pdf) were measured on a `64`-physical-core Dual Intel Xeon Gold 6338 server running Linux `x86_64` (`64` physical cores with OpenMP for Route A; `1` physical core for Route B; $`177.82`$ physical core-hours across $`92`$ single-node instances), with multi-node AMD EPYC CPU clusters (`2,048`–`10,240` cores) used only for the two highest-codimension first-blood challenges (`TII-228` at $`mt = 126`$ and `TII-250` at $`mt = 120`$).
- **Compilers**:
  - A C11 and C++17 compiler (`gcc` / `g++` or `clang` / `clang++`) is required to compile the native Route A helpers (`route_a/minwords_mt.c` and `route_a/rankmask.cpp`). All native binaries are compiled into an ephemeral scratch directory (`mktemp -d`) at runtime, leaving no build artifacts in the repository.
  - **OpenMP (`-fopenmp`)**: Built into GNU `gcc` / `g++` on Linux. On macOS (where `/usr/bin/gcc` and `/usr/bin/g++` invoke Apple Clang without OpenMP enabled by default), install GNU GCC via Homebrew (`brew install gcc`) for multi-threaded $`\mathsf{ISD}`$ sampling. When `$CC` / `$CXX` are unset, the build scripts automatically probe `gcc` / `g++` and Homebrew `gcc-*` / `g++-*` for `-fopenmp`, and fall back cleanly to single-threaded compilation (`#ifdef _OPENMP`) when OpenMP is unavailable.
- **Python & SageMath**:
  - **Required**: Python 3.8+ (standard library only — zero third-party packages are needed to compile and test the native C/C++ helpers and mathematically verify all $`94`$ recovered secret keys via `verify_witnesses_pure.py`).
  - **Optional**: SageMath 10.0+ (`sage` on `$PATH`, bundling `numpy`; reference campaign verified under SageMath `10.9`) to execute live end-to-end Route A and Route B key recovery via `tii_runner.py` and to invoke the unmodified official `tii_repo/check_solution_track2.sage` harness directly.

## Directory Layout

```text
artifacts/
├── README.md                                  # This overview and environment guide
├── run.sh                                     # Top-level entry point (--demo and --full)
├── fetch_third_party.sh                       # Downloads & SHA-256-verifies official TII Track-2 public keys & scripts
└── goppa_key_recovery/
    ├── README.md                              # Detailed module guide for the key-recovery suite
    ├── run.sh                                 # Suite runner invoked by ../run.sh
    ├── verify_witnesses_pure.py               # Zero-dependency Python 3.8+ GF(2)/GF(2^m)[X] verifier
    ├── tii_runner.py                          # Live end-to-end SageMath driver (Route A & Route B)
    ├── run_demo.log                           # Reference log from ./run.sh --demo on 64 Xeon cores
    ├── route_a/
    │   ├── minwords_mt.c                      # Lock-free OpenMP bit-parallel multi-target ISD engine
    │   ├── rankmask.cpp                       # Exhaustive bit-sliced GF(2) symmetric rank-census helper
    │   └── goppa_attack.py                    # Primal-shortened interpolation, T_pub quotient & linear lift
    ├── route_b/
    │   ├── hovlib.py                          # Shortened-dual p=3 HOV operator and kernel construction
    │   ├── hover_direction.py                 # Single-direction Frobenius catalecticant rank selector
    │   ├── selector.py                        # Stacked catalecticant linearization and eigenline ordering
    │   ├── selincr.py                         # Rank-checked incremental E_beta endomorphism solver
    │   ├── bitsolve.py                        # Packed uint64 GF(2) RREF and kernel routines
    │   ├── hover_fibers.py                    # Truncated 450-coordinate Segre-Veronese fiber extraction
    │   ├── labels.py                          # Quadric ideal tangent spaces, fiber frames & conic helper
    │   └── hover_recover.py                   # Squared-multiplier flag, Goppa extraction & O(n) completion
    ├── tii_repo/                              # (Fetched on demand via ../fetch_third_party.sh)
    │   ├── McElieceKeyGen.sage                # Official TII challenge generator & parity-check builder
    │   ├── check_solution_track2.sage         # Unmodified official TII Track-2 key-verification script
    │   ├── key_rec_instances                  # Official parameter catalog for all Track-2 challenges
    │   └── public_keyRec/                     # Official public keys pk_McEliece_<sec>.txt (94 instances)
    └── witnesses/
        ├── sk_McEliece_<sec>.txt              # All 94 recovered Track-2 secret keys in two-line (g, L) format
        ├── verified_tii_<sec>.json            # Per-instance timing, parameter, and RREF certificates
        ├── verified_tii_94_summary.json       # Full 94-challenge campaign & core-hour summary manifest
        ├── results.json                       # SHA-256 digests & SageMath 10.9 check_key=1 audit table
        ├── pass.txt                           # Per-challenge verification summary (94/94 PASS)
        └── run.log                            # Full log of the 94-instance official SageMath verification
```

## Quick Start / Usage

### 1. Automated Verification (`--demo` and `--full`)

From either `artifacts/` or `artifacts/goppa_key_recovery/` (on first run, `run.sh`, `verify_witnesses_pure.py`, and `tii_runner.py` automatically invoke `fetch_third_party.sh --ensure` to fetch and SHA-256-verify `goppa_key_recovery/tii_repo/` with a terminal notice; subsequent runs are completely offline):

```bash
./fetch_third_party.sh   # Optional: explicitly pre-fetch / verify official TII challenge scripts & public keys
./run.sh --demo   # Compiles C/C++ helpers, runs live Sage recovery (if sage is installed), and verifies 13 flagship keys
./run.sh --full   # Runs extended live recovery (Route A + Route B, if sage is installed) and verifies all 94/94 keys
```

- **With `sage` on `$PATH`**: `./run.sh --demo` performs live end-to-end Route A key recovery and official `check_solution_track2.sage` verification on `TII-133, 83, 91, 253` (writing ephemeral outputs to `mktemp -d`), followed by `verify_witnesses_pure.py --demo`. `./run.sh --full` runs live Route A recovery on `TII-133, 83, 91, 125, 253, 254`, live Route B recovery on `TII-172`, and full `94/94` pure-Python RREF and orthogonality verification. You can override the OpenMP thread count via `TII_THREADS=<N> ./run.sh --demo`.
- **Without `sage` (Python 3.8+ & C/C++ compiler only)**: `./run.sh` automatically skips the live SageMath stage and runs `verify_witnesses_pure.py`, which compiles and tests `minwords_mt.c` and `rankmask.cpp` (verifying Theorem 1 and Theorem 2 on `64` live-sampled weight-$`15`$ codewords of `TII-133`), audits all $`94`$ public/secret key SHA-256 digests, and verifies `RREF(H_raw(L, g)) == H_pub` in ~`0.5 s` (`--demo`, $`13`$ keys) or ~`2 s` (`--full`, all $`94`$ keys).

You can also invoke the pure-Python + C/C++ verifier directly from `goppa_key_recovery/` without `run.sh`:

```bash
cd goppa_key_recovery
python3 verify_witnesses_pure.py --demo
python3 verify_witnesses_pure.py --full
```

### 2. Live Single-Instance Key Recovery (Requires SageMath)

To recover a private key from scratch from its public parity-check matrix in `tii_repo/public_keyRec/` without modifying `witnesses/`, point `TII_OUTPUT_DIR` to a scratch directory:

```bash
SAGE_PY="$(sage -c 'import sys; print(sys.executable)')"
OUT_DIR="$(mktemp -d)"

# Route A (Primal-Shortened ISD + Systematic Linear Lift), e.g., TII-133 or TII-253 using 8 OpenMP threads:
TII_OUTPUT_DIR="$OUT_DIR" "$SAGE_PY" tii_runner.py A 133 8
TII_OUTPUT_DIR="$OUT_DIR" "$SAGE_PY" tii_runner.py A 253 8

# Route B (Single-Core Frobenius-Reduced HOVER), e.g., TII-172 (~90 s on 1 core):
TII_OUTPUT_DIR="$OUT_DIR" "$SAGE_PY" tii_runner.py B 172
rm -rf "$OUT_DIR"
```

### 3. Direct Official TII SageMath Checker on Any Shipped Witness

```bash
(cd tii_repo && sage -c 'load("check_solution_track2.sage"); print(check_key(253, "../witnesses/sk_McEliece_253.txt"))')
# Expected output: 1
```

## External Reference Sources (Fetched on Demand)

Official TII McEliece Challenge scripts and Track-2 public keys are not redistributed in this repository; `./fetch_third_party.sh` (invoked automatically when missing) downloads and SHA-256-verifies `McElieceKeyGen.sage`, `check_solution_track2.sage`, `key_rec_instances`, and the 94 Track-2 `public_keyRec/pk_McEliece_<sec>.txt` files from `ElenaKirshanova/tii_decoding_challenge` @ `e502cd59` into the gitignored `goppa_key_recovery/tii_repo/` directory.
