# Artifacts: Algebraic Cryptanalysis of LOCKER and ROLLO Ideal-LRPC KEMs

Modular proof-of-concept solvers, complexity estimators, and phase-transition falsifiers accompanying [../paper/main.pdf](../paper/main.pdf). All attack implementations, symbolic verifiers, and benchmark harnesses are partitioned by target mechanism and paper section.

---

## Directory Layout

```text
artifacts/
├── run.sh                                      # Unified & per-suite runner (--demo / --full / --suite <name>)
├── fetch_third_party.sh                        # Downloads & SHA-256-verifies NIST Round-1 LOCKER reference on demand
├── run_demo.log                                # Reference stdout captured from ./run.sh --demo on x86_64
├── reference/
│   └── native/                                 # (Fetched on demand) NIST Round-1 LOCKER C/C++ reference (KEMSET-I..IX)
├── ciphertext_attacks/                         # §5: Ciphertext Support-Recovery Solvers
│   ├── schur_locker1/                          # §5.2: kemsolve.cpp, verify_5trials_standalone.py, recover.sage
│   ├── parameterized_tsm/                      # §5.1–§5.3 & §5.5: tsm.cpp, tsm_recover.py, gen_instance.cpp,
│   │                                           #       run_scaling.py, gcd_test.py, design_theory.py, *.sage
│   ├── cyclotomic_ring/                        # §5.5: make_phi_variant.py, kemsolve_phi83.cpp, attack_param.cpp,
│   │                                           #       canonical_ring_component.sage, canonical_ring_frobenius.sage
│   └── sparse_krylov/                          # §5.4: gf2n_pclmul.hpp, sparse_bw.hpp, tsm_bw.cpp, test_gf2n.cpp,
│                                               #       gen_instance.py, cost_model.py, verify_s16.sh
├── pubkey_attacks/                             # §6: Public-Key Secret-Recovery Solvers
│   ├── direct_closed_form/                     # §6.1–§6.2: closed_form_rollo_pk_break.py (2019 ROLLO-I-128)
│   └── subspace_window/                        # §6.3: rolloI_window.py, verify_window_attack.py (2020 ROLLO-I)
└── benchmarks_and_barriers/                    # §7–§8: Phase-Transition Falsifier & 27-Set Complexity Table
    ├── falsifier/                              # §4.2 & §7: s19_falsifier.py, rollo_gf2_engine.py, gf2rank.c,
    │                                           #            exp1_census.py .. exp10_degree2_obstruction.py
    └── estimators/                             # §8: compare_all27_master.py, rsd_estimator.py,
                                                #     rollo_structural_attack.py, frobenius_applicability.py,
                                                #     rollo_params_all27_official.json
```

| Subdirectory | Paper Ref. | Core Verification Target |
|---|---|---|
| `ciphertext_attacks/schur_locker1/` | §5.2, §8.2 | Schur-reduced $`2{,}184 \times 2{,}240`$ solver over $`\mathbb{F}_{2^{83}}`$ recovering the exact 64-byte shared secret on fresh full-strength `LOCKER-I` $`(m=71, n=83, r=5, d=7)`$ instances in $`58\text{ MB}`$ RAM |
| `ciphertext_attacks/parameterized_tsm/` | §5.1–§5.3, §5.5 | Runtime-parameterized dense transposed Support-Minors (`TSM`) solver across all $`r=5`$ `LOCKER` sets (`KEMSET-I, II, IV, VII, VIII`), $`\gcd(m, n) > 1`$ split-ring sweeps, and design-point analysis |
| `ciphertext_attacks/cyclotomic_ring/` | §5.5 | Canonical $`X^n - 1 = (X+1)\Phi_n(X)`$ cyclotomic projection solver (`kemsolve_phi83.cpp`, `make_phi_variant.py`) and symbolic SageMath component checks |
| `ciphertext_attacks/sparse_krylov/` | §5.4, §8.2 | Constant-time `PCLMULQDQ` $`\mathbb{F}_{2^n}`$ arithmetic kernel (`gf2n_pclmul.hpp`), `CSR` sparse Wiedemann/Lanczos engine (`tsm_bw.cpp`, $`1.65\text{ GB}`$ vs. $`565\text{ GB}`$ dense on `LOCKER-III`), and non-asymptotic gate cost model |
| `pubkey_attacks/direct_closed_form/` | §6.1–§6.2, §8.2 | Zero-heuristic single-core public-key recovery via principal-minor readout and closed-form Moore-row gauge recursion on `2019 ROLLO-I-128` $`(79, 47, 5, 6)`$ |
| `pubkey_attacks/subspace_window/` | §6.3 | Transposed subspace-window public-key recovery, scale-invariance verification, and exact Gaussian-binomial weak-key census on `2020 ROLLO-I` |
| `benchmarks_and_barriers/falsifier/` | §4.2, §7.1–§7.2 | Exact $`\mathbb{F}_2`$ nullity law $`\nu_2(P) = n\left(2^r \binom{P-1}{r-1} - m\binom{P-2}{r-1}\right)`$ below the $`m \le 2^r`$ threshold (`W1`–`W12`) and quadratic Plücker monomial deficit verification |
| `benchmarks_and_barriers/estimators/` | §8.1 | Side-by-side 2-input gate complexity ledger across all $`27`$ published `LOCKER` (2017), `ROLLO` (2019), `ROLLO` (2020), and `new2ROLLO` (2020) parameter sets, plus `RYDE` non-applicability controls |
| `reference/native/` *(fetched on demand)* | §3, §8.2 | Official NIST PQC Round-1 `LOCKER` C/C++ submission (`KEMSET-I..IX`) fetched via `./fetch_third_party.sh` for cross-checking key/ciphertext serialization and reference decapsulation |

---

## Environment and Dependencies

### Supported Platforms & CPU Architecture

- **Portable across 64-bit Linux and macOS (`x86_64` and `arm64` / Apple Silicon)**:
  - **Public-key key-recovery suites** (`pubkey_attacks/direct_closed_form/`, `pubkey_attacks/subspace_window/`): pure Python 3.
  - **27-set complexity ledger and phase-transition falsifier** (`benchmarks_and_barriers/estimators/`, `benchmarks_and_barriers/falsifier/`): Python 3 (`numpy` for `falsifier/`; optional `libm4ri` C helper).
  - **Schur-reduced `LOCKER-I` and cyclotomic C++ solvers** ($`\Phi_{83}`$, `ciphertext_attacks/schur_locker1/kemsolve.cpp`, `ciphertext_attacks/cyclotomic_ring/kemsolve_phi83.cpp`): compile and verify end-to-end on both `x86_64` and `arm64`.
    - On `x86_64` with `PCLMULQDQ` (or `AVX-512` `VPCLMULQDQ`), hardware carry-less multiplication solves each `LOCKER-I` instance in $`1.15`$–$`2.15\text{ s}`$ on $`12`$ OpenMP threads ($`\sim 12\text{ s}`$ across 5 trials in `run_demo.log`).
    - On `arm64` (`aarch64` Linux / macOS Apple Silicon), `kemsolve.cpp` automatically activates its built-in scalar carry-less fallback (`KEM_NO_PCLMUL`), verifying exact 64-byte shared-secret recovery in $`\sim 68\text{ s}`$ per trial with OpenMP ($`\sim 130\text{ s}`$ single-threaded). Set `TRIALS=1` for a faster single-trial check on `arm64`.
  - **Sparse Krylov Python utilities & SageMath scripts**: `sparse_krylov/gen_instance.py`, `sparse_krylov/cost_model.py`, and all `.sage` scripts run on any 64-bit architecture.
- **Requires `x86_64` with `PCLMULQDQ` and `SSE2`** (optional `AVX-512` `VPCLMULQDQ` acceleration where available):
  - **`ciphertext_attacks/sparse_krylov/` native C++ kernels** (`gf2n_pclmul.hpp`, `test_gf2n.cpp`, `tsm_bw.cpp`, `verify_s16.sh`): use `_mm_clmulepi64_si128` and 128-bit SSE2 intrinsics (`__m128i`) directly for table-driven $`\mathbb{F}_{2^n}`$ arithmetic ($`1 \le n \le 512`$). On non-`x86_64` hosts, `./run.sh --full` and `./run.sh --suite sparse` automatically report `[SKIP]` for the `x86_64`-only C++ binaries while still executing `gen_instance.py` and `cost_model.py`.
  - **Optional cyclotomic dense C++ research driver** (`ciphertext_attacks/cyclotomic_ring/attack_param.cpp`): invokes `_mm_clmulepi64_si128` directly and requires an `x86_64` `PCLMULQDQ` target (`parameterized_tsm/tsm.cpp` includes a portable scalar fallback for `arm64`).

### Compilers & Native Libraries

`run.sh` and `verify_5trials_standalone.py` respect `CXX` (default `g++`) and `CC` (default `gcc`), automatically probing `-std=c++17`, `-march=native`, `-mpclmul`, and OpenMP (`-fopenmp` on GCC/Linux or `-Xpreprocessor -fopenmp -lomp` with Homebrew `libomp` on macOS Apple Clang, falling back cleanly to single-threaded execution if OpenMP is absent).

| Target / Directory | Compiler & ISA | External C/C++ Libraries |
|---|---|---|
| `schur_locker1/kemsolve.cpp`, `cyclotomic_ring/kemsolve_phi83.cpp` | C++17 (`g++` or `clang++`), optional OpenMP; `x86_64` (`PCLMULQDQ` / `VPCLMULQDQ`) or `arm64` (scalar fallback) | None (C++17 standard library only) |
| `sparse_krylov/{test_gf2n.cpp,tsm_bw.cpp}` | C++17 (`g++` or `clang++`), OpenMP, `x86_64` `PCLMULQDQ`/`SSE2` | None (C++17 standard library only) |
| `cyclotomic_ring/attack_param.cpp` *(optional standalone)* | C++17, OpenMP, `x86_64` `PCLMULQDQ`/`SSE2` | None (C++17 standard library only) |
| `benchmarks_and_barriers/falsifier/gf2rank.c` *(optional accelerator)* | C11 (`gcc` or `clang`), `x86_64` or `arm64` | `libm4ri` (optional; `s19_falsifier.py` falls back to pure Python + `numpy` if absent) |
| `parameterized_tsm/tsm.cpp` *(optional multi-set solver)* | C++17, optional OpenMP; `x86_64` (`PCLMULQDQ`) or `arm64` (scalar fallback) | `libm4ri` (`pkg-config --cflags --libs m4ri`) |
| `reference/native/` & `parameterized_tsm/build_gen.sh` *(optional NIST ref)* | C11 / C++17 (`-fcommon`), `x86_64` or `arm64` | `libntl`, `libgmp`, `libssl`, `libcrypto` (`reference/native/cli.c` also uses GNU `<argp.h>`) |

**Ubuntu / Debian (`apt`) installation**:

```bash
# Core C++17 + OpenMP + Python 3 + NumPy (sufficient for ./run.sh --demo, --full, and all --suite targets)
sudo apt update
sudo apt install -y build-essential g++ pkg-config python3 python3-numpy

# Optional libraries for M4RI rank acceleration, parameterized_tsm/tsm.cpp, and NIST reference builds
sudo apt install -y libm4ri-dev libntl-dev libgmp-dev libssl-dev
```

**macOS (`brew`) installation**:

```bash
# Core OpenMP runtime (enables multi-threading with Apple Clang) + Python 3 + NumPy; or install GNU gcc
brew install libomp pkg-config python numpy
# Optional: GNU GCC toolchain (e.g., CXX=g++-14 ./run.sh --demo)
brew install gcc

# Optional libraries for M4RI rank acceleration and NIST reference instance generators
brew install libm4ri ntl gmp openssl@3 argp-standalone
```

### Python & SageMath Dependencies

- **`Python 3.10+`**:
  - **Standard library only**: `verify_5trials_standalone.py`, `make_phi_variant.py`, `sparse_krylov/{gen_instance.py,cost_model.py}`, `closed_form_rollo_pk_break.py`, `pubkey_attacks/subspace_window/`, and `benchmarks_and_barriers/estimators/`.
  - **`numpy`** (`python3 -m pip install numpy`): required by `benchmarks_and_barriers/falsifier/` (`exp2_rank_truth_table.py`, `exp5_sigma_equivariance.py`–`exp10_degree2_obstruction.py`, and `s19_falsifier.py`).
- **`SageMath` (`sage` 9.5+ / 10.x, optional)**:
  - Not required by `./run.sh`.
  - Used only when running the standalone symbolic algebra and Plücker-recovery scripts: `ciphertext_attacks/schur_locker1/recover.sage`, `ciphertext_attacks/cyclotomic_ring/canonical_ring_*.sage`, and `ciphertext_attacks/parameterized_tsm/{synth.sage,verify_model.sage,verify_reformulation.sage,tsm_recover.py,gcd_test.py,run_scaling.py}`.

---

## Quick Start / Usage

Ensure `run.sh` is executable (`chmod +x run.sh`) and run from inside `artifacts/` (if `reference/native/` is not yet present locally, `./run.sh` and `ciphertext_attacks/parameterized_tsm/build_gen.sh` automatically invoke `./fetch_third_party.sh --ensure` and inform the user on `stderr`; subsequent runs are completely offline):

```bash
# Optional: explicitly pre-fetch / verify the NIST Round-1 LOCKER C/C++ reference
./fetch_third_party.sh

# 1. Fast end-to-end verification:
#    - 5 fresh random LOCKER-I ciphertext recoveries (kemsolve.cpp)
#    - 2019 ROLLO-I-128 full public-key recovery
#    - Partial-Frobenius sigma-equivariance check (exp5)
#    - 27-parameter-set gate-complexity ledger (compare_all27_master.py)
./run.sh --demo

# 2. Comprehensive multi-module verification:
#    Adds cyclotomic Phi_83 solver derivation & build, sparse Krylov check,
#    2020 ROLLO-I subspace-window quick suite, Plucker deficit check (exp10),
#    and ROLLO structural dual-LRPC audit
./run.sh --full

# 3. Individual suite runners for targeted benchmarking or extension:
./run.sh --suite ciphertext   # LOCKER-I 5-trial recovery + cyclotomic Phi_83 build check
./run.sh --suite sparse       # Full sparse Krylov & Wiedemann benchmark (verify_s16.sh on x86_64; Python cost model on arm64)
./run.sh --suite pubkey       # Direct closed-form PK break on 2019 ROLLO-I-128
./run.sh --suite window       # Full 2020 ROLLO-I subspace-window & weak-key verification
./run.sh --suite falsifier    # Full W1-W12 phase-transition falsifier (s19_falsifier.py)
./run.sh --suite ledger       # RYDE/random-code non-applicability + dual-LRPC + 27-set master table
```

### Expected Runtimes by Mode / Suite

| Command | `x86_64` Time | `arm64` Time | Scope |
|---|---|---|---|
| `./run.sh --demo` | $`25`$–$`35\text{ s}`$ | $`\approx 6\text{ min}`$ | End-to-end demo across ciphertext, public-key, falsifier, and ledger suites |
| `./run.sh --full` | $`40`$–$`55\text{ s}`$ | $`\approx 6.5\text{ min}`$ | Checks across all 7 attack and barrier modules |
| `./run.sh --suite ciphertext` | $`10`$–$`15\text{ s}`$ | $`\approx 5.5\text{ min}`$ | 5 fresh `LOCKER-I` encapsulations and `kemsolve_phi83.cpp` build |
| `./run.sh --suite sparse` | $`6`$–$`9\text{ min}`$ | $`< 2\text{ s}`$ | `CSR` Wiedemann solve on `LOCKER-I` and `LOCKER-III` SpMV benchmark |
| `./run.sh --suite pubkey` | $`5`$–$`8\text{ s}`$ | $`5`$–$`8\text{ s}`$ | Single-core Python recovery of `2019 ROLLO-I-128` |
| `./run.sh --suite window` | $`2`$–$`3\text{ min}`$ | $`2`$–$`3\text{ min}`$ | Full `verify_window_attack.py` subspace-window sweep |
| `./run.sh --suite falsifier` | $`45`$–$`65\text{ s}`$ | $`45`$–$`65\text{ s}`$ | All `W1`–`W12` algebraic witnesses |
| `./run.sh --suite ledger` | $`20`$–$`30\text{ s}`$ | $`20`$–$`30\text{ s}`$ | Exhaustive $`\mathrm{Gr}(3, 9)_2`$ controls and 27-set complexity table |

Useful environment overrides:

```bash
# Override C++ compiler, thread count, or LOCKER-I trial count:
CXX=g++-14 OMP_NUM_THREADS=8 TRIALS=1 ./run.sh --demo
```

### Standalone Optional Targets (SageMath & Multi-Set C++ Research Drivers)

```bash
# Symbolic verification of canonical X^n - 1 cyclotomic projection (requires SageMath):
sage ciphertext_attacks/cyclotomic_ring/canonical_ring_component.sage
sage ciphertext_attacks/cyclotomic_ring/canonical_ring_frobenius.sage

# Symbolic verification of the transposed subfield-intersection reformulation theorem
# (run with the Python interpreter of a SageMath environment where sage.all is importable):
"$(dirname "$(realpath "$(command -v sage)")")/python3" \
    ciphertext_attacks/parameterized_tsm/verify_reformulation.sage --m 13 --n 17 --r 3 --trials 3

# Analytic design-point sweep across LOCKER parameter sets I..IX (pure Python):
python3 ciphertext_attacks/parameterized_tsm/design_theory.py

# Multi-set dense C++ TSM solver (requires x86_64 PCLMULQDQ, libm4ri, libntl, libgmp, libssl, and SageMath):
g++ -O3 -march=native -mpclmul -fopenmp -std=c++17 \
    ciphertext_attacks/parameterized_tsm/tsm.cpp \
    $(pkg-config --cflags m4ri) -o ciphertext_attacks/parameterized_tsm/tsm $(pkg-config --libs m4ri)
sh ciphertext_attacks/parameterized_tsm/build_gen.sh I
```

---

## External Reference Sources (Fetched on Demand)

The third-party NIST PQC Round-1 `LOCKER` C/C++ reference submission is not redistributed in this repository; `./fetch_third_party.sh` (invoked automatically when missing) downloads and SHA-256-verifies `LOCKER.zip` and CLI bridge files from upstream (`ethz-spylab/cryptanalysis-benchmark` @ `d243a0a0` / NIST PQC Round-1 archive, MIT) into the gitignored `reference/native/` directory.
