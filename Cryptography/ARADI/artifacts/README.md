# $`\mathsf{OPTIC}`$ Artifact Suite: Algebraic Verification, SOTA Audit, and 5R–12R Cryptanalysis of $`\mathsf{ARADI}`$

This directory contains the self-contained, single-core Python 3 and C99 verification and key-recovery artifacts accompanying [../paper/main.pdf](../paper/main.pdf) (*Integral Cryptanalysis of ARADI*). All implementations target the unmodified $`128`$-bit block, $`256`$-bit key $`\mathsf{ARADI}`$ specification published by the U.S. National Security Agency (IACR ePrint 2024/1240) starting at Round $`r_0 = 0`$:

- **Specification & Key-Schedule Inversion (`§2.1`–`§2.4`):** Bit-exact reference implementation verified against all $`21`$ official test vectors, $`17`$ round keys, and $`16`$ intermediate states (IACR ePrint 2024/1240, Figures 3.4, 4.3, and 4.4), together with closed-form $`O(1)`$ master-key inversion $`\mathsf{KS}^{-1}(K^{(2m)}, K^{(2m+1)}, 2m)`$ for every even index $`2m \in \{0, 2, \dots, 14\}`$.
- **Structural Algebra & Prior-Work Audits (`§3.1`–`§3.4`):** Formal computer-assisted verification of linear-layer involutivity ($`\Lambda_j^2 = I_{128}`$), unique orthogonality of $`L_3 = (8, 9, 7)`$ ($`L_3^\top L_3 = I_{32}`$) and its eight $`4`$-column islands, exact $`16`$-bit super-S-box ANF degree spectrum, composite linear-layer canceling-mask census ($`1{,}364`$ products for $`k \le 5`$), $`5`$-round blind-column audit of prior truncated-differential claims (*LATINCRYPT 2025*, IACR ePrint 2025/1443 `§6.2`), analytical zero-sampling derivation of all $`198`$ universal $`5`$-round zero-sum equations on a $`48`$-CP prefix, and live $`6`$-round $`2^{16}`$-CP colliding-cube distinguisher ($`42`$-dimensional null space).
- **Practical Single-Core Key Recovery Across 5–8 Rounds (`§4.1`–`§4.2`):** Unconditional all-key ($`100\%`$ of $`2^{256}`$ keys) subkey and full $`256`$-bit master-key recovery on $`5`$ rounds ($`6`$ CP in $`\approx 11\text{ ms}`$ for $`K^{(5)}`$; $`4`$ CP in $`\approx 16.0\text{ s}`$ for the full master key), $`6`$ rounds ($`48`$ CP in $`\approx 2.4\text{ s}`$; $`32\times`$ less data than the $`1{,}536`$-CP prototype in IACR ePrint 2025/1227), $`7`$ rounds ($`3{,}072`$ CP in $`\approx 3.5\text{ s}`$; $`> 2^{20}\times`$ less data than *LATINCRYPT 2025*), and $`8`$ rounds ($`90{,}112`$ CP in $`\approx 6.0\text{ s}`$ mean, $`< 20\text{ MB}`$ RAM).
- **Higher-Round Kernel & Complexity Verification Across 9–12 Rounds (`§4.3`–`§4.4`, `§5.1`):** Single-core C99 verification of the inner $`9`$-round streaming joint-histogram Walsh–Hadamard ($`\mathsf{JM\text{-}WHT}`$) kernel under a fixed outer $`\kappa`$ ($`D = 2^{22.32}`$ CP, full attack $`T = 2^{54.77}`$ butterflies, $`4.7\text{ MB}`$ RAM, $`\rho = 256`$), and exact $`\mathbb{F}_2`$ key-schedule rank and progressive list-join complexity auditing across $`9`$ to $`12`$ rounds under the packed-integer cost model of IACR ePrint 2025/1227 (and compared against weak-key differential bounds in *ToSC* 2026(3)).

For full mathematical proofs, algorithmic descriptions, and comparison tables, see [../paper/main.pdf](../paper/main.pdf).

---

## Environment and Dependencies

### Supported Platforms & Architecture

- **Operating Systems:** Linux and macOS (`Darwin`).
- **CPU Architectures:** `x86_64` (`amd64`) and `arm64` (`aarch64` / Apple Silicon). All Python and C99 implementations use portable, endian-agnostic fixed-width integer arithmetic with no architecture-specific SIMD intrinsics.
- **Single-Core Execution:** Every experiment is designed for deterministic single-core execution (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `VECLIB_MAXIMUM_THREADS=1`, all exported automatically by `run.sh`). Peak memory footprint across all live executables is under $`20\text{ MB}`$ RAM.
- **Benchmark Environment Note:** Reference single-core wall-clock timings reported in [../paper/main.pdf](../paper/main.pdf) (`§5.1`, Table 3) and archived in `logs/` were recorded on a single `x86_64` Linux CPU core with `pycryptosat` installed (`~30–70 s` for `--demo`, `~2–8 min` for `--full`). On modern `arm64` CPUs (such as Apple Silicon), the vectorized FHT and C99 stages typically run $`1.5\times`$–$`4\times`$ faster (`~5–6 s` for `--higher-only`, `~20–25 s` for `--demo` without `pycryptosat`).

### C Compiler

- **Compiler:** Either `gcc` (GNU GCC) or `clang` (LLVM / Apple Clang) supporting standard ISO C99 (`-std=c99`).
- **Flags & Libraries:** Compiled with `-std=c99 -O3 -Wall -Wextra` and linked against `-lm` (`libm`, used by `attack_9r_nested_fht.c` for `log2`).
- **Automatic Hardware Tuning & Cleanup:** `run.sh` honors the `CC` environment variable (defaulting to `gcc` with automatic fallback to `cc`), probes compiler support for `-march=native` / `-mcpu=native` automatically, and places all compiled binaries inside a temporary `mktemp -d` scratch directory that is removed automatically on exit via `trap`.

### Python Environment

- **Interpreter:** `Python 3.10+` (uses `int.bit_count()`; override the interpreter in `run.sh` via `PYTHON=python3`). Note that `aradi_ref.py` depends exclusively on the Python 3 standard library.
- **Required Package (`numpy`):** Used by `attack_7r_fht.py` ($`12`$-bit FHT partial-sum kernel) and `verify_higher_rounds.py` ($`16`$-variable binary Möbius transform and vectorized $`2^{16}`$-point cube verification), and for accelerated $`48`$-CP rank sampling in `verify_algebra_and_audit.py`. Install via:
  ```bash
  pip install numpy
  ```
- **Optional Package (`pycryptosat`):** Python bindings for the CryptoMiniSat CNF-XOR solver, used only for live SAT solving in two depth-$`1`$ sub-stages: `attack_5r_lowdata.py` (Part B: $`5`$-round $`4`$-CP $`4`$-layer subspace peel) and `attack_6r_48cp.py` (Layer 1: $`6`$-round $`48`$-CP one-hot XOR recovery of $`K^{(6)}`$). Install via:
  ```bash
  pip install pycryptosat
  ```
  **Graceful Pure-Python Fallback (No `pycryptosat`):** All remaining stages—`aradi_ref.py`, `verify_algebra_and_audit.py`, `attack_5r_lowdata.py` (Part A: $`5`$-round $`6`$-CP $`2^{12}`$ $`8`$-island MITM), `attack_7r_fht.py`, `attack_8r_mitm.c`, `verify_higher_rounds.py`, and `attack_9r_nested_fht.c`—run completely unassisted in pure Python / C99 without `pycryptosat`. When `pycryptosat` is not installed, `attack_5r_lowdata.py` and `attack_6r_48cp.py` print an explicit `[Note]` message and execute a deterministic pure-Python fallback:
  - **5R Part B (`attack_5r_lowdata.py`, `4` CP):** Verifies all $`4`$-layer differential/cube subspace annihilators ($`\dim V_k \in \{3, 12, 36, 88\}`$) and all $`256`$ Law-III affine key-schedule consistency equations directly on the witness round keys $`(K^{(2)}, K^{(3)}, K^{(4)}, K^{(5)})`$, then executes live closed-form $`\mathsf{KS}^{-1}(K^{(2)}, K^{(3)}, 2)`$ master-key reconstruction and $`4`$-block re-encryption verification.
  - **6R (`attack_6r_48cp.py`, `48` CP):** Verifies all $`198`$ analytically derived Layer-1 one-hot XOR equations directly on $`K^{(6)}`$, peels Round $`5`$ on the $`48`$ ciphertexts, and runs the live, unassisted Layer-2 $`2^{12}`$ $`8`$-island MITM recovery of $`K^{(5)}`$ (`~28–72 ms`), Layer-3 $`16`$-point subcube + $`2`$-round dead-column recovery of $`K^{(4)}`$ (`~21–108 ms`), and $`O(1)`$ $`\mathsf{KS}^{-1}(K^{(4)}, K^{(5)}, 4)`$ master-key recovery verified across all $`48`$ plaintexts.

---

## Quick Start / Usage

Ensure `run.sh` is executable (`chmod +x run.sh`) and run from inside `artifacts/`:

```bash
# 1. Fast higher-round check: 16-bit super-S-box spectrum, 6R d=16 distinguisher,
#    9R–12R progressive list-join rank/complexity audit, and 2x 9R inner JM-WHT checks (~5–8 s, 1 core)
./run.sh --higher-only

# 2. Fast end-to-end verification across all 7 stages (5R–8R live + 9R–12R kernel/audit; default mode)
#    (~20–25 s without pycryptosat; ~35–70 s with pycryptosat, 1 core)
./run.sh --demo

# 3. Full multi-seed statistical suite matching §5.1 and Table 3 of ../paper/main.pdf
#    (200x 5R 6-CP, 10x 5R 4-CP, 10x 6R, 10x 7R, 10x 8R, plus 9R–12R kernel/complexity checks; ~2–8 min, 1 core)
./run.sh --full
```

Optional environment overrides supported by `run.sh`:

```bash
CC=clang SEED=20260409 PYTHON=python3 ./run.sh --demo
```

### Standalone Script & Binary Invocations

Individual Python modules and C99 solvers can also be run directly (compiling C binaries into a temporary directory keeps `artifacts/` clean):

```bash
# Official NSA ARADI test vectors & O(1) even-pair key-schedule inversion (stdlib only)
python3 -B aradi_ref.py

# Linear-layer algebra, LATINCRYPT 2025 blind-column audit & 198-eq certificate
python3 -B verify_algebra_and_audit.py

# 5-round key recovery: 6-CP 8-island MITM (K^(5)) and 4-CP 4-layer peel (256-bit MK)
python3 -B attack_5r_lowdata.py --trials-6cp 25 --trials-4cp 2 --seed 20260409

# 6-round 48-CP full 256-bit master-key recovery (zero-sampling certificate + 3-layer cascade)
python3 -B attack_6r_48cp.py --trials 2 --seed 20260409

# 7-round full 256-bit master-key recovery (--shifts 6 => 3,072 CP; --shifts 8 => 4,096 CP)
python3 -B attack_7r_fht.py --trials 3 --shifts 6 --seed 20260707

# Higher-round algebraic theorems, live 6R d=16 cube check & 9R–12R complexity audit
python3 -B verify_higher_rounds.py

# 8-round 90,112-CP root-free 3-branch FHT + 2 MB Bloom-MITM full master-key recovery (C99)
TMP=$(mktemp -d)
cc -std=c99 -O3 -Wall -Wextra attack_8r_mitm.c -o "$TMP/attack_8r_mitm" -lm
"$TMP/attack_8r_mitm" --trials 2 --seed 20260409

# 9-round inner streaming JM-WHT kernel & cross-shift MK closure check (fixed outer kappa)
cc -std=c99 -O3 -Wall -Wextra attack_9r_nested_fht.c -o "$TMP/attack_9r_nested_fht" -lm
"$TMP/attack_9r_nested_fht" --trials 2 --kbits 0 --cubes 5 --shifts 16
rm -rf "$TMP"
```

---

## Directory Layout

```text
artifacts/
├── README.md                      # Documentation, dependencies, and reproduction guide
├── run.sh                         # Single-core test runner (--higher-only, --demo, --full)
├── aradi_ref.py                   # Pure-stdlib ARADI reference, KS^-1, and 21 official test vectors
├── verify_algebra_and_audit.py    # L_3 orthogonality proof, blind-column audit & 198-eq zero-sum certificate
├── verify_higher_rounds.py        # Super-S-box ANF spectrum, 1-Round Lift, 6R d=16 cube & 9R–12R audit
├── attack_5r_lowdata.py           # 5-round 6-CP island MITM (K^(5)) & 4-CP 4-layer peel (256-bit MK)
├── attack_6r_48cp.py              # 6-round 48-CP zero-sampling full 256-bit master-key recovery
├── attack_7r_fht.py               # 7-round 3,072-CP colliding-cube FHT full 256-bit master-key recovery
├── attack_8r_mitm.c               # 8-round 90,112-CP 3-branch FHT + Bloom-MITM master-key solver (C99)
├── attack_9r_nested_fht.c         # 9-round streaming JM-WHT peel & cross-shift MK verifier (C99)
└── logs/                          # Archived single-core reference logs matching ../paper/main.pdf
    ├── run_demo.log
    ├── run_full.log
    ├── verify_higher_rounds.log
    ├── attack_9r_20trials.log
    ├── attack_9r_splitmix64_demo.log
    └── distinguishers_6r_nullspace_audit.log
```

### Detailed Component Reference

| File | Paper Ref. | Target & Verified Claims |
| :--- | :--- | :--- |
| `run.sh` | `§5.1` | Single-core orchestrator (`--higher-only`, `--demo`, `--full`); enforces single-thread BLAS/OpenMP variables, compiles C99 sources into an auto-cleaned scratch directory, and reports overall `RESULT: SUCCESS` / `FAIL`. |
| `aradi_ref.py` | `§2.1`–`§2.4` | Dependency-free Python 3 reference implementation (`sbox`, `sbox_inv`, `linear_L`, `expand_key`, `encrypt_rounds`, `decrypt_rounds`, and $`O(1)`$ even-pair inverse key schedule `master_key_from_even_pair`) verified against all $`21`$ official test vectors, $`17`$ round keys, and $`16`$ intermediate states from IACR ePrint 2024/1240. |
| `verify_algebra_and_audit.py` | `§2.2`–`§2.3`, `§3.1`–`§3.2` | Verifies $`\Lambda_j^2 = I_{128}`$ ($`j = 0\dots 3`$), proves $`L_3 = (8, 9, 7)`$ is the unique symmetric orthogonal word map ($`L_3^\top L_3 = I_{32}`$, $`255`$ non-zero canceling masks) with $`8`$ disjoint $`4`$-column islands, proves single-column $`3`$-round truncated differentials leave $`8`$ or $`13`$ blind columns ($`32`$ or $`52`$ unconstrained key bits; *LATINCRYPT 2025*, IACR ePrint 2025/1443 `§6.2`), and derives all $`198`$ universal $`5`$-round zero-sum equations analytically in $`\approx 5\text{ ms}`$ ($`0`$ key samples). |
| `attack_5r_lowdata.py` | `§3.1`, `§4.1` | Implements both $`5`$-round single-key chosen-plaintext attacks: (Part A) $`6`$-CP ($`6.21`$ CP mean with $`1`$-pair adaptive tie-break) $`2^{12}`$ $`8`$-island MITM recovery of $`128`$-bit $`K^{(5)}`$; (Part B) $`4`$-CP $`2`$D affine-plane $`4`$-layer subspace peel ($`\dim V_k \in \{3, 12, 36, 88\}`$) coupled with $`256`$ Law-III key-schedule equations and $`\mathsf{KS}^{-1}(K^{(2)}, K^{(3)}, 2)`$ for full $`256`$-bit master-key recovery. |
| `attack_6r_48cp.py` | `§3.2`, `§4.1` | Recovers the full $`256`$-bit master key of $`6`$-round $`\mathsf{ARADI}`$ from $`48`$ chosen plaintexts ($`2^{5.58}`$ CP; $`32\times`$ less data than the $`1{,}536`$-CP FHT prototype in IACR ePrint 2025/1227) via the $`198`$-equation analytical certificate ($`K^{(6)}`$) followed by a zero-extra-query cascade through $`K^{(5)}`$ ($`2^{12}`$ island MITM) and $`K^{(4)}`$ ($`2`$-round dead-column peel) and $`\mathsf{KS}^{-1}(K^{(4)}, K^{(5)}, 4)`$. |
| `attack_7r_fht.py` | `§3.3`, `§4.2` | Recovers the full $`256`$-bit master key of $`7`$-round $`\mathsf{ARADI}`$ from $`3{,}072`$ chosen plaintexts ($`6 \times 2^9 = 2^{11.58}`$ CP) via $`L_0`$-colliding $`9`$D cubes, progressive $`2^{12}`$-point FHT column overlap for $`K^{(7)}`$, Round-$`6`$ peel onto $`19`$ balanced columns of $`Z_4 = \Lambda_0(Q_5)`$ for $`K^{(6)}`$, and $`\mathsf{KS}^{-1}(K^{(6)}, K^{(7)}, 6)`$. |
| `attack_8r_mitm.c` | `§3.3`, `§4.2` | C99 implementation recovering the full $`256`$-bit master key of $`8`$-round $`\mathsf{ARADI}`$ from $`90{,}112`$ chosen plaintexts ($`22 \times 2^{12} = 2^{16.46}`$ CP, $`1.44\text{ MB}`$ ciphertexts, $`< 20\text{ MB}`$ total RAM) via root-free $`3`$-branch additive-split $`2^{16}`$-entry FHT tables, a $`2\text{ MB}`$ L3-resident $`24`$-bit Bloom filter ($`2^{32}`$ probes at $`b = 25`$), shift-$`8`$ cascade, and $`\mathsf{KS}^{-1}(K^{(6)}, K^{(7)}, 6)`$. |
| `verify_higher_rounds.py` | `§2.3`, `§3.2`–`§3.4`, `§4.4` | Deterministic single-core auditor verifying: (1) the exact $`2^{16}`$-mask ANF degree histogram `{0: 1, 3: 3, 5: 28, 6: 224, 8: 65280}` of the $`16`$-bit $`L_3`$-island super-S-box $`\Psi_G`$; (2) composite linear-layer canceling-mask census across all $`1{,}364`$ products ($`k = 1\dots 5`$, reducing exclusively to $`I_{32}`$ or $`L_3`$: $`98`$ hits); (3) $`\mathbb{Z}_{16}`$ intra-half shift symmetry, the $`1`$-Round Lift identity, and live $`6`$-round $`2^{16}`$-CP cube `(0, 5, 19, 24)` checks ($`Q_5`$ nullity $`128`$; $`Q_6`$ column $`25`$ + $`24`$ lifted masks on $`12`$ ready columns); (4) exact $`\mathbb{F}_2`$ rank $`\rho = 256`$ for $`9`$-round shifted cubes; and (5) progressive list-join sizes $`\log_2 \vert L_t \vert`$ and packed-integer complexities for $`10`$R, $`11`$R, and $`12`$R vs. IACR ePrint 2025/1227. |
| `attack_9r_nested_fht.c` | `§4.3`, `§5.1` | C99 inner-kernel verifier for $`9`$-round $`\mathsf{ARADI}`$ ($`D = 80 \cdot 2^{16} = 2^{22.32}`$ CP, full attack $`T = 2^{54.77}`$ butterflies, $`4.7\text{ MB}`$ RAM, $`\rho = 256`$): encrypts $`16`$ $`\mathbb{Z}_{16}`$-shifted translates of the $`d = 16`$ colliding cube `(0, 5, 19, 24)`, streams ciphertexts into a $`4\text{ KB}`$ $`12`$-bit joint parity histogram $`\mathsf{JM}`$, and verifies a single iteration of the $`2^{16}`$-candidate $`(\gamma_E, g_u)`$ $`\mathsf{WHT}_{12}`$ / $`\mathsf{WHT}_4`$ kernel, cross-shift $`\gamma`$ propagation, and $`\mathsf{KS}^{-1}(K^{(8)}, K^{(9)}, 8)`$ closure under a fixed outer $`\kappa`$. |
| `logs/` | `§5.1` | Pre-generated single-core reference logs: `run_demo.log` (`5R`–`8R` demo), `run_full.log` (`200` $`\times`$ `5R` $`6`$-CP + `10` $`\times`$ `5R`/`6R`/`7R`/`8R` full statistical benchmark), `verify_higher_rounds.log` (`6R` distinguisher + `9R`–`12R` audit), `attack_9r_splitmix64_demo.log` (`9R` inner kernel check at $`D = 2^{22.32}`$ CP), `attack_9r_20trials.log` (`20`-trial `9R` inner kernel check), and `distinguishers_6r_nullspace_audit.log` (`160`-key empirical nullity audit of $`Q_5, Q_6, Q_7, Q_8`$ at $`d \in \{16, 20, 24\}`$). |
