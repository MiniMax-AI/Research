# Reproducible Key-Recovery Artifacts for HAWK-256

Single-core experimental suite verifying our end-to-end key recovery against `HAWK-256` ($`n = 256, \ell = 9`$, *ASIACRYPT 2022*) from the 450-byte reference public key alone (Sections 4–7 of [../paper/main.pdf](../paper/main.pdf)). Starting from the $`\tau`$-cocycle framework (*ePrint 2026/1593*), the suite constructs the index-$`2^{255}`$ modulo-two sublattice $`\Lambda' = 2\Lambda_\tau + \mathbb{Z}w_0`$, applies exact high-precision Cholesky basis preconditioning and `flatter` reduction, recovers the unique shortest cocycle vector $`\pm V_\tau`$ via progressive single-thread `g6k` sieving at maximum dimension $`d = 80`$, descends the cyclotomic $`\tau`$-tower to an equivalent secret basis $`B_{\mathrm{rec}}`$ satisfying $`B_{\mathrm{rec}}^* B_{\mathrm{rec}} = Q`$, and verifies the recovered 592-byte secret key directly against the unmodified HAWK NIST reference C implementation.

---

## Environment and Dependencies

### Supported Platform and Architecture
- **Supported Platforms**: Linux and macOS (`x86_64` and `arm64` / `aarch64`); reference single-core timing benchmarks were measured on Linux `x86_64`.
  - `kat/build_shim.sh` compiles `kat/hawk_shim.c` on the fly into a temporary scratch binary against the unmodified official HAWK NIST submission portable C reference implementation in `kat/hawk_ref/` (`hawk-sign/dev` @ `1b9fef5`), providing native C `sign`/`verify` verification across Linux and macOS on both `x86_64` and `arm64`.
  - All algebraic certificate, sublattice construction, Cholesky preconditioning, sieving, and $`\tau`$-tower descent scripts in `chain/` and `hawk_recover/` are platform-independent SageMath/Python.
- **Single-Core Benchmark Configuration**: All scripts enforce strictly single-threaded execution (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `threads=1` in `fpylll`/`g6k`) using an isolated temporary working directory (`HAWK_WORK`). Peak memory usage is $`\approx 8\text{ GiB}`$ during `--full` (HNF kernel extraction) and $`< 1\text{ GiB}`$ during `--demo`.

### Required Software (All Modes)
- **C Compiler** (`cc`, `gcc`, or `clang`, configurable via `CC`): compiles the portable C `sign`/`verify` oracle (`kat/hawk_shim.c` + `kat/hawk_ref/*.c`).
- **SageMath 10.x** (`sage` on `$PATH`, or `PY=/path/to/sage-python`): provides `sage.all`, `numpy`, `scipy`, `fpylll` ($`0.6.\text{x}`$), and `cypari2`/PARI. Used by all modes (`--demo`, `--kat`, and `--full`).

### Full-Pipeline Lattice Tools (`--full` Mode)
- **`flatter`** on `$PATH` (*CRYPTO 2023*, required for `--full`): used in Stage 2 to compress the $`5{,}483`$-bit HNF integer kernel of $`\Lambda_\tau`$ to a $`12`$-bit basis, and in Stage 4 to reduce the $`s = 61`$ Cholesky basis $`B_{\mathrm{chol}}`$ before exact unimodular transform recovery $`G_1 = U G_p U^T`$.
  - *Build / install*: requires CMake, a C++17 compiler, `gmp`, `mpfr`, `eigen3`, and `fplll`:
    ```bash
    git clone https://github.com/keeganryan/flatter.git
    cmake -S flatter -B flatter/build -DCMAKE_BUILD_TYPE=Release
    cmake --build flatter/build -j
    sudo cmake --install flatter/build
    ```
- **`g6k` 0.1.2** (General Sieve Kernel, *EUROCRYPT 2019*, optional for `--full`): used in Stage 5 (`chain/g6k_tour.py`) to execute the single-thread progressive `pump_n_jump_bkz_tour` across sieve dimensions $`d \in \{40, 50, 60, 70, 80\}`$ on the $`s = 21`$ Cholesky embedding `Bs.txt`.
  - *Install into SageMath's Python*:
    ```bash
    git clone --branch v0.1.2 https://github.com/fplll/g6k.git
    cd g6k && "$(sage -c 'import sys; print(sys.executable)')" -m pip install -e .
    ```
    (Note: `g6k` emits a benign `MAX_SIEVING_DIM` warning because ambient lattice rank is $`256`$, whereas active sieving never exceeds $`d = 80`$.)
  - *Automatic fallback when `g6k` is absent*: if `import g6k` is not available in SageMath's Python environment, `./run.sh --full` still rebuilds $`\Lambda_\tau`$, $`\Lambda'`$, and the Cholesky-preconditioned Gram matrix $`G_1`$ from `data/pk.bin` via `flatter`, announces the fallback, reuses the precomputed reduced basis `data/Bred.sobj`, verifies all algebraic certificates, runs the full $`\tau`$-tower descent ($`\ell = 9 \to 3`$) to generate a fresh 592-byte secret key, and verifies it against `data/pk.bin`.

---

## Quick Start / Usage

On first run, `./run.sh`, `kat/build_shim.sh`, and `chain/*.py` automatically invoke `./fetch_third_party.sh --ensure` (notifying the user on `stderr`) to download and SHA-256-verify the external HAWK C reference (`kat/hawk_ref/`), `hawk_recover/` package, and AGPS20 classical sieve cost table (`agps20/`); all subsequent runs execute completely offline and single-threaded (`OMP_NUM_THREADS=1`) in a temporary scratch directory cleaned up automatically on exit:

```bash
chmod +x run.sh
./fetch_third_party.sh   # Optional: explicitly pre-fetch / verify third-party reference sources

./run.sh --demo    # ~4 s   (default; alias: --kat): verifies all 4 exact algebraic sublattice,
                   #        Cholesky pullback, and cocycle certificates on data/, verifies GSA block
                   #        sizes & AGPS20 gate complexities across HAWK-256/512/1024, then verifies
                   #        data/recovered_skd.bin against data/pk.bin via C sign/verify (kat/)
./run.sh --kat     # ~4 s   : alias for --demo (fast certificate + security estimator + C sign/verify check)
./run.sh --full    # ~31 min: rebuilds Lambda_tau & Lambda' from data/pk.bin (HNF + flatter + Cholesky),
                   #        runs single-thread g6k d = 80 SVP tour (~7 min total if falling back to
                   #        data/Bred.sobj when g6k is absent), verifies certificates & security estimates,
                   #        runs tau-tower descent (ell = 9 -> 3), encodes 592B key, and verifies via C sign/verify
./run.sh --help    # prints CLI usage summary
```

To run the security and AGPS20 gate-complexity estimator (including the DDGR20 / HAWK BKZ simulator) or the cyclotomic $`\tau`$-tower descent standalone:

```bash
"$(sage -c 'import sys; print(sys.executable)')" -B -u chain/estimate_security.py --simulator
HAWK_WORK="$(mktemp -d)" "$(sage -c 'import sys; print(sys.executable)')" -B -u chain/recover.py
```

---

## Key Files

| Path | Description |
| :--- | :--- |
| `run.sh` | Unified single-core verification entry point (`--demo`, `--kat`, `--full`, `--help`; auto-fetches external reference files on first run). |
| `fetch_third_party.sh` | Downloads and SHA-256-verifies `kat/hawk_ref/` (`hawk-sign/dev` @ `1b9fef52`), `kat/hawk_shim.c`, `hawk_recover/` (`anthropics/cryptography-research-demo` @ `fa01c398`), and `agps20/cost-estimate-list_decoding-classical.csv` (`jschanck/eprint-2019-1161` @ `1ec49e5b`) into gitignored directories. |
| `chain/common.py` | Path-independent matrix I/O (`.gz`-transparent), automatic `fetch_third_party.sh --ensure` hook, and POSIX child-inclusive CPU-time accounting. |
| `chain/estimate_security.py` | Concrete security and AGPS20 gate-complexity estimator (Sections 3.5, 4.1, 6.2–6.3, Table 4): verifies GSA uSVP block sizes $`\beta(\Lambda_\tau)=(92, 205, 432)`$ and $`\beta(\Lambda')=(58, 156, 351)`$, effective sieve dimensions $`\beta_{\mathrm{eff}}`$, classical list-decoding gate counts (exact spherical-wedge evaluation at $`d=40`$ + tabulated CSV lookup for $`d \ge 64`$), provable $`\mathrm{SVP}_{n/2+1}`$ bounds, and (`--simulator`) DDGR20 BKZ simulator block sizes. |
| `chain/build_lambda.py` | Stages 1–4 (Section 5): 450-byte public-key decode, $`1024 \times 512`$ constraint matrix assembly, HNF + `flatter` saturated kernel extraction for $`\Lambda_\tau`$, $`\mathbb{F}_2`$ parity solve for $`\Lambda'`$, exact $`s = 61`$ Cholesky preconditioning $`G_1 = U G_p U^T`$, and $`s = 21`$ sieve embedding (`Bs.txt`). |
| `chain/g6k_tour.py` | Stage 5 (Section 6.1): single-thread (`threads=1`) `g6k` progressive `pump_n_jump_bkz_tour` across sieve dimensions $`d \in \{40, 50, 60, 70, 80\}`$. |
| `chain/certificates.py` | Stage 6 (Sections 4.4 & 6.2): verifies $`\log_2\det\mathrm{Gram}(\Lambda') = 2941.00`$ ($`[\Lambda_\tau:\Lambda']=2^{255}`$), shortest-row norm ratio $`0.999625`$, exact triangular pullback $`v = uU`$ ($`u G_1 u^T = v G_p v^T = 512`$), $`V_\tau \equiv I_2 \pmod 2`$, and $`V_\tau \tau(V_\tau) = I_2`$; also provides `--verify-sk` for exact algebraic verification of 592-byte secret keys. |
| `chain/recover.py` | Stages 7–8 (Section 6.3): $`\tau`$-tower descent ($`\ell = 9 \to 3`$) from $`V_\tau`$ to $`B_{\mathrm{rec}}`$ ($`B_{\mathrm{rec}}^* B_{\mathrm{rec}} = Q`$) and encoding to the 592-byte reference secret-key format. |
| `hawk_recover/` *(fetched on demand)* | Reference implementation of the $`\tau`$-cocycle and $`\tau`$-tower descent (*ePrint 2026/1593*, `anthropics/cryptography-research-demo` @ `fa01c398`), with `infra/shim.py` wired to `kat/build_shim.sh` for cross-platform builds. |
| `kat/` | Portable C `sign`/`verify` build wrapper (`build_shim.sh`); compiles `hawk_shim.c` and the unmodified HAWK NIST submission C reference (`hawk_ref/`, fetched on demand from `hawk-sign/dev` @ `1b9fef52`). |
| `data/` | Genuine 450-byte reference public key (`pk.bin`), original secret key (`skd.bin`, never read by the attack), recovered secret key (`recovered_skd.bin`), preconditioned Gram/Cholesky matrices (`Gp_k1.txt.gz`, `G1_k1.txt.gz`, `Bc.txt.gz`, `Bcr.txt.gz`, `Bs.txt.gz`, `Kflat.txt.gz`), reduced basis (`Bred.sobj`), and $`V_\tau`$ coordinates (`vtau_coords.sobj`). |
| `logs/` | Recorded single-core execution logs for Stages 1–4 (`stage1_4_precond.log`, $`204.5\text{ s}`$), Stage 5 (`stage5_g6k_svp.log`, $`1{,}469.5\text{ s}`$), and Stages 6–8 (`stage6_8_tower.log`, $`177.8\text{ s}`$). |

---

## Measured Single-Core Stage Breakdown (`1,851.8` CPU-Core-Seconds Total)

| Stage | Operation | Single-Core Time | Share |
| :--- | :--- | ---: | ---: |
| Stage 1 | Decode 450-byte `pk.bin` $`\to Q`$ and assemble $`1024 \times 512`$ constraint matrix | $`0.3\text{ s}`$ | $`<0.1\%`$ |
| Stage 2 | Extract saturated $`\mathbb{Z}`$-basis of $`\Lambda_\tau`$ (HNF $`17.1\text{ s}`$ + `flatter` $`125.9\text{ s}`$) | $`143.0\text{ s}`$ | $`7.7\%`$ |
| Stage 3 | Compute $`\mathrm{Gram}(\Lambda_\tau)`$, solve $`\mathbb{F}_2`$ parity for $`\Lambda'`$, form $`G_p = \mathrm{Gram}(\Lambda')`$ | $`1.9\text{ s}`$ | $`0.1\%`$ |
| Stage 4 | High-precision Cholesky $`G_p = LL^T`$, `flatter` on $`B_{\mathrm{chol}}`$, pullback $`G_1 = U G_p U^T`$ | $`59.3\text{ s}`$ | $`3.2\%`$ |
| Subtotal | All pre-SVP setup and preconditioning stages (Stages 1–4) | $`204.5\text{ s}`$ | $`11.0\%`$ |
| Stage 5 | Single-thread `g6k` progressive sieve tour ($`d \in \{40, 50, 60, 70, 80\}`$) | $`1{,}469.5\text{ s}`$ | $`79.4\%`$ |
| Stages 6–7 | Exact coordinate pullback $`v = uU`$ and $`\tau`$-tower descent ($`\ell = 9 \to 3`$) | $`174.2\text{ s}`$ | $`9.4\%`$ |
| Stage 8 | Encode 592-byte secret key and run C reference `sign`/`verify` check | $`3.6\text{ s}`$ | $`0.2\%`$ |
| Total | Complete end-to-end HAWK-256 key recovery ($`690\times`$ vs. $`1{,}278{,}720\text{ s}`$ baseline) | $`1{,}851.8\text{ s}`$ | $`100.0\%`$ |

---

## External Reference Sources (Fetched on Demand)

Third-party reference implementations are not redistributed in this repository; `./fetch_third_party.sh` (invoked automatically when missing) downloads and SHA-256-verifies them locally into gitignored directories:

- `hawk_recover/` and `kat/hawk_shim.c` — Fetched from `anthropics/cryptography-research-demo` @ `fa01c398` (Apache-2.0).
- `kat/hawk_ref/` — Fetched from `hawk-sign/dev` @ `1b9fef52` (MIT).
- `agps20/cost-estimate-list_decoding-classical.csv` — Fetched from `jschanck/eprint-2019-1161` @ `1ec49e5b` (MIT).
