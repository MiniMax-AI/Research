# Verification Artifacts for Single-Signature Key Recovery on $`\mathsf{SBC}_{\mathsf{VOLE}}`$

This directory contains the self-contained cryptanalysis, simulator-distinguisher, and protocol-repair verification suites accompanying the paper ([`../paper/main.pdf`](../paper/main.pdf)). All experiments verify against both standalone transcriptions of the specification and the official ASIACRYPT 2025 C reference implementation (`asiacrypt/2025/a10`, fetched and SHA-256-verified on demand via `./fetch_third_party.sh`) across all six published $`(D, \tau)`$ parameter sets and both commitment variants ($`\mathsf{SBC}_{\mathsf{VOLE}}`$ and $`\mathsf{SBC}_{\mathsf{VOLE}}^{\mathsf{M}}`$).

## Environment and Dependencies

- **Supported Platforms & Architecture**: Portable C and pure Python 3; runs out-of-the-box on 64-bit Linux and macOS (`x86_64` and `arm64`). On first run, `./run.sh` (as well as sub-suite `run.sh` / `build.sh` / `make_masked.py`) automatically invokes `./fetch_third_party.sh --ensure` to download and SHA-256-verify the external IACR reference archive and portable AES/SHA-256 routines (printing a terminal notice); all subsequent runs execute completely offline. Hardware-specific x86 AES-NI intrinsics in the upstream reference code are bypassed automatically via `-DGENERICAES -DGENERICSHA`.
- **C Compiler & System Utilities**: Standard C compiler (`gcc` or `clang`, configurable via `CC`), `bash`, `curl` (optional; falls back to Python `urllib`), and `sha256sum` / `shasum` (used by `key_recovery/` and `protocol_repair/` to build and verify against the ASIACRYPT 2025 C reference implementation; portable AES and SHA-256 routines are fetched into `vendor/crypto-algorithms/` and `protocol_repair/crypto-algorithms/`, requiring zero external C libraries).
- **Python Environment**: `Python 3.8+` (standard library only; zero third-party packages required).

## Directory Layout

| Directory / File | Paper Reference | Contents |
| :--- | :--- | :--- |
| [`run.sh`](run.sh) | All sections | Top-level verification runner executing all three sub-suites sequentially (`--demo` or `--full`). |
| [`fetch_third_party.sh`](fetch_third_party.sh) | Infrastructure | Downloads and SHA-256-verifies `asiacrypt-2025-a10.zip` and `B-Con/crypto-algorithms` into gitignored paths on demand. |
| [`key_recovery/`](key_recovery/) | Section 4 & Section 7 (Lemma 1, Theorem 1, Algorithm 1, Table 3) | Concrete bit-operation cost model, closed-form $`\mathbb{F}_2`$-affine identity verification, zero-patch C transcript extractor (`forge.c`), single-signature $`257 \times 256`$ $`\mathbb{F}_2`$ key-recovery solver, universal signature forgery, and $`1{,}200`$-instance sweep across all $`12`$ configurations. |
| [`hvzk_distinguisher/`](hvzk_distinguisher/) | Section 4.2 & Section 5 (Lemma 1, Theorem 2) | Pure-Python verification of characteristic-$`2`$ quadratic cross-term cancellation under both padding conventions ($`\eta \in \{0, 1\}`$), zero conditional entropy of $`B`$, and the $`1`$-query Game 2 simulator distinguisher $`\mathcal{D}`$. |
| [`protocol_repair/`](protocol_repair/) | Section 6 (Theorem 3, Table 2) | Correlated $`\mathsf{VOLE}`$ one-time-pad protocol repair (`make_masked.py`), exact byte-overhead verification against Table 2, multi-signature pooled rank-deficiency tests ($`m \in \{1, 2, 4, 16\}`$), and pre-challenge $`H_0`$-binding soundness controls (`forge_noh0.inc`). |

## Quick Start / Usage

### Run All Verification Suites

From `artifacts/`, run all three verification suites end-to-end:

```bash
./run.sh --demo        # Fast end-to-end verification across all 3 suites (~40 s)
./run.sh --full        # Extended verification including headline (D, tau) = (13, 10) sets (~90 s)
```

Each suite prints `RESULT: SUCCESS` and exits with status `0` iff every check passes.

### 1. Single-Signature Key Recovery (`key_recovery/`)

```bash
cd key_recovery
./run.sh --demo        # (D, tau) = (9, 15) basic + MinMax, cost model, and Fig. 5 checks (~18 s)
./run.sh --full        # Adds headline (D, tau) = (13, 10) 2,978 B / 2,658 B sets (~36 s)
./run.sh --sweep 100   # Full 1,200-instance sweep: 6 (D, tau) sets x 2 variants x 100 keys
```

Standalone algebraic and cost verifiers can also be invoked directly:

```bash
python3 cost.py
python3 verify_lemma1.py
python3 structural.py
python3 sbcvole.py --demo
```

### 2. $`\mathsf{HVZK}`$ & Game 2 Simulator Distinguisher (`hvzk_distinguisher/`)

```bash
cd hvzk_distinguisher
./run.sh --demo        # 30 instances per padding convention & simulator strategy (~6 s)
./run.sh --full        # 120 instances per padding convention & simulator strategy (~25 s)
```

### 3. Correlated $`\mathsf{VOLE}`$-Masking Protocol Repair (`protocol_repair/`)

```bash
cd protocol_repair
./run.sh --demo        # Repair, zero-leakage, & H_0-binding checks at (D, tau) = (9, 15) (~15 s)
./run.sh --full        # Same checks on headline (D, tau) = (13, 10) 2,978 B / 2,658 B sets
```

See [`protocol_repair/README.md`](protocol_repair/README.md) (and Table 2 of [`../paper/main.pdf`](../paper/main.pdf)) for the complete signature-size comparison across all $`12`$ configurations.

## External Reference Sources (Fetched on Demand)

Third-party reference implementations are not redistributed in this repository; `./fetch_third_party.sh` (invoked automatically when missing) downloads and SHA-256-verifies them locally into gitignored directories:

- `asiacrypt-2025-a10.zip` (`https://artifacts.iacr.org/files/asiacrypt/2025/asiacrypt-2025-a10.zip`, SHA-256 `e275d529...`, GNU AGPL-3.0) — Official ASIACRYPT 2025 `SBCvole` C artifact, placed into `key_recovery/vendor/` and `protocol_repair/vendor/`.
- `aes.c`, `aes.h`, `sha256.c`, `sha256.h` (`B-Con/crypto-algorithms` @ `cfbde484`, Public Domain) — Portable AES and SHA-256 implementations for `-DGENERICAES -DGENERICSHA`, placed into `key_recovery/vendor/crypto-algorithms/` and `protocol_repair/crypto-algorithms/`.
