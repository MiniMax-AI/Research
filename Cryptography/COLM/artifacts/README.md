# Verification Artifacts: Two-Query Nonce-Respecting Forgeries on $`\mathsf{COLM}_0`$

This directory contains the self-contained experimental and algebraic verification suite accompanying the paper [Two-Query Nonce-Respecting Forgeries on $`\mathsf{COLM}_0`$](../paper/main.pdf).

The artifacts verify all main technical results of the paper:
1. **End-to-End Two-Query Encryption-Only Splice Forgeries (`splice_forgery/`):** Demonstrates deterministic (success probability $`1`$) ciphertext-splicing forgeries using $`2`$ encryption queries and $`0`$ decryption queries under distinct nonces ($`\mathsf{INT\text{-}CTXT}`$) against both the unmodified 2016 CAESAR SUPERCOP implementations (`colm0v1/ref`, `colm0v1/aesni4`, `colm0v1/aesni8`) and the 2024 *Journal of Cryptology* specification of $`\mathsf{AES\text{-}COLM}_0\text{-}128`$ and $`\mathsf{AES\text{-}COLM}_0\text{-}256`$, together with length-floor boundary sweeps, $`\mathsf{COLM}_{127}`$ negative controls, and one-query baseline checks.
2. **Algebraic Derivations and Security-Proof Audit (`algebra_and_proof_audit/`):** Pure Python 3 scripts deriving the splice index sets $`S_{\mathrm{fade}} = \{1, 2, 3, 8, 129\}`$ and $`S_{\mathrm{trans}}(\ell, m) = (\mathrm{supp}(3^{-1}x^{m-1}) + 1) \cup \{m, \dots, \ell+1\}`$ ($`129 \le m \le \ell`$, e.g. $`\{1, 3, 4, 5, 6, 7\} \cup \{129, \dots, \ell+1\}`$ at $`m = 129`$) directly from $`\mathbb{F}_{2^{128}} = \mathbb{F}_2[x]/(x^{128} + x^7 + x^2 + x + 1)`$, proving the tight length floors ($`\ell_{\min} = 130`$ and $`\ell_{\min} = 129`$) and the $`\mathsf{COLM}_\tau`$ ($`1 \le \tau \le 127`$) degree obstruction, and evaluating Theorem 2 of both the 2024 *Journal of Cryptology* proof ($`2^{-105.92}`$) and the 2016 CAESAR Round-3 security proof addendum ($`2^{-107.56}`$) at the exact adversary parameters.

---

## Environment and Dependencies

- **Supported Platforms & Architecture:**
  - **Portable C & Python suites (`x86_64` and `arm64`):** All portable C targets—`colm0v1/ref`, `colm127v1/ref`, the *Journal of Cryptology* (2024) implementations ($`\mathsf{AES\text{-}COLM}_0\text{-}128`$ and $`\mathsf{AES\text{-}COLM}_0\text{-}256`$), cross-implementation equivalence checks (`xcheck`), and one-query baselines—as well as all Python scripts run out of the box on any 64-bit Linux or macOS machine (`x86_64` or `arm64` / Apple Silicon).
  - **Optional `x86_64` hardware-intrinsic cross-checks:** Cross-verification against the parallel SUPERCOP implementations (`colm0v1/aesni4` and `colm0v1/aesni8`) requires an `x86_64` CPU with `AES-NI` and `AVX2` support. On `arm64` hosts (or `x86_64` CPUs without `AES-NI`/`AVX2`), `splice_forgery/run.sh` automatically skips those two hardware-specific builds with an explicit status message while executing all portable targets. Symbol isolation for `xcheck` uses standard C preprocessor renaming (`-Dcrypto_aead_encrypt=...`), so GNU `objcopy` / `binutils` is not required on any platform.
- **C Compiler & OpenSSL:**
  - **Compiler:** Any standard C99/C11 compiler (`gcc` or `clang`) with standard `libc`/`libm` (honors `CC` if set; defaults to `gcc` or `cc`).
  - **OpenSSL / AES primitive:** `splice_forgery/openssl/aes.h` provides a self-contained FIPS-197 `AES-128`/`AES-256` header shim implementing the OpenSSL `AES_KEY` interface (`AES_set_encrypt_key`, `AES_set_decrypt_key`, `AES_encrypt`, `AES_decrypt`), validated at runtime against FIPS-197 Appendix C.1 by `aes_kat.c`. As a result, `./run.sh` builds out of the box with zero external library dependencies. If compiling the upstream SUPERCOP sources standalone against system OpenSSL `libcrypto` outside this harness, install OpenSSL development headers (`sudo apt install libssl-dev` on Ubuntu/Debian; `brew install openssl@3` on macOS).
- **Python Environment:**
  - `Python 3.8+` using the standard library only (zero external `pip` packages or computer-algebra dependencies) for `algebra_and_proof_audit/`.

---

## Directory Layout

```text
artifacts/
├── README.md                        # This overview and reproducibility guide
├── run.sh                           # Top-level runner executing both verification suites
├── fetch_third_party.sh             # Streams & SHA-256-verifies CAESAR SUPERCOP COLM sources on demand
├── splice_forgery/                  # End-to-end C splice forgery & cross-implementation suite
│   ├── README.md                    # Documentation for the four splice variants
│   ├── run.sh                       # Build & test runner (--demo and --full modes)
│   ├── splice.c                     # Consolidated driver deriving splice sets over GF(2^128) at runtime
│   ├── joc2024_colm0.c              # Direct C implementation of J. Cryptol. 2024 Fig. 2 & Fig. 3 (128/256-bit)
│   ├── xcheck.c                     # Cross-implementation ciphertext & forgery equivalence checker
│   ├── onequery.c                   # Baseline check rejecting 8 single-query manipulation strategies
│   ├── aes_kat.c                    # FIPS-197 Appendix C.1 Known-Answer Test for AES-128
│   ├── crypto_aead.h                # SUPERCOP crypto_aead API declarations
│   ├── openssl/aes.h                # Self-contained FIPS-197 AES-128 & AES-256 header shim
│   ├── colm0v1_ref/                 # (Fetched on demand) CAESAR SUPERCOP colm0v1/ref implementation
│   ├── colm0v1_aesni4/              # (Fetched on demand) CAESAR SUPERCOP colm0v1/aesni4 (4-way AES-NI/AVX2)
│   ├── colm0v1_aesni8/              # (Fetched on demand) CAESAR SUPERCOP colm0v1/aesni8 (8-way AES-NI/AVX2)
│   ├── colm127v1_ref/               # (Fetched on demand) CAESAR SUPERCOP colm127v1/ref (intermediate-tag control)
│   ├── witnesses/                   # Archived binary test vectors & forged ciphertexts (trial 0)
│   └── run_demo.log                 # Archived reference log of ./run.sh --demo on x86_64 w/ AVX2
└── algebra_and_proof_audit/         # Pure Python 3 algebraic derivations & security-proof audit
    ├── README.md                    # Documentation for the Python verification scripts
    ├── run.sh                       # Runner executing all three Python verification scripts
    ├── algebra.py                   # GF(2^128) splice sets, rank/length-floor proofs & COLM_127 obstruction
    ├── derivation.py                # Case A / B / C differential equations & exhaustive GF(2^m) checks
    ├── thm2.py                      # Quantitative evaluation of J. Cryptol. 2024 & 2016 Addendum Theorem 2
    └── run_demo.log                 # Archived reference output log
```

---

## Quick Start / Usage

All scripts compile temporary binaries into an ephemeral `mktemp -d` directory cleaned up automatically on exit, leaving the working tree clean. On the first run, `./run.sh` (or `splice_forgery/run.sh`) automatically calls `./fetch_third_party.sh --ensure` to stream and SHA-256-verify the CAESAR SUPERCOP `colm0v1` / `colm127v1` sources (~200 KB extracted, streamed early from the SUPERCOP archive) with a terminal notice; subsequent runs execute completely offline.

### Run Everything (Top-Level Runner)

From `artifacts/`:

```bash
./fetch_third_party.sh   # Optional: explicitly pre-fetch / verify SUPERCOP COLM files
./run.sh                 # Fast verification (~10-15 s; runs algebra_and_proof_audit + splice_forgery --demo)
./run.sh --full          # Full verification (~2 min; runs algebra_and_proof_audit + splice_forgery --full)
```

### 1. End-to-End Two-Query Splice Forgeries (`splice_forgery/`)

```bash
cd splice_forgery
./run.sh --demo   # Fast verification: 25 trials/variant + 16 boundary sweeps + 200 random controls (~10 s)
./run.sh --full   # Full verification: 1,000 trials/variant + 16 boundary sweeps + 2,000 random controls (~2 min)
```

`splice_forgery/run.sh` executes seven verification stages and exits with `RESULT: SUCCESS` (exit code `0`) iff every assertion holds:
1. Builds portable binaries (`splice0`, `splice127`, `joc128`, `joc256`, `aes_kat`, `onequery`).
2. Runs the FIPS-197 Appendix C.1 AES-128 Known-Answer Test (`aes_kat`).
3. Verifies $`100\%`$ forgery acceptance (and $`0\%`$ acceptance across drop-one, add-one, and random subset negative controls) for all four variants (`fade`, `transplant`, `misuse`, `prefix`) on SUPERCOP `colm0v1/ref` and on the *Journal of Cryptology* (2024) specifications $`\mathsf{AES\text{-}COLM}_0\text{-}128`$ and $`\mathsf{AES\text{-}COLM}_0\text{-}256`$.
4. Sweeps $`16`$ boundary configurations confirming exact rejection at $`\ell = 129`$ (`fade`), $`\ell = 128`$ (`transplant`), $`\ell = 130`$ (`misuse` with $`s = 1`$), and $`128`$-block shared prefix (`prefix`), as well as acceptance under empty/unaligned associated data and incomplete final blocks.
5. Verifies as a control that $`\mathsf{COLM}_{127}`$ (`colm127v1/ref`, stride $`\tau = 127`$) rejects $`100\%`$ of splice attempts.
6. Cross-checks byte-for-byte ciphertext identity and joint forgery acceptance between `colm0v1/ref` and *Journal of Cryptology* (2024) $`\mathsf{AES\text{-}COLM}_0\text{-}128`$ (plus `colm0v1/aesni4` and `colm0v1/aesni8` when run on `x86_64` with `AES-NI`/`AVX2`).
7. Verifies that all $`8`$ single-query ciphertext manipulation strategies in `onequery.c` are rejected ($`0`$ acceptances).

### 2. $`\mathbb{F}_{2^{128}}`$ Algebraic Derivations and Proof Audit (`algebra_and_proof_audit/`)

```bash
cd algebra_and_proof_audit
./run.sh
```

`algebra_and_proof_audit/run.sh` runs `algebra.py`, `derivation.py`, and `thm2.py` (`python3 -B`) and prints `RESULT: SUCCESS` (exit code `0`) after verifying all algebraic identities, subspace ranks, small-field ($`m \in \{4, \dots, 12\}`$) exhaustive subset checks, and security-proof bound evaluations.

---

## External Reference Sources (Fetched on Demand)

Third-party CAESAR SUPERCOP reference implementations are not redistributed in this repository; `./fetch_third_party.sh` (invoked automatically when missing) streams and SHA-256-verifies `crypto_aead/colm0v1/{ref,aesni4,aesni8}` and `crypto_aead/colm127v1/ref` from SUPERCOP (`https://bench.cr.yp.to/supercop.html`, Public Domain) into gitignored `splice_forgery/colm*` directories.
