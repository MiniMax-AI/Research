# Correlated $`\mathsf{VOLE}`$-Masking Protocol Repair (`protocol_repair`)

This directory contains the C implementation and experimental verification for Section 6 (Theorem 3 and Table 2) of the paper ([`../../paper/main.pdf`](../../paper/main.pdf)).

## Overview

`make_masked.py` applies our 5-step correlated $`\mathsf{VOLE}`$-masking repair (`-DMASKB`) to a clean unpacked copy of the official ASIACRYPT 2025 C reference implementation (`SBCvole.c` and `SBCunified.h`, fetched on demand via `../fetch_third_party.sh`):
1. Expands each $`\mathsf{cGGM}`$ leaf with an additional field mask share $`\rho_j^{[\![\ell]\!]} \in \mathbb{F}_{2^{257}}`$.
2. Defines the global mask $`\rho \in \mathbb{F}_{2^{257}}`$ from tree $`T_{\vec{x}_0}`$ and transmits $`(\tau - 1)`$ per-tree alignment offsets $`(\delta_{\rho_1}, \dots, \delta_{\rho_{\tau-1}}) \in (\mathbb{F}_{2^{257}})^{\tau-1}`$.
3. Binds $`(\delta_{\rho_j})_{j \in [\tau]^*}`$ inside the first Fiat–Shamir hash $`H_0`$ prior to challenge derivation.
4. Masks both verification coefficients as $`\tilde{A} = A \oplus A_{F_\rho}`$ (hashed into $`H_1`$) and $`\tilde{B} = B \oplus \rho`$ (transmitted).
5. Updates the verifier to reconstruct $`F_\rho(\Delta^*) = A_{F_\rho} \oplus \Delta^* \rho`$ and check $`\tilde{A}' = \Phi(\Delta^*) \oplus F_\rho(\Delta^*) \oplus \tilde{B}\Delta^*`$.

In addition to verifying exact byte-size overheads (Table 2) and zero single- or multi-signature key leakage ($`m \in \{1, 2, 4, 16\}`$), `run.sh` compiles four soundness-control binaries (`forge_noh0.inc`) confirming that omitting $`(\delta_{\rho_j})_{j \in [\tau]^*}`$ from $`H_0`$ enables universal key-only forgery whereas binding it in $`H_0`$ rejects false witnesses.

## Requirements

- **C Compiler & Utilities**: `gcc` or `clang` (override via `CC`), `bash`, and `sha256sum` / `shasum` (`run.sh` and `make_masked.py` automatically invoke `../fetch_third_party.sh --ensure` on first run to fetch `vendor/asiacrypt-2025-a10.zip` and `crypto-algorithms/`, notifying the user on `stderr`; zero external C libraries are required).
- **Python**: `Python 3.8+` (standard library only; zero third-party packages required; subsequent runs after first fetch are completely offline).
- **Platforms**: 64-bit Linux and macOS (`x86_64` and `arm64`).

## Usage

```bash
./run.sh --demo   # Verifies repair, zero leakage, & H_0 binding on (D, tau) = (9, 15) (~15 s)
./run.sh --full   # Verifies repair, zero leakage, & H_0 binding on headline (D, tau) = (13, 10)
```

## Signature-Size Comparison Across All Parameter Sets

Measured directly on the official C reference implementation (`q = 2, k = 257, n = 130, lambda = 128`), where the exact bit-packed overhead is `delta_bytes(D, tau) = 32*(tau - 1) + ceil(2*tau/8) - ceil((tau + 1)/8)`:

| `(D, tau)` | `SBC_VOLE` (unmodified) | `SBC_VOLE` (repaired) | Overhead | `SBC_VOLE^M` (unmodified) | `SBC_VOLE^M` (repaired) | Overhead |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `(9, 15)` | `3458 B` | `3908 B` | `+450 B` | `2978 B` | `3428 B` | `+450 B` |
| `(10, 13)` | `3218 B` | `3604 B` | `+386 B` | `2802 B` | `3188 B` | `+386 B` |
| `(11, 12)` | `3170 B` | `3523 B` | `+353 B` | `2786 B` | `3139 B` | `+353 B` |
| `(12, 11)` | `3090 B` | `3411 B` | `+321 B` | `2738 B` | `3059 B` | `+321 B` |
| **`(13, 10)`** | **`2978 B`** | **`3267 B`** | **`+289 B`** | **`2658 B`** | **`2947 B`** | **`+289 B`** |
| `(15, 9)` | `2978 B` | `3235 B` | `+257 B` | `2690 B` | `2947 B` | `+257 B` |
