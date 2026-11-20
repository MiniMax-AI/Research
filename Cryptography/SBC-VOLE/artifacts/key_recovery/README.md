# Single-Signature Key Recovery Artifact (`key_recovery`)

This directory contains the self-contained verification suite for Section 4 (Lemma 1, Theorem 1, and Algorithm 1) and Section 7 (Table 3) of the paper ([`../../paper/main.pdf`](../../paper/main.pdf)):
- Verification of the characteristic-$`2`$ quadratic cross-term cancellation and closed-form $`\mathbb{F}_2`$-affine identity for $`B`$ (Lemma 1).
- Concrete bit-operation cost evaluation under headline, conservative, and Karatsuba models (`cost.py`, Theorem 1).
- End-to-end single-signature secret-key recovery and universal signature forgery against the official ASIACRYPT 2025 C reference implementation (`asiacrypt/2025/a10`, fetched on demand via `../fetch_third_party.sh`) across all six published $`(D, \tau)`$ parameter sets and both variants ($`\mathsf{SBC}_{\mathsf{VOLE}}`$ and $`\mathsf{SBC}_{\mathsf{VOLE}}^{\mathsf{M}}`$).

## Requirements

- **C Compiler & Utilities**: `gcc` or `clang` (override via `CC`), `bash`, and `sha256sum` / `shasum` (`run.sh` and `build.sh` automatically invoke `../fetch_third_party.sh --ensure` on first run to fetch `vendor/asiacrypt-2025-a10.zip` and `vendor/crypto-algorithms/`, notifying the user on `stderr`; zero external C libraries are required).
- **Python**: `Python 3.8+` (standard library only; zero third-party packages required; subsequent runs after first fetch are completely offline).
- **Platforms**: 64-bit Linux and macOS (`x86_64` and `arm64`).

## Usage

```bash
./run.sh --demo        # (D, tau) = (9, 15) basic + MinMax, plus printed-spec keys (~18 s)
./run.sh --full        # Adds headline (D, tau) = (13, 10) 2,978 B / 2,658 B sets (~36 s)
./run.sh --sweep 100   # All 12 published (scheme x parameter) configurations x 100 keys (1,200 instances)
```

You can also run the standalone formula, reference-protocol, and cost verifiers directly:

```bash
python3 cost.py
python3 verify_lemma1.py
python3 structural.py
python3 sbcvole.py --demo
```

## Key Files

| File | Description |
| :--- | :--- |
| `run.sh` | Main verification driver (`--demo`, `--full`, `--sweep`). |
| `build.sh` | Unpacks `vendor/asiacrypt-2025-a10.zip` into a temporary directory and compiles unmodified `SBCvole.c` alongside `forge.c`. |
| `cost.py` | Evaluates the bit-operation complexity formula (Theorem 1) under headline, conservative, and Karatsuba models. |
| `verify_lemma1.py` | Verifies Lemma 1's closed-form expression for $`B`$ against Figure 5 Step 4 across $`40`$ random instances. |
| `structural.py` | Standalone transcription of Figure 5 and Figure 10 verifying vanishing second differences and recovering $`20 / 20`$ keys. |
| `recover_key.py` | Implements Algorithm 1 ($`257 \times 256`$ $`\mathbb{F}_2`$ Gaussian elimination and $`\mathsf{NSBC}`$ disambiguation). |
| `rank_sweep.py` | Measures the rank-deficiency distribution of $`M_\sigma`$ across $`N`$ fresh keys (matching uniform $`257 \times 256`$ $`\mathbb{F}_2`$ matrices: $`\Pr[d = 0, 1, 2] \approx 0.578, 0.385, 0.037`$; also prints a $`\mathrm{Poisson}(1/2)`$ reference). |
| `sbcvole.py` | Pure-Python $`\mathbb{F}_{2^{257}}`$ arithmetic ($`\theta^{257} + \theta^{12} + 1`$), standalone protocol transcription, and $`\mathbb{F}_2`$ linear solver. |
| `forge.c` | Zero-patch C driver linking against 100% unmodified `SBCvole.c`; extracts $`(\Delta^*, B, \vec{F}_x(\Delta^*), \vec{F}_y(\Delta^*))`$ directly from raw binary signature and public-key files via `CheckVOLESignature`. |
| `vendor/` *(fetched on demand)* | Populated automatically via `../fetch_third_party.sh` with `asiacrypt-2025-a10.zip` and `crypto-algorithms/{aes,sha256}.{c,h}`. |
