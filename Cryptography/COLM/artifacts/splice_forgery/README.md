# Two-Query Encryption-Only Splice Forgery on $`\mathsf{COLM}_0`$

This directory contains the self-contained C verification harness, unmodified official CAESAR SUPERCOP implementations (`crypto_aead/colm0v1/ref`, `colm0v1/aesni4`, `colm0v1/aesni8`, and `colm127v1/ref`), and direct 2024 *Journal of Cryptology* specification implementations ($`\mathsf{AES\text{-}COLM}_0\text{-}128`$ and $`\mathsf{AES\text{-}COLM}_0\text{-}256`$) demonstrating all four variants of the two-query encryption-only splice forgery on $`\mathsf{COLM}_0`$:

1. **Difference Fade** (`fade`, $`\ell \ge 130`$): Splices blocks $`S_{\mathrm{fade}} = \{1, 2, 3, 8, 129\} = \mathrm{supp}(f) + 1`$ from a second encryption query (under a distinct nonce $`N_2 \neq N_1`$ or distinct associated data $`A_2 \neq A_1`$) into the first ciphertext, retaining the first query's final block and tag block verbatim.
2. **State Transfer and Tag Transplant** (`transplant`, $`\ell \ge 129`$): Splices blocks $`S_{\mathrm{trans}}(\ell, m) = (\mathrm{supp}(3^{-1}x^{m-1}) + 1) \cup \{m, \dots, \ell+1\}`$ ($`129 \le m \le \ell`$, e.g. $`\{1, 3, 4, 5, 6, 7\} \cup \{129, \dots, \ell+1\}`$ at $`m = 129`$) to transfer the decryption chaining state onto the second query's trajectory and transplant its tag block.
3. **Common-Prefix Relaxation** (`prefix`, $`\ell_1, \ell_2 \ge 130`$): Demonstrates that the two queried plaintexts need only agree on their first $`129`$ blocks ($`2{,}064`$ bytes) while differing arbitrarily in subsequent blocks, total message length, nonce, and associated-data length.
4. **Nonce-Misuse Variant** (`misuse`, $`\ell \ge s + 130`$): Splices the shifted support $`S_{\mathrm{fade}} + s`$ when two messages differing at block $`s`$ are encrypted under a repeated nonce.

## Key Files

- `splice.c`: Consolidated test driver deriving every splice index set at runtime over $`\mathbb{F}_{2^{128}} = \mathbb{F}_2[x]/(x^{128} + x^7 + x^2 + x + 1)`$ and running positive trials alongside drop-one, add-one, and random-subset negative controls.
- `joc2024_colm0.c`: Line-by-line C implementation of the 2024 *Journal of Cryptology* specification (Figures 2 and 3), supporting both $`\mathsf{AES\text{-}COLM}_0\text{-}128`$ (`-DCOLM_KEY_BITS=128`) and $`\mathsf{AES\text{-}COLM}_0\text{-}256`$ (`-DCOLM_KEY_BITS=256`).
- `xcheck.c`: Cross-implementation consistency harness verifying byte-identical ciphertexts and joint forgery acceptance across independent $`\mathsf{COLM}_0`$ implementations.
- `onequery.c`: Single-query baseline check confirming $`0`$ acceptances across $`8`$ single-ciphertext manipulation strategies.
- `aes_kat.c` & `openssl/aes.h`: Self-contained FIPS-197 `AES-128`/`AES-256` header shim and Appendix C.1 Known-Answer Test (no external `libssl-dev` / `-lcrypto` installation required).
- `witnesses/` & `run_demo.log`: Archived trial-0 binary inputs/forged ciphertexts ($`C^{\mathrm{mix}}`$) and reference `--demo` execution log.

## Usage

Requires a standard C compiler (`gcc` or `clang`, honoring `CC`). If the SUPERCOP `colm0v1_*` / `colm127v1_ref` directories are not yet present locally, `./run.sh` automatically invokes `../fetch_third_party.sh --ensure` and notifies the user on `stderr` before compiling. Portable targets run on 64-bit Linux and macOS (`x86_64` and `arm64`); the optional `colm0v1/aesni4` and `colm0v1/aesni8` cross-checks run automatically on `x86_64` CPUs supporting `AES-NI` and `AVX2` (and print a skip notice otherwise):

```bash
./run.sh --demo   # Fast verification (25 trials per variant + 16 boundary sweeps, ~10 s)
./run.sh --full   # Full verification (1,000 trials per variant + 2,000 random controls, ~2 min)
```
