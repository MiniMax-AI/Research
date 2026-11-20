# $`\mathsf{HVZK}`$ & Game 2 Simulator Distinguisher (`hvzk_distinguisher`)

This directory contains the pure-Python experimental verification for Lemma 1 / Theorem 1 (Section 4.2–4.3) and Theorem 2 (Section 5.2) of the paper ([`../../paper/main.pdf`](../../paper/main.pdf)):

1. **Closed-Form Affine Identity Across Padding Conventions** (`[A]`, over $`\mathbb{F}_2`$):
   Verifies that $`\mathrm{vec}(B) = M_\sigma w \oplus \mathrm{vec}(\gamma_\sigma)`$ holds with $`0`$ mismatches across both the homogeneous $`\Delta`$-scaled padding convention ($`\eta = 0`$, used in `SBCvole.c`) and the literal printed padding convention ($`\eta = 1`$, Figure 5 Step 4), and recovers $`100\%`$ of planted $`\mathsf{NSBC}`$ keys from a single transcript.
2. **Distinguisher Against Witness-Free Simulators in Game 2** ($`\mathcal{D}`$, `[B]`):
   - **Real Signing Oracle** ($`\mathcal{O}(\mathsf{sk}, \mathsf{pk})`$): Accepted with probability $`1`$ ($`120 / 120`$ in `--full` mode).
   - **Simulator Strategy 1** ($`\mathcal{S}_1`$: uniform $`B \in \mathbb{F}_{2^{257}}`$, $`A \leftarrow \Phi(\Delta^*) \oplus B\Delta^*`$): Rejected $`120 / 120`$ times ($`\ge 1/2`$ fail the linear span test $`\mathrm{vec}(B \oplus \gamma_\sigma) \in \mathrm{im}(M_\sigma)`$ alone).
   - **Simulator Strategy 2** ($`\mathcal{S}_2`$: random dummy witness $`(\tilde{\vec{x}}, \tilde{\vec{y}})`$ as in Appendix A, p. 35): Rejected $`120 / 120`$ times via the public $`\mathsf{NSBC}`$ relation check.
   - **Repaired Protocol** ($`\tilde{B} = B \oplus \rho`$): Confirms that masking $`B`$ with an independent uniform one-time pad $`\rho \in \mathbb{F}_{2^{257}}`$ makes the transcript compatible with every candidate key $`w \in \mathbb{F}_2^{256}`$ ($`120 / 120`$).

## Requirements

- **Python**: `Python 3.8+` (standard library only; zero third-party packages or network access required) and `bash`.
- **Platforms**: 64-bit Linux and macOS (`x86_64` and `arm64`).

## Usage

```bash
./run.sh --demo   # 30 instances per padding convention & simulator strategy (~6 s)
./run.sh --full   # 120 instances per padding convention & simulator strategy (~25 s)
```
