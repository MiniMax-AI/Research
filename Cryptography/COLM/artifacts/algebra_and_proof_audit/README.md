# $`\mathbb{F}_{2^{128}}`$ Algebraic Derivations and Security Proof Audit

This directory contains pure Python 3 scripts (`Python 3.8+`, standard library only, zero external dependencies) verifying all algebraic claims and security-proof evaluations in the paper:

- `algebra.py`: Derives the fade splice set $`S_{\mathrm{fade}} = \{1, 2, 3, 8, 129\}`$ and transplant splice set $`S_{\mathrm{trans}} = \{1, 3, 4, 5, 6, 7, 129, 130\}`$ directly from $`\mathbb{F}_{2^{128}} = \mathbb{F}_2[x]/(x^{128} + x^7 + x^2 + x + 1)`$, proves the exact length floors $`\ell_{\min} = 130`$ (`fade`) and $`\ell_{\min} = 129`$ (`transplant`) via $`\mathbb{F}_2`$-subspace rank computation, and proves the $`\mathsf{COLM}_\tau`$ obstruction $`\deg(3^{-1}x^\tau \bmod f(x)) = 127 > \tau - 1`$ for all $`1 \le \tau \le 127`$.
- `derivation.py`: Verifies Case A (pure $`\Delta\mathsf{IV}`$ difference), Case B (pure block-$`1`$ message difference under nonce reuse, $`\vert S_B \vert = 123`$, $`\max(S_B) = 129`$), and Case C (simultaneous $`\Delta\mathsf{IV}`$ and block-$`1`$ message difference forces $`S = \emptyset`$, confirmed both algebraically and by exhaustive subset enumeration over $`\mathbb{F}_{2^m}`$ for $`m \in \{4, \dots, 12\}`$).
- `thm2.py`: Evaluates Theorem 2 of both the 2024 *Journal of Cryptology* full proof ($`\mathsf{Adv}^{\mathrm{oAE}} \le 2^{-105.92}`$ at $`\sigma \le 795`$) and the 2016 CAESAR Round-3 security proof addendum ($`\mathsf{Adv}^{\mathrm{auth}} \le 2^{-107.56}`$ at $`\sigma_{\mathrm{auth}} \le 398`$) at the exact parameters of our two-query adversary ($`q = 2`$, $`s = 1`$, $`\ell \in \{129, 130\}`$), comparing both claimed bounds against the measured advantage $`1`$.
- `run_demo.log`: Archived reference output from `./run.sh`.

## Usage

```bash
./run.sh
```
