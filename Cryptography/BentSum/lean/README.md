# Lean formalization

Lean 4 + Mathlib formalization of the main results of `../paper/main.tex`.
The development follows the mathematics of the paper but not its structure or
numbering.

`lake build` succeeds with no warnings and no `sorry`; no `native_decide` is
used, and the main theorems depend only on the standard axioms:

```
'BentSum.IsBentSum.two_pow_dvd_wt'      depends on axioms: [propext, Classical.choice, Quot.sound]
'BentSum.not_bentSumConjecture'         depends on axioms: [propext, Classical.choice, Quot.sound]
'BentSum.not_derivativeConjecture'      depends on axioms: [propext, Classical.choice, Quot.sound]
```

## Build

Toolchain `leanprover/lean4:v4.35.0-rc3`, Mathlib tag `v4.35.0-rc3`.

```
cd lean
lake exe cache get     # first time only
lake build
```

To print the axioms:

```
TMP=$(mktemp -d)
cat > "$TMP/ax.lean" <<'EOF'
import BentSum
#print axioms BentSum.IsBentSum.two_pow_dvd_wt
#print axioms BentSum.not_bentSumConjecture
#print axioms BentSum.not_derivativeConjecture
EOF
lake env lean "$TMP/ax.lean"
rm -rf "$TMP"
```

## Mathlib style

The code follows the Mathlib conventions, although it is not intended for
submission:

- `lakefile.toml` enables the Mathlib linter set (`linter.mathlibStandardSet`,
  including the line-length, flexible-`simp` and `show` linters; only the
  copyright-header linter is disabled) and
  `linter.missingDocs`; the build produces no warnings.
- `#lint in BentSum` (Batteries' environment linters: `docBlame`, `simpNF`,
  `unusedArguments`, …) passes.
- Every file has a module docstring with *Main
  definitions*, *Main statements* and *References*; imports are minimal
  (no `import Mathlib`); there is no `set_option` in the sources.
- Names follow the Mathlib naming convention; Mathlib's own notions are used
  where they exist (`⬝ᵥ`, `Pi.single`, `Function.Periodic`, `CharTwo.*`), and no
  name clashes with a root-level Mathlib declaration.

## Main statements

`BentSumConjecture` and `DerivativeConjecture` are the bent-sum decomposition
conjecture and the derivative conjecture of `../paper/main.tex`; the Lean names
describe the statements, and the attributions are in the `References` sections of
the source files.

```lean
def BentSumConjecture (n : ℕ) : Prop :=
  ∀ f : BoolFn n, InRM (n / 2) f → IsBentSum f

def DerivativeConjecture (n : ℕ) : Prop :=
  ∀ F : BoolFn n, IsBalanced F → InRM (n / 2 - 1) F →
    (∃ y : Vec n, y ≠ 0 ∧ Function.Periodic F y) →
      ∃ (y : Vec n) (g : BoolFn n), IsBent g ∧ F = dirDeriv y g

theorem IsBentSum.two_pow_dvd_wt {f : BoolFn (2 * m)} (hf : IsBentSum f) (hk : k < m)
    {S : Finset (Fin (2 * m))} (hS : #S = k) {c : Vec (2 * m)}
    (hsupp : ∀ x, f x = 1 → ∀ i ∈ S, x i = c i) : 2 ^ (k + 1) ∣ wt f
theorem not_isBentSum_cylQ3 (m : ℕ) (hm : 4 ≤ m) : ¬ IsBentSum (cylQ3 m)
theorem not_bentSumConjecture (m : ℕ) (hm : 4 ≤ m) : ¬ BentSumConjecture (2 * m)
theorem liftCylQ3_ne_dirDeriv (q : ℕ) (hq : 4 ≤ q) (y) (hg : IsBent g) :
    liftCylQ3 q ≠ dirDeriv y g
theorem not_derivativeConjecture (q : ℕ) (hq : 4 ≤ q) : ¬ DerivativeConjecture (2 + 2 * q)
```

`not_bentSumConjecture_eight` (`n = 8`) and `not_derivativeConjecture_ten`
(`n = 10`) are the smallest instances; `exists_isBent_eight` and
`exists_isBent_ten` show that `IsBent` is not vacuous in those dimensions.

## Files

```
BentSum.lean                root module
BentSum/Basic.lean          F₂ⁿ, sgn, walshHadamard, walsh, IsBent, wt, IsBentSum, InRM,
                            dirDeriv, IsBalanced, the two conjectures, q3, cylQ3; the
                            splitting x = (x₀, x₁, z) (ins, tail2, e0, liftFn);
                            IsBent.eq_zero_of_periodic; IsBent.hyperbolicAdd; isBent_q3
BentSum/Congruence.lean     Fourier analysis on a coordinate subspace (zeroOn); the
                            weight congruence IsBentSum.two_pow_dvd_wt
BentSum/Counterexample.lean weight, degree and support of cylQ3 m; not_isBentSum_cylQ3;
                            the constant derivatives of cylQ3 m (cylQ3_linearStructure)
BentSum/Derivative.lean     the bent-sum reduction (sliceCofactor,
                            IsBentSum.of_dirDeriv_e0); IsBent.comp_transvection;
                            liftCylQ3 q = x₁ + cylQ3 q and liftCylQ3_ne_dirDeriv
```

## Definitions to review

These are the only places where the formal statement can differ from the
paper; everything else is checked by Lean.

- `Vec n := Fin n → ZMod 2`, `BoolFn n := Vec n → ZMod 2`. Indices are
  0-based: the paper's `x₁, …, x_n` are `x 0, …, x (n-1)`.
- `sgn a := if a = 0 then 1 else -1 : ℤ`,
  `walsh f ω := ∑ x, sgn (f x) * sgn (x ⬝ᵥ ω)`.
- `IsBent f := ∀ ω, walsh f ω ^ 2 = 2 ^ n`.
- `wt f := #{x | f x = 1}`; `IsBalanced f := 2 * wt f = 2 ^ n`.
- `IsBentSum f := ∃ g h, IsBent g ∧ IsBent h ∧ f = g + h`.
- `InRM r f`: `f` is the evaluation of an `MvPolynomial (Fin n) (ZMod 2)` of
  total degree `≤ r`, i.e. `deg f ≤ r` for the algebraic normal form.
- `dirDeriv y g x := g x + g (x + y)`.
- `cylQ3 m`: the product of the coordinates `0, …, m-3` times `q3` of the
  coordinates `m-2, …, m+3` (the paper's `f_{2m}`); `liftCylQ3 q := liftFn (cylQ3 q)`,
  i.e. `x ↦ x 1 + cylQ3 q (x 2, …)` (the paper's `x₂ + f_{n-2}`).

## Differences from the paper

- **The weight congruence is proved for coordinate cosets**
  `{x | ∀ i ∈ S, x i = c i}` with `#S = k < m`, which is all the counterexample
  needs.  The paper states it for an arbitrary affine subspace of codimension
  `k`; the reduction by an affine change of variables is not formalized.
- **The derivative counterexample uses a single transvection.**  If
  `liftCylQ3 q = D_y g` with `g` bent, `y` is a period of `liftCylQ3 q`, so
  `w = e₀ + y` does not affect it (`cylQ3_linearStructure`).  Since `y ≠ 0`,
  some linear form `ℓ` has `ℓ(e₀) = 1` and `ℓ(w) = 0`; the transvection
  `T x = x + ℓ(x) w` preserves bentness, carries `e₀` to `y` and fixes
  `liftCylQ3 q`, so `D_{e₀} (g ∘ T) = liftCylQ3 q`, and the bent-sum reduction
  makes `cylQ3 q` a sum of two bent functions.
- Only the direction of the reduction lemma that is needed is proved, and it
  needs no parity hypothesis on the dimension.
- Exact degrees and the bentness of the general quadratics `q_ℓ` are not
  formalized; only `q3` is (derived from `x₀x₁` by `IsBent.hyperbolicAdd`), and
  its weight `28` is derived the same way (`wt_hyperbolicAdd`).
