import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.MvPolynomial.Degrees
import Mathlib.Algebra.Ring.Periodic
import Mathlib.Data.Matrix.Mul

/-!
# Boolean functions, the Walsh transform and bent functions

This file sets up Boolean functions on `F₂ⁿ = Fin n → ZMod 2` and the facts about bent
functions used by the rest of the development.  Indices are 0-based.

## Main definitions

* `BentSum.walshHadamard`: the Walsh–Hadamard transform of an integer-valued function.
* `BentSum.walsh`: the Walsh transform `W_f(ω) = ∑ₓ (-1)^{f(x) + x ⬝ ω}`.
* `BentSum.IsBent`: `f` is bent if `W_f(ω)² = 2ⁿ` for every `ω`.
* `BentSum.IsBentSum`: `f` is a sum of two bent functions.
* `BentSum.InRM`: `f` has algebraic degree at most `r`.
* `BentSum.dirDeriv`: the derivative `D_y g (x) = g(x) + g(x + y)`.
* `BentSum.BentSumConjecture`, `BentSum.DerivativeConjecture`: the bent-sum conjecture
  and the derivative conjecture.
* `BentSum.ins`, `BentSum.tail2`: the splitting `x = (x₀, x₁, z)` of `F₂^{2+k}`.
* `BentSum.q3`, `BentSum.cylQ3`: the functions of the counterexamples.

## Main statements

* `BentSum.IsBent.eq_zero_of_periodic`: a bent function has no nonzero period.
* `BentSum.IsBent.hyperbolicAdd`: if `b` is bent, so is `x₀x₁ + b(z)`.
* `BentSum.isBent_q3`: `q3` is bent.

## References

* [N. Tokareva, *On decomposition of a Boolean function into sum of bent
  functions*][tokareva2014]
* [N. Tokareva, *On the set of derivatives of a Boolean bent function*][tokareva2016]

## Tags

bent function, Walsh transform, Boolean function
-/

open Finset

namespace BentSum

/-! ## Definitions -/

/-- Vectors in `F₂ⁿ`. -/
abbrev Vec (n : ℕ) := Fin n → ZMod 2

/-- Boolean functions on `F₂ⁿ`. -/
abbrev BoolFn (n : ℕ) := Vec n → ZMod 2

/-- The sign character `(-1)^a`, with values in `ℤ`. -/
def sgn (a : ZMod 2) : ℤ := if a = 0 then 1 else -1

/-- The Walsh–Hadamard transform `φ̂(ω) = ∑ₓ φ(x) (-1)^{x ⬝ ω}` of `φ : F₂ⁿ → ℤ`. -/
def walshHadamard {n : ℕ} (φ : Vec n → ℤ) (ω : Vec n) : ℤ :=
  ∑ x, φ x * sgn (x ⬝ᵥ ω)

/-- The Walsh transform `W_f(ω) = ∑ₓ (-1)^{f(x) + x ⬝ ω}`. -/
def walsh {n : ℕ} (f : BoolFn n) (ω : Vec n) : ℤ :=
  walshHadamard (fun x ↦ sgn (f x)) ω

/-- `f` is bent if `W_f(ω)² = 2ⁿ` for every `ω`. -/
def IsBent {n : ℕ} (f : BoolFn n) : Prop := ∀ ω, walsh f ω ^ 2 = 2 ^ n

/-- The Hamming weight of `f`. -/
def wt {n : ℕ} (f : BoolFn n) : ℕ := #{x | f x = 1}

/-- `f` is a sum of two bent functions. -/
def IsBentSum {n : ℕ} (f : BoolFn n) : Prop :=
  ∃ g h : BoolFn n, IsBent g ∧ IsBent h ∧ f = g + h

/-- `f ∈ RM(r, n)`: `f` is the evaluation of a polynomial of total degree `≤ r`. -/
def InRM {n : ℕ} (r : ℕ) (f : BoolFn n) : Prop :=
  ∃ p : MvPolynomial (Fin n) (ZMod 2),
    p.totalDegree ≤ r ∧ ∀ x, f x = MvPolynomial.eval x p

/-- The bent-sum conjecture in dimension `n`: every Boolean
function on `F₂ⁿ` of algebraic degree at most `n/2` is a sum of two bent functions. -/
def BentSumConjecture (n : ℕ) : Prop :=
  ∀ f : BoolFn n, InRM (n / 2) f → IsBentSum f

/-- The derivative `D_y g (x) = g(x) + g(x + y)` in the direction `y`. -/
def dirDeriv {n : ℕ} (y : Vec n) (g : BoolFn n) : BoolFn n := fun x ↦ g x + g (x + y)

/-- `f` is balanced: `wt f = 2^{n-1}`, written without subtraction. -/
def IsBalanced {n : ℕ} (f : BoolFn n) : Prop := 2 * wt f = 2 ^ n

/-- The derivative conjecture in dimension `n`: a balanced function of degree at
most `n/2 - 1` with a nonzero period is a derivative of a bent function. -/
def DerivativeConjecture (n : ℕ) : Prop :=
  ∀ F : BoolFn n, IsBalanced F → InRM (n / 2 - 1) F →
    (∃ y : Vec n, y ≠ 0 ∧ Function.Periodic F y) →
      ∃ (y : Vec n) (g : BoolFn n), IsBent g ∧ F = dirDeriv y g

/-- `q3 z = z₀z₁ + z₂z₃ + z₄z₅` on `F₂⁶`. -/
def q3 : BoolFn 6 := fun z ↦ z 0 * z 1 + z 2 * z 3 + z 4 * z 5

/-- Coordinate `j` of `x`, or `0` if `j` is out of range. -/
def coord {n : ℕ} (x : Vec n) (j : ℕ) : ZMod 2 := if h : j < n then x ⟨j, h⟩ else 0

/-- `cylQ3 m` on `F₂^{2m}`: the product of the coordinates `0, …, m-3` times `q3` of the
coordinates `m-2, …, m+3`. -/
def cylQ3 (m : ℕ) : BoolFn (2 * m) := fun x ↦
  (∏ j ∈ range (m - 2), coord x j) *
    (coord x (m - 2) * coord x (m - 1) + coord x m * coord x (m + 1)
      + coord x (m + 2) * coord x (m + 3))

/-! ### Elementary facts about `ZMod 2` and `sgn` -/

theorem sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  fin_cases a <;> fin_cases b <;> decide

theorem sgn_eq_one_or_neg_one (a : ZMod 2) : sgn a = 1 ∨ sgn a = -1 := by
  unfold sgn; split_ifs <;> simp

@[simp]
theorem sgn_sq (a : ZMod 2) : sgn a ^ 2 = 1 := by
  rcases sgn_eq_one_or_neg_one a with h | h <;> rw [h] <;> norm_num

theorem zmod2_eq_one_of_ne_zero {a : ZMod 2} (h : a ≠ 0) : a = 1 := by
  revert a; decide

theorem zmod2_eq_zero_or_eq_one (a : ZMod 2) : a = 0 ∨ a = 1 := by revert a; decide

theorem zmod2_mul_eq_one (a b : ZMod 2) : a * b = 1 ↔ a = 1 ∧ b = 1 := by
  revert a b; decide

theorem zmod2_eq_add_of_add_eq {x y s : ZMod 2} (h : x + y = s) : y = x + s := by
  revert x y s; decide

theorem zmod2_add_eq_of_add_eq {u v w : ZMod 2} (h : u + v = w) : w + v = u := by
  revert u v w; decide

theorem prod_range_eq_one_iff (k : ℕ) (g : ℕ → ZMod 2) :
    ∏ j ∈ range k, g j = 1 ↔ ∀ j < k, g j = 1 := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [prod_range_succ, zmod2_mul_eq_one, ih]
    constructor
    · rintro ⟨h1, h2⟩ j hj
      rcases Nat.lt_succ_iff_lt_or_eq.1 hj with h | rfl
      exacts [h1 j h, h2]
    · exact fun h ↦ ⟨fun j hj ↦ h j (by omega), h k (by omega)⟩

theorem coord_of_lt {n : ℕ} (x : Vec n) (j : ℕ) (h : j < n) : coord x j = x ⟨j, h⟩ := by
  simp [coord, h]

theorem vec_add_self {n : ℕ} (v : Vec n) : v + v = 0 :=
  funext fun _ ↦ CharTwo.add_self_eq_zero _

/-!
## Splitting off the first two coordinates

`x = (x₀, x₁, z) ∈ F₂ × F₂ × F₂^k`, with the index forms `Fin.castAdd k 0`,
`Fin.castAdd k 1` and `Fin.natAdd 2` produced by `Fin.sum_univ_add`.
-/

section

variable {k : ℕ}

/-- The coordinate `0` of `F₂^{2+k}`. -/
abbrev j0 (k : ℕ) : Fin (2 + k) := Fin.castAdd k 0

/-- The coordinate `1` of `F₂^{2+k}`. -/
abbrev j1 (k : ℕ) : Fin (2 + k) := Fin.castAdd k 1

/-- Drop the first two coordinates. -/
def tail2 (x : Vec (2 + k)) : Vec k := fun i ↦ x (Fin.natAdd 2 i)

/-- Insert two coordinates in front: `ins a b z = (a, b, z)`. -/
def ins (a b : ZMod 2) (z : Vec k) : Vec (2 + k) := fun i ↦
  if h0 : i.val = 0 then a
  else if h1 : i.val = 1 then b
  else z ⟨i.val - 2, by have := i.isLt; omega⟩

/-- The standard basis vector `e₀ = (1, 0, 0)`. -/
def e0 (k : ℕ) : Vec (2 + k) := ins 1 0 0

@[simp]
theorem ins_j0 (a b : ZMod 2) (z : Vec k) : ins a b z (j0 k) = a := by
  simp [ins, j0]

@[simp]
theorem ins_j1 (a b : ZMod 2) (z : Vec k) : ins a b z (j1 k) = b := by
  simp [ins, j1]

@[simp]
theorem ins_natAdd (a b : ZMod 2) (z : Vec k) (i : Fin k) :
    ins a b z (Fin.natAdd 2 i) = z i := by
  have h : (Fin.natAdd 2 i).val = 2 + i.val := rfl
  simp only [ins, h]
  split_ifs with h0 h1
  · exact absurd h0 (by omega)
  · exact absurd h1 (by omega)
  · exact congrArg z (Fin.ext (by simp))

@[simp]
theorem tail2_ins (a b : ZMod 2) (z : Vec k) : tail2 (ins a b z) = z := by
  funext i; simp [tail2]

/-- Every vector is an insertion: `x = (x₀, x₁, z)`. -/
theorem ins_tail2 (x : Vec (2 + k)) : ins (x (j0 k)) (x (j1 k)) (tail2 x) = x := by
  funext ⟨i, hi⟩
  match i with
  | 0 => rfl
  | 1 => rfl
  | i + 2 => simp [ins, tail2, Nat.add_comm]

@[simp]
theorem ins_add (a b a' b' : ZMod 2) (z z' : Vec k) :
    ins a b z + ins a' b' z' = ins (a + a') (b + b') (z + z') := by
  funext i
  simp only [ins, Pi.add_apply]
  split_ifs <;> rfl

@[simp]
theorem tail2_add (x y : Vec (2 + k)) : tail2 (x + y) = tail2 x + tail2 y := by
  funext i; simp [tail2]

/-- The function `x₁ + h₀(z)`. -/
def liftFn (h0 : BoolFn k) : BoolFn (2 + k) := fun x ↦ x (j1 k) + h0 (tail2 x)

@[simp]
theorem liftFn_ins (h0 : BoolFn k) (a b : ZMod 2) (z : Vec k) :
    liftFn h0 (ins a b z) = b + h0 z := by
  simp [liftFn]

/-- `ins` is a bijection onto `F₂^{2+k}`. -/
def insEquiv (k : ℕ) : ZMod 2 × ZMod 2 × Vec k ≃ Vec (2 + k) where
  toFun p := ins p.1 p.2.1 p.2.2
  invFun x := (x (j0 k), x (j1 k), tail2 x)
  left_inv p := by simp
  right_inv x := by simp [ins_tail2]

/-- Summation split over the first two coordinates. -/
theorem sum_ins {M : Type*} [AddCommMonoid M] (φ : Vec (2 + k) → M) :
    ∑ x, φ x = ∑ a : ZMod 2, ∑ b : ZMod 2, ∑ z : Vec k, φ (ins a b z) := by
  rw [← (insEquiv k).sum_comp φ, Fintype.sum_prod_type]
  simp only [Fintype.sum_prod_type]
  rfl

theorem ins_dotProduct_ins (a b : ZMod 2) (z : Vec k) (c d : ZMod 2) (ω : Vec k) :
    ins a b z ⬝ᵥ ins c d ω = a * c + b * d + z ⬝ᵥ ω := by
  simp only [dotProduct, Fin.sum_univ_add, Fin.sum_univ_two, ins_natAdd]
  rw [show Fin.castAdd k (0 : Fin 2) = j0 k from rfl,
    show Fin.castAdd k (1 : Fin 2) = j1 k from rfl, ins_j0, ins_j0, ins_j1, ins_j1]

/-- The function `x₀x₁ + b(z)` on `F₂^{2+k}`. -/
def hyperbolicAdd (b : BoolFn k) : BoolFn (2 + k) := fun x ↦ x (j0 k) * x (j1 k) + b (tail2 x)

end

/-! ### Reed–Muller membership -/

theorem inRM_comp_tail2 {k r : ℕ} {h : BoolFn k} (hh : InRM r h) :
    InRM r (fun x : Vec (2 + k) ↦ h (tail2 x)) := by
  obtain ⟨p, hpdeg, hp⟩ := hh
  exact ⟨MvPolynomial.rename (Fin.natAdd 2) p,
    (MvPolynomial.totalDegree_rename_le _ _).trans hpdeg,
    fun x ↦ by rw [MvPolynomial.eval_rename, ← hp]; rfl⟩

theorem InRM.add {n r : ℕ} {f g : BoolFn n} (hf : InRM r f) (hg : InRM r g) :
    InRM r (f + g) := by
  obtain ⟨p, hp, hpf⟩ := hf
  obtain ⟨q, hq, hqg⟩ := hg
  exact ⟨p + q, (MvPolynomial.totalDegree_add p q).trans (max_le hp hq),
    fun x ↦ by rw [Pi.add_apply, MvPolynomial.eval_add, hpf, hqg]⟩

/-! ## Facts about bent functions -/

/-- A bent function has no nonzero period. -/
theorem IsBent.eq_zero_of_periodic {n : ℕ} {g : BoolFn n} (hg : IsBent g) {a : Vec n}
    (h : Function.Periodic g a) : a = 0 := by
  by_contra hne
  obtain ⟨j, hj⟩ := Function.ne_iff.mp hne
  have hdot : a ⬝ᵥ Pi.single j 1 = 1 := by
    rw [dotProduct_single_one]; exact zmod2_eq_one_of_ne_zero hj
  have key : walsh g (Pi.single j 1) = - walsh g (Pi.single j 1) := by
    calc walsh g (Pi.single j 1)
        = ∑ x, sgn (g (x + a)) * sgn ((x + a) ⬝ᵥ Pi.single j 1) :=
          (Fintype.sum_equiv (Equiv.addRight a) _ _ fun _ ↦ rfl).symm
      _ = - walsh g (Pi.single j 1) := by
          simp only [show ∀ x, g (x + a) = g x from h, add_dotProduct, hdot, sgn_add, walsh,
            walshHadamard, ← sum_neg_distrib]
          exact sum_congr rfl fun x _ ↦ by simp [sgn]
  have h2 : (2 : ℤ) ^ n = 0 := by rw [← hg (Pi.single j 1), show walsh g _ = 0 by omega]; ring
  exact pow_ne_zero n two_ne_zero h2

/-! ### Adding a hyperbolic pair preserves bentness

If `b` is bent on `F₂^k`, then `x₀x₁ + b(z)` is bent on `F₂^{2+k}`: its Walsh
transform is `W(c, d, β) = 2 (-1)^{cd} W_b(β)`.
-/

private theorem sgn_split (a b c d u v : ZMod 2) :
    sgn (a * b + v) * sgn (a * c + b * d + u)
      = sgn (a * b + a * c + b * d) * (sgn v * sgn u) := by
  revert a b c d u v; decide

private theorem sum_sq_char (c d : ZMod 2) :
    ∑ a : ZMod 2, ∑ b : ZMod 2, sgn (a * b + a * c + b * d) = 2 * sgn (c * d) := by
  revert c d; decide

theorem walsh_hyperbolicAdd {k : ℕ} (b : BoolFn k) (c d : ZMod 2) (β : Vec k) :
    walsh (hyperbolicAdd b) (ins c d β) = 2 * sgn (c * d) * walsh b β := by
  simp only [walsh, walshHadamard, hyperbolicAdd, sum_ins, ins_j0, ins_j1, tail2_ins,
    ins_dotProduct_ins, sgn_split, ← mul_sum, ← sum_mul, sum_sq_char]

/-- If `b` is bent on `F₂^k`, then `x₀x₁ + b(z)` is bent on `F₂^{2+k}`. -/
theorem IsBent.hyperbolicAdd {k : ℕ} {b : BoolFn k} (hb : IsBent b) :
    IsBent (hyperbolicAdd b) := by
  intro ω
  rw [← ins_tail2 ω, walsh_hyperbolicAdd, mul_pow, mul_pow, sgn_sq, hb, pow_add]
  norm_num

theorem isBent_x0x1 : IsBent (fun x : Vec 2 ↦ x 0 * x 1) := by
  unfold IsBent walsh walshHadamard dotProduct sgn
  decide

theorem q3_eq : q3 = hyperbolicAdd (hyperbolicAdd fun x : Vec 2 ↦ x 0 * x 1) := by
  funext x; simp only [q3, hyperbolicAdd, add_assoc]; rfl

theorem isBent_q3 : IsBent q3 := q3_eq ▸ isBent_x0x1.hyperbolicAdd.hyperbolicAdd

/-- Bent functions exist on `F₂⁸`, so `not_bentSumConjecture_eight` is not vacuous. -/
theorem exists_isBent_eight : ∃ g : BoolFn 8, IsBent g :=
  ⟨_, isBent_q3.hyperbolicAdd⟩

/-- Bent functions exist on `F₂¹⁰`, so `not_derivativeConjecture_ten` is not vacuous. -/
theorem exists_isBent_ten : ∃ g : BoolFn 10, IsBent g :=
  ⟨_, isBent_q3.hyperbolicAdd.hyperbolicAdd⟩

end BentSum
