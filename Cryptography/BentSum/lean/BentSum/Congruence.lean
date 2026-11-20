import BentSum.Basic
import Mathlib.Tactic.LinearCombination

/-!
# The weight congruence for sums of two bent functions

If `f = g + h` with `g, h` bent on `F₂^{2m}` and the support of `f` lies in a
coordinate coset `{x | ∀ i ∈ S, x i = c i}` with `#S = k < m`, then `2^{k+1}`
divides `wt f`.

The signed indicator `φ = 1_{f = 1} · (-1)^g` has Fourier values in `{0, ±2^m}`.
Parseval and Fourier inversion on the coordinate subspace `{ω | ∀ i ∈ S, ω i = 0}`
then give the congruence modulo `2^{m+1}`.

## Main definitions

* `BentSum.signedInd`: the signed indicator `1_{f = 1} · (-1)^g`.
* `BentSum.zeroOn`: the coordinate subspace `{ω | ∀ i ∈ S, ω i = 0}`.

## Main statements

* `BentSum.walshHadamard_signedInd_eq`: the Fourier values of the signed indicator.
* `BentSum.IsBentSum.two_pow_dvd_wt`: the weight congruence.
-/

open Finset

namespace BentSum

/-- The signed indicator `1_{f = 1} · (-1)^g`. -/
def signedInd {n : ℕ} (f g : BoolFn n) (x : Vec n) : ℤ := if f x = 1 then sgn (g x) else 0

private theorem sgn_sub_sgn (a b : ZMod 2) :
    sgn a - sgn b = 2 * (if a + b = 1 then sgn a else 0) := by
  fin_cases a <;> fin_cases b <;> decide

theorem two_mul_walshHadamard_signedInd {n : ℕ} {g h : BoolFn n} (ω : Vec n) :
    2 * walshHadamard (signedInd (g + h) g) ω = walsh g ω - walsh h ω := by
  simp only [signedInd, walsh, walshHadamard, mul_sum, ← sum_sub_distrib, ← sub_mul, sgn_sub_sgn,
    Pi.add_apply, mul_assoc]

/-- The Walsh coefficients of a bent function on `F₂^{2m}` are `± 2^m`. -/
theorem IsBent.walsh_eq_or {m : ℕ} {g : BoolFn (2 * m)} (hg : IsBent g) (ω : Vec (2 * m)) :
    walsh g ω = 2 ^ m ∨ walsh g ω = -2 ^ m := by
  have := hg ω
  rw [show (2 : ℤ) ^ (2 * m) = (2 ^ m) ^ 2 by rw [← pow_mul, mul_comm]] at this
  exact sq_eq_sq_iff_eq_or_eq_neg.mp this

/-- The signed indicator of a sum of two bent functions has Fourier values in
`{0, ±2^m}`. -/
theorem walshHadamard_signedInd_eq {m : ℕ} {g h : BoolFn (2 * m)} (hg : IsBent g) (hh : IsBent h)
    (ω : Vec (2 * m)) :
    walshHadamard (signedInd (g + h) g) ω = 0 ∨ walshHadamard (signedInd (g + h) g) ω = 2 ^ m ∨
      walshHadamard (signedInd (g + h) g) ω = -2 ^ m := by
  have e := two_mul_walshHadamard_signedInd (g := g) (h := h) ω
  generalize walshHadamard (signedInd (g + h) g) ω = F at e ⊢
  rcases hg.walsh_eq_or ω with h1 | h1 <;> rcases hh.walsh_eq_or ω with h2 | h2 <;>
    rw [h1, h2] at e <;> generalize (2 : ℤ) ^ m = t at e ⊢ <;> omega

/-- The coordinate subspace `{ω | ∀ i ∈ S, ω i = 0}`. -/
def zeroOn {n : ℕ} (S : Finset (Fin n)) : Finset (Vec n) := {ω | ∀ i ∈ S, ω i = 0}

theorem card_zeroOn {n : ℕ} (S : Finset (Fin n)) : #(zeroOn S) = 2 ^ (n - #S) := by
  have : zeroOn S = Fintype.piFinset
      (fun i ↦ if i ∈ S then ({0} : Finset (ZMod 2)) else univ) := by
    ext ω
    simp only [zeroOn, mem_filter, mem_univ, true_and, Fintype.mem_piFinset]
    refine ⟨fun h i ↦ ?_, fun h i hi ↦ by simpa [hi] using h i⟩
    split_ifs with hi
    · simp [h i hi]
    · simp
  rw [this, Fintype.card_piFinset]
  have : ∀ i, #(if i ∈ S then ({0} : Finset (ZMod 2)) else univ)
      = 2 ^ (if i ∈ S then 0 else 1) := by
    intro i; split_ifs <;> simp
  simp only [this, prod_pow_eq_pow_sum]
  congr 1
  simp [Finset.sum_ite, Finset.filter_not, card_univ_sdiff]

theorem sgn_add_one (a : ZMod 2) : sgn (a + 1) = - sgn a := by
  fin_cases a <;> decide

/-- Orthogonality on `zeroOn S`: the sum of the characters `x ⬝ ·` over `zeroOn S` is
`#(zeroOn S)` if `x` vanishes outside `S`, and `0` otherwise. -/
theorem sum_zeroOn_sgn_dotProduct {n : ℕ} (S : Finset (Fin n)) (x : Vec n) :
    ∑ ω ∈ zeroOn S, sgn (x ⬝ᵥ ω) = if ∀ i ∉ S, x i = 0 then (#(zeroOn S) : ℤ) else 0 := by
  split_ifs with hx
  · rw [Finset.card_eq_sum_ones, Nat.cast_sum]
    refine sum_congr rfl fun ω hω ↦ ?_
    simp only [zeroOn, mem_filter, mem_univ, true_and] at hω
    have : x ⬝ᵥ ω = 0 := sum_eq_zero fun i _ ↦ by
      by_cases hi : i ∈ S
      · simp [hω i hi]
      · simp [hx i hi]
    simp [this, sgn]
  · push Not at hx
    obtain ⟨j, hjS, hj⟩ := hx
    refine sum_involution (fun ω _ ↦ ω + Pi.single j 1) (fun ω _ ↦ ?_) (fun ω _ _ h ↦ ?_)
      (fun ω hω ↦ ?_) (fun ω _ ↦ ?_)
    · rw [dotProduct_add, dotProduct_single_one, zmod2_eq_one_of_ne_zero hj, sgn_add_one]
      ring
    · simpa using congrFun h j
    · simp only [zeroOn, mem_filter, mem_univ, true_and] at hω ⊢
      intro i hi
      simp [hω i hi, Pi.single_eq_of_ne (fun e : i = j ↦ hjS (e ▸ hi))]
    · rw [add_assoc, vec_add_self, add_zero]

/-- Fourier inversion on `zeroOn S` for a function supported on the coset
`{x | ∀ i ∈ S, x i = c i}`. -/
theorem walshHadamard_inversion_zeroOn {n : ℕ} (S : Finset (Fin n)) (c : Vec n) (φ : Vec n → ℤ)
    (hφ : ∀ x, φ x ≠ 0 → ∀ i ∈ S, x i = c i) (y : Vec n) (hy : ∀ i ∈ S, y i = c i) :
    ∑ ω ∈ zeroOn S, walshHadamard φ ω * sgn (y ⬝ᵥ ω) = (#(zeroOn S) : ℤ) * φ y := by
  have : ∀ ω, walshHadamard φ ω * sgn (y ⬝ᵥ ω) = ∑ x, φ x * sgn ((x + y) ⬝ᵥ ω) := fun ω ↦ by
    simp only [walshHadamard, sum_mul, add_dotProduct, sgn_add, mul_assoc]
  simp only [this]
  rw [sum_comm]
  simp only [← mul_sum, sum_zeroOn_sgn_dotProduct]
  rw [sum_eq_single y]
  · rw [ite_eq_left (fun i _ ↦ by rw [vec_add_self]; rfl)]
    ring
  · intro x _ hxy
    by_cases h0 : φ x = 0
    · simp [h0]
    suffices hn : ¬ ∀ i ∉ S, (x + y) i = 0 by rw [ite_eq_right hn, mul_zero]
    refine fun hc ↦ hxy (funext fun i ↦ CharTwo.add_eq_zero.mp ?_)
    by_cases hi : i ∈ S
    · rw [hφ x h0 i hi, hy i hi]; exact CharTwo.add_self_eq_zero _
    · exact hc i hi
  · simp

/-- Parseval's identity on `zeroOn S`. -/
theorem parseval_zeroOn {n : ℕ} (S : Finset (Fin n)) (c : Vec n) (φ : Vec n → ℤ)
    (hφ : ∀ x, φ x ≠ 0 → ∀ i ∈ S, x i = c i) :
    ∑ ω ∈ zeroOn S, walshHadamard φ ω ^ 2 = (#(zeroOn S) : ℤ) * ∑ x, φ x ^ 2 := by
  have : ∀ ω, walshHadamard φ ω ^ 2 = ∑ x, φ x * (walshHadamard φ ω * sgn (x ⬝ᵥ ω)) := by
    intro ω
    rw [sq]
    nth_rewrite 2 [walshHadamard]
    rw [mul_sum]
    exact sum_congr rfl fun x _ ↦ by ring
  simp only [this]
  rw [sum_comm, mul_sum]
  refine sum_congr rfl fun x _ ↦ ?_
  rw [← mul_sum]
  by_cases h0 : φ x = 0
  · simp [h0]
  rw [walshHadamard_inversion_zeroOn S c φ hφ x (hφ x h0)]
  ring

/-- **The weight congruence.** If `f` is a sum of two bent functions on `F₂^{2m}` and its
support lies in the coset `{x | ∀ i ∈ S, x i = c i}` with `#S = k < m`, then `2^{k+1}`
divides `wt f`. -/
theorem IsBentSum.two_pow_dvd_wt {m k : ℕ} {f : BoolFn (2 * m)} (hf : IsBentSum f)
    (hk : k < m) {S : Finset (Fin (2 * m))} (hS : #S = k) {c : Vec (2 * m)}
    (hsupp : ∀ x, f x = 1 → ∀ i ∈ S, x i = c i) : 2 ^ (k + 1) ∣ wt f := by
  obtain ⟨g, h, hg, hh, rfl⟩ := hf
  set φ := signedInd (g + h) g
  have hφ : ∀ x, φ x ≠ 0 → ∀ i ∈ S, x i = c i := fun x hx ↦
    hsupp x (by_contra fun h1 ↦ hx (by simp only [φ, signedInd, h1, ↓reduceIte]))
  by_cases hw : wt (g + h) = 0
  · simp [hw]
  obtain ⟨y, hy⟩ : ∃ y, (g + h) y = 1 := by
    by_contra hne
    push Not at hne
    exact hw (card_eq_zero.mpr (filter_eq_empty_iff.mpr fun x _ ↦ hne x))
  have hcard : (#(zeroOn S) : ℤ) = 2 ^ (2 * m - k) := by
    rw [card_zeroOn, hS]; push_cast; rfl
  have hvals := walshHadamard_signedInd_eq hg hh
  set N := #{ω ∈ zeroOn S | walshHadamard φ ω ≠ 0}
  -- Parseval: `wt f = N · 2^k`.
  have hsq : ∀ ω, walshHadamard φ ω ^ 2 =
      if walshHadamard φ ω ≠ 0 then (2 : ℤ) ^ (2 * m) else 0 := by
    intro ω
    rcases hvals ω with h0 | h0 | h0 <;> rw [h0] <;> simp <;> ring
  have hφsq : ∀ x, φ x ^ 2 = if (g + h) x = 1 then 1 else 0 := by
    intro x
    simp only [φ, signedInd]
    split_ifs <;> simp
  have hpar := parseval_zeroOn S c φ hφ
  simp only [hsq, hφsq, ← sum_filter, sum_const, nsmul_eq_mul, hcard, mul_one] at hpar
  change (N : ℤ) * _ = _ * (wt (g + h) : ℤ) at hpar
  have hwt : (wt (g + h) : ℤ) = N * 2 ^ k := by
    rw [show (2 : ℤ) ^ (2 * m) = 2 ^ (2 * m - k) * 2 ^ k by rw [← pow_add]; congr 1; omega]
      at hpar
    apply mul_left_cancel₀ (pow_ne_zero (2 * m - k) (two_ne_zero : (2 : ℤ) ≠ 0))
    linear_combination -hpar
  -- Inversion at `y`: `N` is even.
  have hinv := walshHadamard_inversion_zeroOn S c φ hφ y (hsupp y hy)
  rw [hcard] at hinv
  have hterm : ∀ ω, (2 : ℤ) ^ (m + 1) ∣
      walshHadamard φ ω * sgn (y ⬝ᵥ ω) - if walshHadamard φ ω ≠ 0 then 2 ^ m else 0 := by
    intro ω
    rcases hvals ω with h0 | h0 | h0 <;>
      rcases sgn_eq_one_or_neg_one (y ⬝ᵥ ω) with h1 | h1 <;> rw [h0, h1]
    all_goals split_ifs with hF
    all_goals first
      | exact absurd rfl hF
      | exact (hF (pow_ne_zero _ two_ne_zero)).elim
      | exact (hF (neg_ne_zero.mpr (pow_ne_zero _ two_ne_zero))).elim
      | (refine ⟨0, ?_⟩; ring1)
      | (refine ⟨-1, ?_⟩; ring1)
  have hsum := dvd_sum (s := zeroOn S) fun ω _ ↦ hterm ω
  rw [sum_sub_distrib, hinv, ← sum_filter, sum_const, nsmul_eq_mul] at hsum
  have hpow : (2 : ℤ) ^ (m + 1) ∣ 2 ^ (2 * m - k) * φ y :=
    dvd_mul_of_dvd_left (pow_dvd_pow 2 (by omega)) _
  have hNdiv := dvd_sub hpow hsum
  rw [sub_sub_cancel, pow_succ, mul_comm (N : ℤ)] at hNdiv
  obtain ⟨t, ht⟩ : 2 ∣ N := by
    exact_mod_cast (mul_dvd_mul_iff_left (pow_ne_zero m (two_ne_zero : (2 : ℤ) ≠ 0))).mp hNdiv
  rw [show wt (g + h) = N * 2 ^ k by exact_mod_cast hwt, ht]
  exact ⟨t, by ring⟩

end BentSum
