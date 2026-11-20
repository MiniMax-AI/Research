import BentSum.Counterexample

/-!
# A counterexample to the derivative conjecture

For `q ≥ 4` the function `liftCylQ3 q = x₁ + cylQ3 q (x₂, …)` on `F₂^{2+2q}` is
balanced, has degree at most `q = n/2 - 1` and the nonzero period `e₀`, but it is
not the derivative of a bent function in any direction.

The proof has three steps.

* **The bent-sum reduction** (`BentSum.IsBentSum.of_dirDeriv_e0`): if `x₁ + h₀(z)` is
  the derivative of a bent function in the direction `e₀`, then `h₀` is a sum of two
  bent functions.
* **Transvections** `x ↦ x + (x ⬝ ℓ) w` with `w ⬝ ℓ = 0` are involutions and preserve
  bentness (`BentSum.IsBent.comp_transvection`).
* If `liftCylQ3 q = D_y g`, then `y` is a period of `liftCylQ3 q`, which forces
  `w = e₀ + y` to be invisible to it (`BentSum.cylQ3_linearStructure`).  A transvection
  along `w` that carries `e₀` to `y` turns `g` into a bent function whose derivative in
  the direction `e₀` is `liftCylQ3 q`, so `cylQ3 q` is a sum of two bent functions,
  contradicting `BentSum.not_isBentSum_cylQ3`.

## Main definitions

* `BentSum.transvection`: the transvection `x ↦ x + (x ⬝ ℓ) w`.
* `BentSum.liftCylQ3`: the counterexample.

## Main statements

* `BentSum.liftCylQ3_ne_dirDeriv`: `liftCylQ3 q` is not the derivative of a bent function.
* `BentSum.not_derivativeConjecture`: the derivative conjecture fails in every
  even dimension `2 + 2q ≥ 10`.

## References

* [N. Tokareva, *On the set of derivatives of a Boolean bent function*][tokareva2016]
* [A. S. Shaporenko, *On decomposition of bent functions in 8 variables into the sum of
  two bent functions*][shaporenko2022]
-/

open Finset

namespace BentSum

/-! ## The bent-sum reduction -/

private theorem sgn_aux (G a b t c d u : ZMod 2) :
    sgn (G + a * (b + t)) * sgn (a * c + b * d + u)
      = sgn (a * (b + (t + c))) * (sgn G * sgn (b * d) * sgn u) := by
  revert G a b t c d u; decide

private theorem sgn_aux2 (G t c d u : ZMod 2) :
    (2 : ℤ) * (sgn G * sgn ((t + c) * d) * sgn u)
      = sgn (c * d) * (2 * (sgn (G + d * t) * sgn u)) := by
  revert G t c d u; decide

/-- Only the term `b = t` survives the sum over the first coordinate. -/
private theorem sum_sgn_collapse (t : ZMod 2) (Ψ : ZMod 2 → ℤ) :
    ∑ b : ZMod 2, (∑ a : ZMod 2, sgn (a * (b + t))) * Ψ b = 2 * Ψ t := by
  have h (F : ZMod 2 → ℤ) : ∑ b : ZMod 2, F b = F 0 + F 1 := Fin.sum_univ_two F
  rcases zmod2_eq_zero_or_eq_one t with rfl | rfl <;> simp only [h] <;> simp +decide [sgn]

/-- `k_c(z) = g(0, h₀(z) + c, z)`, the cofactor of `g` along the slice `x₀ = 0`,
`x₁ = h₀(z) + c`. -/
def sliceCofactor {k : ℕ} (g : BoolFn (2 + k)) (h0 : BoolFn k) (c : ZMod 2) : BoolFn k :=
  fun z ↦ g (ins 0 (h0 z + c) z)

/-- If `D_{e₀} g = x₁ + h₀(z)` then `g(a, b, z) = g(0, b, z) + a (b + h₀(z))`. -/
theorem apply_ins_of_dirDeriv_e0 {k : ℕ} {g : BoolFn (2 + k)} {h0 : BoolFn k}
    (hd : dirDeriv (e0 k) g = liftFn h0) (a b : ZMod 2) (z : Vec k) :
    g (ins a b z) = g (ins 0 b z) + a * (b + h0 z) := by
  rcases zmod2_eq_zero_or_eq_one a with rfl | rfl
  · simp
  · have h := congrFun hd (ins 0 b z)
    simp only [dirDeriv, e0, ins_add, zero_add, add_zero, liftFn_ins] at h
    rw [one_mul]
    exact zmod2_eq_add_of_add_eq h

/-- If `D_{e₀} g = x₁ + h₀(z)` then `W_g(c, d, β) = (-1)^{cd} · 2 · W_{k_c + d h₀}(β)`. -/
theorem walsh_ins_of_dirDeriv_e0 {k : ℕ} {g : BoolFn (2 + k)} {h0 : BoolFn k}
    (hd : dirDeriv (e0 k) g = liftFn h0) (c d : ZMod 2) (β : Vec k) :
    walsh g (ins c d β) =
      sgn (c * d) * (2 * walsh (fun z ↦ sliceCofactor g h0 c z + d * h0 z) β) := by
  have key : ∀ (b : ZMod 2) (z : Vec k),
      ∑ a : ZMod 2, sgn (g (ins a b z)) * sgn (ins a b z ⬝ᵥ ins c d β)
        = (∑ a : ZMod 2, sgn (a * (b + (h0 z + c)))) *
            (sgn (g (ins 0 b z)) * sgn (b * d) * sgn (z ⬝ᵥ β)) := by
    intro b z
    rw [sum_mul]
    exact sum_congr rfl fun a _ ↦ by
      rw [apply_ins_of_dirDeriv_e0 hd, ins_dotProduct_ins]; exact sgn_aux ..
  rw [walsh, walshHadamard, sum_ins, sum_comm,
    sum_congr rfl fun b _ ↦ (sum_comm.trans (sum_congr rfl fun z _ ↦ key b z)), sum_comm]
  simp only [walsh, walshHadamard, mul_sum, sliceCofactor]
  refine sum_congr rfl fun z _ ↦ ?_
  rw [sum_sgn_collapse]
  exact sgn_aux2 ..

/-- **The bent-sum reduction.** If `x₁ + h₀(z)` is the derivative of a bent function in
the direction `e₀`, then `h₀` is a sum of two bent functions (namely `k₀` and
`k₀ + h₀`). -/
theorem IsBentSum.of_dirDeriv_e0 {k : ℕ} {h0 : BoolFn k} {g : BoolFn (2 + k)} (hg : IsBent g)
    (hd : dirDeriv (e0 k) g = liftFn h0) : IsBentSum h0 := by
  have main : ∀ d : ZMod 2, IsBent (fun z ↦ sliceCofactor g h0 0 z + d * h0 z) := by
    intro d β
    have h := hg (ins 0 d β)
    rw [walsh_ins_of_dirDeriv_e0 hd, zero_mul, show sgn 0 = 1 from rfl, one_mul, pow_add] at h
    refine mul_left_cancel₀ (a := (4 : ℤ)) (by norm_num) ?_
    linear_combination h
  refine ⟨_, _, main 0, main 1, funext fun z ↦ ?_⟩
  simp only [Pi.add_apply, zero_mul, add_zero, one_mul]
  exact (CharTwo.add_cancel_left _ _).symm

/-! ## Transvections -/

section

variable {n : ℕ}

/-- The transvection `x ↦ x + (x ⬝ ℓ) w`. -/
def transvection (ℓ w x : Vec n) : Vec n := x + (x ⬝ᵥ ℓ) • w

theorem transvection_add (ℓ w x y : Vec n) :
    transvection ℓ w (x + y) = transvection ℓ w x + transvection ℓ w y := by
  simp only [transvection, add_dotProduct, add_smul]; abel

theorem transvection_involutive {ℓ w : Vec n} (h : w ⬝ᵥ ℓ = 0) :
    Function.Involutive (transvection ℓ w) := by
  intro x
  simp only [transvection, add_dotProduct, smul_dotProduct, h, smul_eq_mul, mul_zero,
    add_zero, add_assoc, ← add_smul, CharTwo.add_self_eq_zero, zero_smul]

/-- A transvection with `w ⬝ ℓ = 0` preserves bentness: it transports the Walsh
transform to the dual transvection. -/
theorem IsBent.comp_transvection {ℓ w : Vec n} {g : BoolFn n} (hg : IsBent g)
    (h : w ⬝ᵥ ℓ = 0) : IsBent (fun x ↦ g (transvection ℓ w x)) := by
  intro ω
  have hT := transvection_involutive h
  have : walsh (fun x ↦ g (transvection ℓ w x)) ω = walsh g (ω + (w ⬝ᵥ ω) • ℓ) := by
    rw [walsh, walshHadamard, ← Equiv.sum_comp (hT.toPerm _)]
    refine sum_congr rfl fun u _ ↦ ?_
    simp only [Function.Involutive.coe_toPerm, hT u]
    simp only [transvection, add_dotProduct, smul_dotProduct, dotProduct_add, dotProduct_smul,
      smul_eq_mul, mul_comm (u ⬝ᵥ ℓ)]
  rw [this, hg]

/-- A linear form separating the basis vector `e_i` from any other vector `w`. -/
theorem exists_separating_form (i : Fin n) {w : Vec n} (hw : w ≠ Pi.single i 1) :
    ∃ ℓ, Pi.single i 1 ⬝ᵥ ℓ = 1 ∧ w ⬝ᵥ ℓ = 0 := by
  rcases zmod2_eq_zero_or_eq_one (w i) with h | h
  · exact ⟨Pi.single i 1, by simp, by rw [dotProduct_single_one, h]⟩
  · obtain ⟨j, hj⟩ := Function.ne_iff.mp hw
    have hji : j ≠ i := by rintro rfl; simp [h] at hj
    rw [Pi.single_eq_of_ne hji] at hj
    refine ⟨Pi.single i 1 + Pi.single j 1, ?_, ?_⟩ <;>
      simp only [dotProduct_add, dotProduct_single_one, Pi.single_eq_same,
        Pi.single_eq_of_ne hji, h, zmod2_eq_one_of_ne_zero hj]
    all_goals decide

end

/-! ## The counterexample -/

/-- The counterexample `x₁ + cylQ3 q (x₂, …)` on `F₂^{2+2q}`. -/
def liftCylQ3 (q : ℕ) : BoolFn (2 + 2 * q) := liftFn (cylQ3 q)

/-- `liftCylQ3 q` is invariant under vectors that vanish on `x₁` and on the coordinates
of `cylQ3 q`. -/
theorem liftCylQ3_add_of_vanish {q : ℕ} (hq : 4 ≤ q) {v : Vec (2 + 2 * q)}
    (hv1 : v (j1 _) = 0) (hv : ∀ j : Fin (2 * q), j.val < q + 4 → tail2 v j = 0)
    (x : Vec (2 + 2 * q)) : liftCylQ3 q (x + v) = liftCylQ3 q x := by
  simp only [liftCylQ3, liftFn, Pi.add_apply, hv1, add_zero, tail2_add,
    cylQ3_add_of_vanish q hq _ _ hv]

theorem e0_eq_single (k : ℕ) : e0 k = Pi.single (j0 k) 1 := by
  rw [← ins_tail2 (Pi.single (j0 k) 1), e0, Pi.single_eq_same,
    Pi.single_eq_of_ne (by simp [Fin.ext_iff])]
  congr 1; funext i
  exact (Pi.single_eq_of_ne (M := fun _ ↦ ZMod 2)
    (show Fin.natAdd 2 i ≠ j0 k by simp [Fin.ext_iff]) 1).symm

theorem e0_ne_zero (k : ℕ) : e0 k ≠ 0 := fun h ↦ by simpa [e0] using congrFun h (j0 k)

theorem isBalanced_liftCylQ3 (q : ℕ) : IsBalanced (liftCylQ3 q) := by
  have h (c : ZMod 2) : ∑ b : ZMod 2, (if b + c = 1 then 1 else 0 : ℕ) = 1 := by
    revert c; decide
  rw [IsBalanced, wt, card_filter, sum_ins]
  simp only [liftCylQ3, liftFn_ins]
  rw [sum_comm (s := univ (α := ZMod 2)) (t := univ (α := Vec (2 * q)))]
  simp only [h, sum_const, card_univ, Fintype.card_fun, ZMod.card, Fintype.card_fin,
    smul_eq_mul, mul_one, pow_add]
  ring

theorem inRM_liftCylQ3 (q : ℕ) (hq : 4 ≤ q) : InRM q (liftCylQ3 q) :=
  InRM.add (f := fun x ↦ x (j1 (2 * q))) ⟨MvPolynomial.X (j1 (2 * q)), by simp; omega, by simp⟩
    (inRM_comp_tail2 (inRM_cylQ3 q hq))

theorem periodic_liftCylQ3_e0 (q : ℕ) (hq : 4 ≤ q) : Function.Periodic (liftCylQ3 q) (e0 _) :=
  liftCylQ3_add_of_vanish hq (by simp [e0]) fun _ _ ↦ by simp [e0]

theorem liftCylQ3_ne_zero (q : ℕ) (hq : 4 ≤ q) : liftCylQ3 q ≠ 0 := fun h ↦ by
  simpa [liftCylQ3, cylQ3_one q hq] using congrFun h (ins 0 0 1)

/-- A derivative is invariant under translation by its own direction. -/
theorem periodic_dirDeriv {n : ℕ} (y : Vec n) (g : BoolFn n) :
    Function.Periodic (dirDeriv y g) y := fun x ↦ by
  simp only [dirDeriv, add_assoc, vec_add_self, add_zero, add_comm (g (x + y))]

/-- `liftCylQ3 q` is not the derivative of a bent function, in any direction. -/
theorem liftCylQ3_ne_dirDeriv (q : ℕ) (hq : 4 ≤ q) (y : Vec (2 + 2 * q))
    {g : BoolFn (2 + 2 * q)} (hg : IsBent g) : liftCylQ3 q ≠ dirDeriv y g := by
  intro h
  -- `y` is a period of `liftCylQ3 q`, so `D_{tail2 y} (cylQ3 q)` is the constant `y₁` …
  have key : ∀ z, dirDeriv (tail2 y) (cylQ3 q) z = y (j1 _) := fun z ↦ by
    have h1 := periodic_dirDeriv y g (ins 0 0 z)
    rw [← h, ← ins_tail2 y, ins_add] at h1
    simp only [liftCylQ3, liftFn_ins, zero_add] at h1
    exact zmod2_add_eq_of_add_eq h1
  -- … hence `w = e₀ + y` is invisible to `liftCylQ3 q`.
  obtain ⟨hy1, ha⟩ := cylQ3_linearStructure q hq key
  set w := e0 (2 * q) + y with hw
  have hfree : ∀ (c : ZMod 2) x, liftCylQ3 q (x + c • w) = liftCylQ3 q x := fun c x ↦
    liftCylQ3_add_of_vanish hq (by simp [hw, e0, hy1])
      (fun j hj ↦ by have := ha j hj; simp only [tail2] at this ⊢; simp [hw, e0, this]) x
  -- `y ≠ 0`, so `w ≠ e₀` and some transvection along `w` carries `e₀` to `y`.
  have hwe : w ≠ Pi.single (j0 _) 1 := by
    rw [← e0_eq_single, hw]
    intro he
    apply liftCylQ3_ne_zero q hq
    rw [h, show y = 0 by simpa using he]
    funext x
    simp [dirDeriv, CharTwo.add_self_eq_zero]
  obtain ⟨ℓ, hℓ1, hℓw⟩ := exists_separating_form _ hwe
  have hTe : transvection ℓ w (e0 (2 * q)) = y := by
    rw [transvection, e0_eq_single, hℓ1, one_smul, ← e0_eq_single, hw, ← add_assoc,
      vec_add_self, zero_add]
  -- `g ∘ T` is bent and `D_{e₀} (g ∘ T) = liftCylQ3 q`.
  have hd : dirDeriv (e0 (2 * q)) (fun x ↦ g (transvection ℓ w x)) = liftFn (cylQ3 q) := by
    funext x
    rw [dirDeriv, transvection_add, hTe, ← liftCylQ3, ← hfree (x ⬝ᵥ ℓ) x]
    exact (congrFun h _).symm
  exact not_isBentSum_cylQ3 q hq (.of_dirDeriv_e0 (hg.comp_transvection hℓw) hd)

/-- The derivative conjecture fails in every even dimension `2 + 2q ≥ 10`. -/
theorem not_derivativeConjecture (q : ℕ) (hq : 4 ≤ q) : ¬ DerivativeConjecture (2 + 2 * q) := by
  intro hconj
  obtain ⟨y, g, hg, hFy⟩ := hconj (liftCylQ3 q) (isBalanced_liftCylQ3 q)
    (by rw [show (2 + 2 * q) / 2 - 1 = q by omega]; exact inRM_liftCylQ3 q hq)
    ⟨e0 _, e0_ne_zero _, periodic_liftCylQ3_e0 q hq⟩
  exact liftCylQ3_ne_dirDeriv q hq y hg hFy

/-- The derivative conjecture fails in dimension `10`. -/
theorem not_derivativeConjecture_ten : ¬ DerivativeConjecture 10 :=
  not_derivativeConjecture 4 le_rfl

end BentSum
