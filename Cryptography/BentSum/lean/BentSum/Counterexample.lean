import BentSum.Congruence

/-!
# A counterexample to the bent-sum conjecture

For `m ≥ 4`, `cylQ3 m` has degree at most `m` and weight `7 · 2^{m-2}`, and its support
lies in the coset `{x₀ = … = x_{m-3} = 1}` of codimension `m - 2`.  The weight
congruence would force `2^{m-1} ∣ wt (cylQ3 m)`, so `cylQ3 m` is not a sum of two bent
functions.

The file also determines the constant derivatives of `cylQ3 m`, which is used for the
derivative counterexample.

## Main statements

* `BentSum.wt_cylQ3`: `wt (cylQ3 m) = 7 · 2^{m-2}`.
* `BentSum.not_isBentSum_cylQ3`: `cylQ3 m` is not a sum of two bent functions.
* `BentSum.not_bentSumConjecture`: the bent-sum conjecture fails in every even
  dimension `2m ≥ 8`.
* `BentSum.cylQ3_linearStructure`: a constant derivative `D_a (cylQ3 m)` is `0`, and `a`
  vanishes on the coordinates `0, …, m+3`.

## References

* [N. Tokareva, *On decomposition of a Boolean function into sum of bent
  functions*][tokareva2014]
-/

open Finset

namespace BentSum

/-! ## The weight of `cylQ3 m` -/

theorem wt_hyperbolicAdd {k : ℕ} (b : BoolFn k) : wt (hyperbolicAdd b) = 2 ^ k + 2 * wt b := by
  have h (t : ZMod 2) : ∑ a : ZMod 2, ∑ c : ZMod 2, (if a * c + t = 1 then 1 else 0 : ℕ)
      = 1 + 2 * if t = 1 then 1 else 0 := by revert t; decide
  have hs (F : ZMod 2 → ZMod 2 → Vec k → ℕ) :
      ∑ a, ∑ c, ∑ z, F a c z = ∑ z, ∑ a, ∑ c, F a c z :=
    (sum_congr rfl fun _ _ ↦ sum_comm).trans sum_comm
  rw [wt, wt, card_filter, card_filter, sum_ins]
  simp only [hyperbolicAdd, ins_j0, ins_j1, tail2_ins]
  rw [hs]
  simp only [h, sum_add_distrib, ← mul_sum, sum_const, card_univ, Fintype.card_fun, ZMod.card,
    Fintype.card_fin, smul_eq_mul, mul_one]

theorem wt_q3 : wt q3 = 28 := by
  rw [q3_eq, wt_hyperbolicAdd, wt_hyperbolicAdd, show wt (fun x : Vec 2 ↦ x 0 * x 1) = 1 by decide]
  rfl

/-- The support of `cylQ3 (t+4)`: the coordinates `0, …, t+1` are `1`, and the `q3`
block is `1`. -/
theorem cylQ3_eq_one_iff (t : ℕ) (x : Vec (2 * (t + 4))) :
    cylQ3 (t + 4) x = 1 ↔ (∀ j < t + 2, coord x j = 1) ∧
      q3 (fun i : Fin 6 ↦ coord x (t + 2 + i)) = 1 := by
  unfold cylQ3
  rw [zmod2_mul_eq_one, prod_range_eq_one_iff]
  have h1 : t + 4 - 2 = t + 2 := by omega
  have h2 : t + 4 - 1 = t + 3 := by omega
  rw [h1, h2]
  simp only [q3]
  rfl

/-- The equivalence between the support of `cylQ3 (t+4)` and `supp q3 × Vec t`. -/
def cylQ3SuppEquiv (t : ℕ) :
    {x : Vec (2 * (t + 4)) // cylQ3 (t + 4) x = 1} ≃ {b : Vec 6 // q3 b = 1} × Vec t where
  toFun x := (⟨fun i ↦ coord x.1 (t + 2 + i), ((cylQ3_eq_one_iff t x.1).1 x.2).2⟩,
    fun i ↦ coord x.1 (t + 8 + i))
  invFun p := ⟨fun j ↦ if j.1 < t + 2 then 1 else if h : j.1 < t + 8 then
      p.1.1 ⟨j.1 - (t + 2), by omega⟩ else p.2 ⟨j.1 - (t + 8), by omega⟩, by
    rw [cylQ3_eq_one_iff]
    refine ⟨fun j hj ↦ ?_, ?_⟩
    · rw [coord_of_lt _ _ (by omega)]; simp [hj]
    · convert p.1.2 using 2
      funext i
      rw [coord_of_lt _ _ (by omega)]
      simp only [show ¬ (t + 2 + (i : ℕ) < t + 2) by omega, ↓reduceIte,
        show t + 2 + (i : ℕ) < t + 8 by omega, ↓reduceDIte]
      exact congrArg p.1.1 (Fin.ext (by simp))⟩
  left_inv x := Subtype.ext <| funext fun j ↦ by
    have := ((cylQ3_eq_one_iff t x.1).1 x.2).1 j.1
    dsimp only
    split_ifs with h1 h2
    · simpa [coord_of_lt _ _ j.2] using (this h1).symm
    all_goals rw [coord_of_lt _ _ (by omega)]; exact congrArg x.1 (Fin.ext (by simp; omega))
  right_inv p := by
    obtain ⟨⟨b, hb⟩, c⟩ := p
    refine Prod.ext (Subtype.ext (funext fun i ↦ ?_)) (funext fun i ↦ ?_) <;> dsimp only <;>
      rw [coord_of_lt _ _ (by omega)]
    · simp only [show ¬ (t + 2 + (i : ℕ) < t + 2) by omega, ↓reduceIte,
        show t + 2 + (i : ℕ) < t + 8 by omega, ↓reduceDIte]
      exact congrArg b (Fin.ext (by simp))
    · simp only [show ¬ (t + 8 + (i : ℕ) < t + 2) by omega, ↓reduceIte,
        show ¬ (t + 8 + (i : ℕ) < t + 8) by omega, ↓reduceDIte]
      exact congrArg c (Fin.ext (by simp))

/-- `wt (cylQ3 m) = 28 · 2^{m-4} = 7 · 2^{m-2}` for `m ≥ 4`. -/
theorem wt_cylQ3 (m : ℕ) (hm : 4 ≤ m) : wt (cylQ3 m) = 7 * 2 ^ (m - 2) := by
  obtain ⟨t, rfl⟩ : ∃ t, m = t + 4 := ⟨m - 4, by omega⟩
  have hQ : Fintype.card {b : Vec 6 // q3 b = 1} = 28 := by
    rw [Fintype.card_subtype]; exact wt_q3
  have : wt (cylQ3 (t + 4)) = Fintype.card {x : Vec (2 * (t + 4)) // cylQ3 (t + 4) x = 1} := by
    rw [Fintype.card_subtype]; rfl
  rw [this, Fintype.card_congr (cylQ3SuppEquiv t), Fintype.card_prod, hQ]
  simp only [Fintype.card_fun, ZMod.card, Fintype.card_fin]
  rw [show t + 4 - 2 = t + 2 by omega, pow_add]
  ring

/-! ### Weights, Reed–Muller membership and supports -/

/-- The coordinates `0, …, m-3` of the product factor of `cylQ3 m`. -/
def prodCoords (m : ℕ) : Finset (Fin (2 * m)) := {i | i.val < m - 2}

/-- The product factor of `cylQ3 m` has `m - 2` coordinates. -/
theorem card_prodCoords (m : ℕ) (hm : 4 ≤ m) : #(prodCoords m) = m - 2 := by
  rw [prodCoords, Fin.card_filter_val_lt]; omega

/-- The coordinate polynomial of the coordinate `j` of `F₂^{2m}`. -/
noncomputable def coordPoly (m j : ℕ) : MvPolynomial (Fin (2 * m)) (ZMod 2) :=
  if h : j < 2 * m then MvPolynomial.X ⟨j, h⟩ else 0

/-- The coordinate polynomial has degree at most `1`. -/
theorem totalDegree_coordPoly_le (m j : ℕ) : (coordPoly m j).totalDegree ≤ 1 := by
  unfold coordPoly; split_ifs <;> simp [MvPolynomial.totalDegree_X]

/-- The coordinate polynomial evaluates to the coordinate. -/
theorem eval_coordPoly (m j : ℕ) (x : Vec (2 * m)) :
    MvPolynomial.eval x (coordPoly m j) = coord x j := by
  unfold coordPoly coord; split_ifs <;> simp

/-- `cylQ3 m` has degree at most `m = n/2`. -/
theorem inRM_cylQ3 (m : ℕ) (hm : 4 ≤ m) : InRM m (cylQ3 m) := by
  set C := coordPoly m
  have h1 : ∀ j, (C j).totalDegree ≤ 1 := totalDegree_coordPoly_le m
  have h2 (i j : ℕ) : (C i * C j).totalDegree ≤ 2 :=
    (MvPolynomial.totalDegree_mul _ _).trans (by have := h1 i; have := h1 j; omega)
  refine ⟨(∏ j ∈ range (m - 2), C j) *
    (C (m - 2) * C (m - 1) + C m * C (m + 1) + C (m + 2) * C (m + 3)), ?_, ?_⟩
  · have hp : (∏ j ∈ range (m - 2), C j).totalDegree ≤ m - 2 :=
      (MvPolynomial.totalDegree_finsetProd _ _).trans
        ((sum_le_sum fun j _ ↦ h1 j).trans (by simp))
    have h3 := (MvPolynomial.totalDegree_add _ _).trans
      (max_le ((MvPolynomial.totalDegree_add _ _).trans
        (max_le (h2 (m - 2) (m - 1)) (h2 m (m + 1)))) (h2 (m + 2) (m + 3)))
    have := MvPolynomial.totalDegree_mul (∏ j ∈ range (m - 2), C j)
      (C (m - 2) * C (m - 1) + C m * C (m + 1) + C (m + 2) * C (m + 3))
    omega
  · intro x
    simp [C, cylQ3, map_prod, eval_coordPoly]

/-! ## The structure of `cylQ3 m` and its constant derivatives

`cylQ3 m` depends only on the coordinates `0, …, m+3`; a constant derivative
`D_a (cylQ3 m)` is `0`, and `a` vanishes on these coordinates.
-/

/-- The product part of `cylQ3 m`: the coordinates `0, …, m-3`. -/
def cylQ3Prod (m : ℕ) (x : Vec (2 * m)) : ZMod 2 := ∏ j ∈ range (m - 2), coord x j

/-- The `q3` block of `cylQ3 m`: the six coordinates `m-2, …, m+3`, renumbered as `Fin 6`. -/
def cylQ3Block (m : ℕ) (x : Vec (2 * m)) : Vec 6 := fun i ↦ coord x (m - 2 + i.val)

/-- `cylQ3 m` is the product of the product part and of `q3` on the block. -/
theorem cylQ3_eq_mul (m : ℕ) (hm : 4 ≤ m) (x : Vec (2 * m)) :
    cylQ3 m x = cylQ3Prod m x * q3 (cylQ3Block m x) := by
  obtain ⟨t, rfl⟩ : ∃ t, m = t + 2 := ⟨m - 2, by omega⟩
  simp [cylQ3, cylQ3Prod, cylQ3Block, q3, add_assoc]

theorem cylQ3Prod_eq_one_iff {m : ℕ} {x : Vec (2 * m)} :
    cylQ3Prod m x = 1 ↔ ∀ j : Fin (2 * m), j.val < m - 2 → x j = 1 := by
  rw [cylQ3Prod, prod_range_eq_one_iff]
  refine ⟨fun h j hj ↦ by simpa [coord_of_lt _ _ j.isLt] using h j hj, fun h j hj ↦ ?_⟩
  rw [coord_of_lt _ _ (by omega)]
  exact h _ hj

theorem cylQ3Prod_add {m : ℕ} {x y : Vec (2 * m)}
    (hy : ∀ j : Fin (2 * m), j.val < m - 2 → y j = 0) :
    cylQ3Prod m (x + y) = cylQ3Prod m x :=
  prod_congr rfl fun j hj ↦ by
    have hj : j < m - 2 := mem_range.mp hj
    rw [coord_of_lt _ _ (by omega), coord_of_lt _ _ (by omega), Pi.add_apply, hy _ hj, add_zero]

theorem cylQ3Block_add (m : ℕ) (hm : 4 ≤ m) (x y : Vec (2 * m)) :
    cylQ3Block m (x + y) = cylQ3Block m x + cylQ3Block m y := by
  funext i
  simp only [cylQ3Block, Pi.add_apply, coord_of_lt _ _ (show m - 2 + i.val < 2 * m by omega)]

/-- The support of `cylQ3 m` lies in the coset `{x₀ = … = x_{m-3} = 1}`. -/
theorem cylQ3_supp (m : ℕ) (hm : 4 ≤ m) :
    ∀ x, cylQ3 m x = 1 → ∀ i ∈ prodCoords m, x i = (fun _ ↦ 1) i := fun x hx i hi ↦
  cylQ3Prod_eq_one_iff.1 ((zmod2_mul_eq_one _ _).1 ((cylQ3_eq_mul m hm x).symm.trans hx)).1 i
    (by simpa [prodCoords] using hi)

theorem cylQ3_zero (m : ℕ) (hm : 3 ≤ m) : cylQ3 m 0 = 0 := by
  rw [cylQ3, prod_eq_zero (i := 0) (mem_range.mpr (by omega)) (by simp [coord]), zero_mul]

theorem cylQ3_one (m : ℕ) (hm : 4 ≤ m) : cylQ3 m 1 = 1 := by
  rw [cylQ3_eq_mul m hm, (cylQ3Prod_eq_one_iff (x := 1)).2 fun _ _ ↦ rfl,
    show cylQ3Block m 1 = 1 from funext fun i ↦ coord_of_lt _ _ (by omega)]
  rfl

/-- The point whose product part is `1` and whose `q3` block is `u`. -/
def cosetPoint (m : ℕ) (u : Vec 6) : Vec (2 * m) := fun j ↦
  if j.val < m - 2 then 1 else coord u (j.val - (m - 2))

theorem cylQ3Prod_cosetPoint (m : ℕ) (u : Vec 6) : cylQ3Prod m (cosetPoint m u) = 1 :=
  cylQ3Prod_eq_one_iff.2 fun j hj ↦ by simp [cosetPoint, hj]

theorem cylQ3Block_cosetPoint (m : ℕ) (hm : 4 ≤ m) (u : Vec 6) :
    cylQ3Block m (cosetPoint m u) = u := by
  funext i
  have := i.isLt
  simp [cylQ3Block, cosetPoint, coord, show m - 2 + i.val < 2 * m by omega]

/-- A constant derivative `D_a (cylQ3 m) = c` forces `c = 0` and `a` to vanish on the
coordinates `0, …, m+3` that `cylQ3 m` depends on. -/
theorem cylQ3_linearStructure (m : ℕ) (hm : 4 ≤ m) {a : Vec (2 * m)} {c : ZMod 2}
    (h : ∀ z, dirDeriv a (cylQ3 m) z = c) :
    c = 0 ∧ ∀ j : Fin (2 * m), j.val < m + 4 → a j = 0 := by
  simp only [dirDeriv] at h
  have hF := cylQ3_eq_mul m hm
  have hprod_one : ∀ {x}, cylQ3 m x = 1 → cylQ3Prod m x = 1 := fun hx ↦
    ((zmod2_mul_eq_one _ _).1 ((hF _).symm.trans hx)).1
  -- `c = cylQ3 m a`; if it were `1`, then `D_a (cylQ3 m)` would vanish at `e₀`.
  have hc0 : c = 0 := by
    by_contra hc
    have ha := cylQ3Prod_eq_one_iff.1 <| hprod_one <| by
      simpa [cylQ3_zero m (by omega), zmod2_eq_one_of_ne_zero hc] using h 0
    have h0 : cylQ3Prod m (Pi.single ⟨0, by omega⟩ 1) = 0 :=
      prod_eq_zero (i := 1) (mem_range.mpr (by omega))
        (by rw [coord_of_lt _ _ (by omega)]; exact Pi.single_eq_of_ne (by simp) _)
    have h1 : cylQ3Prod m (Pi.single ⟨0, by omega⟩ 1 + a) = 0 :=
      prod_eq_zero (i := 0) (mem_range.mpr (by omega)) (by
        rw [coord_of_lt _ _ (by omega), Pi.add_apply, ha _ (by simp; omega)]
        exact (congrArg (· + 1) (Pi.single_eq_same _ _)).trans rfl)
    have := h (Pi.single ⟨0, by omega⟩ 1)
    rw [hF, hF, h0, h1] at this
    exact hc (by simpa using this.symm)
  have hper : ∀ z, cylQ3 m (z + a) = cylQ3 m z := fun z ↦
    (CharTwo.add_eq_zero.1 (hc0 ▸ h z)).symm
  -- The product coordinates of `a` vanish: compare at the all-ones vector.
  have hprod : ∀ j : Fin (2 * m), j.val < m - 2 → a j = 0 := fun j hj ↦ by
    have := cylQ3Prod_eq_one_iff.1 (hprod_one ((hper 1).trans (cylQ3_one m hm))) j hj
    simpa using this
  -- The block of `a` is a period of `q3`, hence `0`.
  have hblock : cylQ3Block m a = 0 := isBent_q3.eq_zero_of_periodic fun u ↦ by
    have := hper (cosetPoint m u)
    rwa [hF, hF, cylQ3Prod_add hprod, cylQ3Block_add m hm, cylQ3Prod_cosetPoint,
      cylQ3Block_cosetPoint m hm, one_mul, one_mul] at this
  refine ⟨hc0, fun j hj ↦ ?_⟩
  rcases lt_or_ge j.val (m - 2) with h' | h'
  · exact hprod j h'
  · have := congrFun hblock ⟨j.val - (m - 2), by omega⟩
    simpa [cylQ3Block, coord, show m - 2 + (j.val - (m - 2)) = j.val by omega] using this

/-- Adding a vector supported on the coordinates `m+4, …, 2m-1` does not change `cylQ3 m`. -/
theorem cylQ3_add_of_vanish (m : ℕ) (hm : 4 ≤ m) (z a : Vec (2 * m))
    (ha : ∀ j : Fin (2 * m), j.val < m + 4 → a j = 0) : cylQ3 m (z + a) = cylQ3 m z := by
  have hb : cylQ3Block m a = 0 := funext fun i ↦ by
    rw [cylQ3Block, coord_of_lt _ _ (by omega)]; exact ha _ (by simp; omega)
  rw [cylQ3_eq_mul m hm, cylQ3_eq_mul m hm, cylQ3Prod_add fun j hj ↦ ha j (by omega),
    cylQ3Block_add m hm, hb,
    add_zero]

/-! ## The counterexample -/

/-- For `m ≥ 4`, `cylQ3 m` is not a sum of two bent functions. -/
theorem not_isBentSum_cylQ3 (m : ℕ) (hm : 4 ≤ m) : ¬ IsBentSum (cylQ3 m) := by
  intro h
  have hdvd := h.two_pow_dvd_wt (k := m - 2) (by omega) (card_prodCoords m hm)
    (cylQ3_supp m hm)
  rw [wt_cylQ3 m hm, pow_succ, mul_comm 7] at hdvd
  have : 2 ∣ 7 := Nat.dvd_of_mul_dvd_mul_left (Nat.two_pow_pos (m - 2)) hdvd
  omega

/-- The bent-sum conjecture fails in every even dimension `2m ≥ 8`. -/
theorem not_bentSumConjecture (m : ℕ) (hm : 4 ≤ m) : ¬ BentSumConjecture (2 * m) :=
  fun h ↦ not_isBentSum_cylQ3 m hm
    (h _ (by rw [show 2 * m / 2 = m by omega]; exact inRM_cylQ3 m hm))

/-- The bent-sum conjecture fails in dimension `8`. -/
theorem not_bentSumConjecture_eight : ¬ BentSumConjecture 8 := not_bentSumConjecture 4 le_rfl

end BentSum
