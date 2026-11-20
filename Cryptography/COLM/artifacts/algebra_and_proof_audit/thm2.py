#!/usr/bin/env python3
"""
thm2.py -- Quantitative evaluation of BOTH official security proofs for COLM_0:
  (1) The 2024 Journal of Cryptology full paper ("The COLM Authenticated
      Encryption Scheme", J. Cryptol. 37:15, 2024, Section 4, Theorem 2 &
      Lemma 2, pp. 15-25);
  (2) The 2016 CAESAR Round-3 security proof addendum ("Security of COLM",
      Section 2, Theorem 2 & Lemma 3, p. 3);
evaluated at the exact parameters of our two-query nonce-respecting,
encryption-only splice forgery adversary (n = 128, q = 2 forward encryption
queries, s = 1 forgery attempt).

1. 2024 JOURNAL OF CRYPTOLOGY PROOF (J. Cryptol. 37:15, 2024, p. 15):
   Theorem 2 states verbatim:
     "Let us consider an adversary A that makes q encryption queries and s
      forging attempts at an overall aggregate of total a associated data and
      m message blocks and tries to distinguish COLM_0 from a random online
      function. Let sigma = a + 2m + 3(q + s). The oAE advantage of COLM_0 is
      given by
          Adv^oAE_{COLM_0}(A) <= Adv^SPRP_E(A')
                               + (7*sigma^2 + 3*sigma + q*s + 3*s + 2*s*sigma) / 2^n"
   where Lemma 2 (p. 23) bounds bad transcripts by
          epsilon_bad <= (7*sigma^2 + 3*sigma + q*s) / 2^n
   and Section 4.4 (pp. 24-25) bounds good-transcript forgery ratio by
          epsilon_ratio <= (3*s + 2*s*sigma) / 2^n.

2. 2016 CAESAR ROUND-3 ADDENDUM PROOF (p. 3):
   Theorem 2 states verbatim:
     "Let sigma = sum_{i=1}^{q+s} (a_i + l_i) and sigma_auth = sigma + q. Then
          Adv^auth_{COLM_Pi}(A) <= 9 * sigma_auth^2 / 2^n + s / 2^n"
   with Lemma 3 (p. 3) bounding bad views by <= 2 * sigma_auth^2 / 2^n.
"""
import math

n = 128
ok = True


def joc2024_thm2(q, s, a_list, l_list):
    """2024 Journal of Cryptology Theorem 2 & Lemma 2 bounds (pp. 15, 23, 25)."""
    a = sum(a_list)
    m = sum(l_list)
    sigma = a + 2 * m + 3 * (q + s)
    eps_bad = (7.0 * sigma ** 2 + 3.0 * sigma + q * s) / 2.0 ** n
    eps_ratio = (3.0 * s + 2.0 * s * sigma) / 2.0 ** n
    total = eps_bad + eps_ratio
    return a, m, sigma, eps_bad, eps_ratio, total


def addendum2016_thm2(q, s, a_list, l_list):
    """2016 CAESAR Round-3 Addendum Theorem 2 & Lemma 3 bounds (p. 3)."""
    sigma = sum(a_list) + sum(l_list)
    sigma_auth = sigma + q
    t1 = (q + sum(l_list[:q])) ** 2 / 2.0 ** (n + 1)
    t2 = s * (q + s + sigma) / 2.0 ** n
    lem3 = t1 + t2
    thm2 = (9.0 * sigma_auth ** 2 + s) / 2.0 ** n
    return sigma, sigma_auth, lem3, thm2


print("=" * 80)
print("Evaluation of the 2024 Journal of Cryptology Proof (Theorem 2, J. Cryptol. 2024)")
print("and 2016 CAESAR Addendum Proof at q = 2 encryption queries, s = 1 forgery:")
print("=" * 80)
print()

joc_rows = []
add_rows = []
for L, name in ((130, "fade       l = 130 blocks"),
                (129, "transplant l = 129 blocks")):
    for ab in (0, 1, 2):
        # In J. Cryptol. 2024 Sec. 4, the first AD block is (npub || param), so
        # a_i = max(1, ab) + 1 or ab blocks depending on whether (npub||param) is counted in a_i.
        # Under Section 4.1 (p. 16, line 837), each query has a_i AD blocks plus (npub||param)
        # accounted in +3(q+s) of sigma = a + 2m + 3(q+s).
        a_list = [ab] * 3
        l_list = [L] * 3
        a_tot, m_tot, sig24, ebad24, eratio24, tot24 = joc2024_thm2(2, 1, a_list, l_list)
        sig16, sa16, lem3_16, tot16 = addendum2016_thm2(2, 1, a_list, l_list)
        joc_rows.append((name, ab, sig24, ebad24, eratio24, tot24))
        add_rows.append((name, ab, sig16, sa16, lem3_16, tot16))
        print("%-26s  AD = %d block(s)/query" % (name, ab))
        print("  [2024 J. Cryptol. Thm 2] sigma = %d (a=%d, m=%d)" % (sig24, a_tot, m_tot))
        print("      eps_bad (Lemma 2)   = 2^%.2f,  eps_ratio (Sec 4.4) = 2^%.2f"
              % (math.log2(ebad24), math.log2(eratio24)))
        print("      Theorem 2 bound     = 2^%.2f   vs. measured Adv^oAE  = 1.0 (exceeded by 2^%.1f)"
              % (math.log2(tot24), -math.log2(tot24)))
        print("  [2016 Addendum Thm 2]    sigma = %d, sigma_auth = %d" % (sig16, sa16))
        print("      Theorem 2 bound     = 2^%.2f   vs. measured Adv^auth = 1.0 (exceeded by 2^%.1f)"
              % (math.log2(tot16), -math.log2(tot16)))
        print()

best24 = max(joc_rows, key=lambda r: r[5])
best16 = max(add_rows, key=lambda r: r[5])
print("Summary at most favourable parameter choice (fade, l = 130, AD = 2 blocks/query):")
print("  2024 J. Cryptol. Theorem 2 bound : 2^%.2f  (exceeded by 2^%.2f)"
      % (math.log2(best24[5]), -math.log2(best24[5])))
print("  2016 Addendum Theorem 2 bound    : 2^%.2f  (exceeded by 2^%.2f)"
      % (math.log2(best16[5]), -math.log2(best16[5])))
print()

print("=" * 80)
print("Locating the exact flaws in the 2024 Journal of Cryptology proof (pp. 20-25)")
print("=" * 80)
print("  Flaw 1 (Section 4.4, p. 24, BadSamp & eps_ratio):")
print("    In Section 4.4 (p. 24), for a good transcript where {j : X_i^-[j] is fresh}")
print("    is non-empty with smallest index j_i < l_i, the proof samples MM_i^-[j_i]")
print("    randomly and claims that Pr[BadSamp] = Pr[MM_i^-[l_i] not fresh] <= 1/(2^n-2*sigma)")
print("    'due to the randomness of MM_i^-[j_i]', and then samples CC_i^-[l_i+1] uniformly.")
print("    In reality, during COLM_0.dec (Fig. 3, lines 6-18), MM_i^-[l_i] = E_K^-1(X_i^-[l_i])")
print("    is determined solely by C_i^-[l_i] and W_i^-[l_i-1]; it does NOT depend on MM_i^-[j_i]!")
print("    When D[129] = 0 (Fade), W_i^-[l_i-1] = W_1^+[l_i-1] and C_i^-[l_i] = C_1^+[l_i], so")
print("    X_i^-[l_i] == X_1^+[l_i], MM_i^-[l_i] == MM_1^+[l_i], and CC_i^-[l_i+1] == CC_1^+[l_i+1]")
print("    hold deterministically with probability 1 regardless of MM_i^-[j_i].")
print()
print("  Flaw 2 (Lemma 1, p. 20 & Bad Event B6, pp. 22-23):")
print("    In Transplant (l = 129), every block of C_i^- comes from either Query 1 or Query 2,")
print("    so no CC_i^-[j] is fresh. Bad Event B6 (p. 23) invokes Lemma 1 (p. 20) to claim that")
print("    L_X(W[0], Y[1..j+1]) == X'[j+1] forces W[0] == W'[0], overlooking that when two")
print("    queries share a 129-block plaintext prefix (X_1^+[r] == X_2^+[r]) and j+1 >= 129,")
print("    the GF(2^128) relation 3 * sum_{r in {1,3,4,5,6,7}} 2^(r-1) = 2^128 mod f(x)")
print("    transfers W_i^-[128] onto W_2^+[128] even though W_1^+[0] != W_2^+[0]!")
print()

checks = [
    ("2024 J. Cryptol. sigma = a + 2m + 3(q+s) = 795 at a=6, m=390, q=2, s=1",
     joc2024_thm2(2, 1, [2, 2, 2], [130, 130, 130])[2] == 795),
    ("2024 J. Cryptol. Theorem 2 bound at sigma=795 is 2^-105.92",
     abs(math.log2(joc2024_thm2(2, 1, [2, 2, 2], [130, 130, 130])[5]) + 105.92) < 0.01),
    ("2016 Addendum sigma_auth = 398 and Theorem 2 bound is 2^-107.56",
     abs(math.log2(addendum2016_thm2(2, 1, [2, 2, 2], [130, 130, 130])[3]) + 107.56) < 0.01),
]
for name, cond in checks:
    print("  %-72s %s" % (name, "OK" if cond else "FAILED"))
    ok = ok and cond

print()
print("ARITHMETIC: %s" % ("ALL CHECKS PASS" if ok else "FAILURE"))
raise SystemExit(0 if ok else 1)
