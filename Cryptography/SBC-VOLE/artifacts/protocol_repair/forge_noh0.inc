/* forge_noh0.inc -- the cheating prover for the MIS-BOUND masked build.
 *
 * Compiled in only with -DFORGE_NOH0, and only meaningful together with
 * -DMASKB -DMASKB_NOH0 (the mask offsets dw_j are transmitted and enter the
 * verifier's check, but are NOT hashed into h_0).
 *
 * Why the scheme collapses.  dw_j enters the verifier's reconstruction of the
 * mask row only through the side-1 shares it patches, so it contributes
 *
 *     Q_j(dw_j) = Q_j(0) + dw_j * RoundDelta_j ,   RoundDelta_j = sum_{i: b_{i,j}=0} alpha_{j,i}
 *
 * to Qtot, and Qtot enters A'' linearly with coefficient 1.  dw_j enters
 * nothing else: not h_0 (by hypothesis), not h_com, not A~ = A + Mtot (which
 * depends only on the leaf expansions), not B~ = B + W (which depends on
 * dw_0 = 0 alone).  So the challenge (i*_j) and hence Delta* and RoundDelta_j
 * do not depend on dw_1..dw_{tau-1}, and after seeing the challenge the prover
 * solves one linear equation in F_{2^257}:
 *
 *     dw_j  <-  dw_j + (A~ + A''(0)) * RoundDelta_j^{-1}
 *
 * for any j >= 1 with RoundDelta_j != 0.  One field inversion; the Fiat-Shamir
 * check becomes vacuous.  Residual failure probability is exactly the
 * probability that every RoundDelta_j (j >= 1) vanishes, i.e. 2^{-D(tau-1)}
 * (2^-117 at the paper's (13,10)).
 *
 * Everything this block uses is prover-side knowledge: the prover holds all the
 * leaves, so it can evaluate precisely what the verifier will evaluate.  No
 * secret of the *victim* is used -- FORGE_NOH0 is run with a FALSE witness
 * against the victim's public key.
 */
#ifdef FORGE_NOH0
  {
    u64 *FHyp, *FXO, *FYO, *FYR, *FDelta;
    u64 FScalars[35];
    u64 resid[5], rdel[5], rinv[5], prod[5];
    SHA256_CTX fctx;
    int jj, bb, done, q;

    FHyp = (u64*)malloc(sizeof(u64)*2*PLAYER_SHARE_MEM*NbRounds*LogPlayers);
    FXO  = (u64*)malloc(sizeof(u64)*5*Lambda);
    FYO  = (u64*)malloc(sizeof(u64)*5*Lambda);
    FYR  = (u64*)malloc(sizeof(u64)*5*Lambda);
    FDelta = (u64*)malloc(sizeof(u64)*5*NbRounds);

    /* 1. Replicate the verifier's share patching on a scratch copy.
     *    (CheckVOLESignature, "patch hypershares so that VOLE is correctly
     *    evaluated": dy for every round, dz and dw for every round but 0.) */
    memcpy(FHyp,HyperMemorys,sizeof(u64)*2*PLAYER_SHARE_MEM*NbRounds*LogPlayers);
    for(jj=0;jj<NbRounds;jj++) {
      for(bb=0;bb<LogPlayers;bb++) {
        for(q=0;q<2;q++)
          FHyp[PLAYER_SHARE_MEM*(2*(jj*LogPlayers+bb)+1)+yOFFSET+q]^=
            Offsets[jj*PLAYER_OFFSET_MEM+yOFFSET+q];
        if (jj!=0) {
          for(q=0;q<5;q++)
            FHyp[PLAYER_SHARE_MEM*(2*(jj*LogPlayers+bb)+1)+zOFFSET+q]^=
              Offsets[jj*PLAYER_OFFSET_MEM+zOFFSET+q];
          for(q=0;q<5;q++)
            FHyp[PLAYER_SHARE_MEM*(2*(jj*LogPlayers+bb)+1)+mOFFSET+q]^=
              Offsets[jj*PLAYER_OFFSET_MEM+mOFFSET+q];
        }
      }
    }

    /* 2. Evaluate the verifier's own arithmetic to get A''(0). */
    for(q=0;q<35;q++) FScalars[q]=Scalars[q];
    SHA256Init(&fctx);
    DoOLEevaluateCheck(&fctx,Queries,UV,FHyp,FXO,FYO,FYR,FScalars,Tvec);

    /* 3. residual = A~ + A''(0).  A~ is MaskAtilde, stashed by the signer. */
    for(q=0;q<5;q++) resid[q]=MaskAtilde[q]^FScalars[20+q];

    /* 4. RoundDelta_j for every j, then one inversion. */
    done=0;
    for(jj=0;jj<NbRounds;jj++) {
      for(q=0;q<5;q++) FDelta[5*jj+q]=0;
      for(bb=0;bb<LogPlayers;bb++)
        if ((((Queries[jj]>>bb)&1)^1))
          for(q=0;q<5;q++) FDelta[5*jj+q]^=Tvec[5*Lambda+5*(jj*LogPlayers+bb)+q];
    }
    for(jj=1;jj<NbRounds && !done;jj++) {
      int nonzero=0;
      for(q=0;q<5;q++) if (FDelta[5*jj+q]) nonzero=1;
      if (!nonzero) continue;
      FFInv(rinv,FDelta+5*jj);
      FFMul(prod,resid,rinv);
      for(q=0;q<5;q++) Offsets[jj*PLAYER_OFFSET_MEM+mOFFSET+q]^=prod[q];
      done=1;
    }
    if (getenv("FORGE_NOH0"))
      fprintf(stderr,"FORGE_NOH0: forged via tree j=%d of %d ; every Delta*_j (j>=1) zero: %s\n",
              done?jj-1:-1,NbRounds,done?"no":"YES (this instance is unfixable)");
    free(FHyp); free(FXO); free(FYO); free(FYR); free(FDelta);
  }
#endif
