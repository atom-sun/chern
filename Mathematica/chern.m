chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; 
      Assert[HermitianMatrixQ[hk @@ {0, 0}], 
       "Not Hermittian Matrix at {0,0}"]; Assert[AllTrue[hk @@ {0, 0}, 
        NumericQ, 2], "NumericQ fails at {0,0}"]; 
      {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
Attributes[Assert] = {HoldAllComplete}
 
Assert /: Assert::asrtf = "Assertion `1` failed."
 
Assert /: Assert::asrtfe = "Assertion `1` in `2` failed."
 
Assert /: Assert::asrtfl = "Assertion `1` at line `2` in `3` failed."
 
Assert /: Assert::asrttf = 
     "Assertion test `1` evaluated to `2` that is neither True nor False."
 
Assert /: Assert::usage = "\!\(\*RowBox[{\"Assert\", \"[\", \
StyleBox[\"test\", \"TI\"], \"]\"}]\) represents the assertion that \
\!\(\*StyleBox[\"test\", \"TI\"]\) is True. If assertions have been enabled, \
\!\(\*StyleBox[\"test\", \"TI\"]\) is evaluated when the assertion is \
encountered. If \!\(\*StyleBox[\"test\", \"TI\"]\) is not True, then an \
assertion failure is generated.\n\!\(\*RowBox[{\"Assert\", \"[\", \
RowBox[{StyleBox[\"test\", \"TI\"], \",\", StyleBox[\"tag\", \"TI\"]}], \
\"]\"}]\) specifies a tag that will be used to identify the assertion if it \
fails."
 
diagn[hk_, qqx_, qqy_] := Module[{nqqx = Length[qqx], nqqy = Length[qqy], kx, 
      ky, hErg, hVec, hamilt, erg, vec, where, Null}, 
     hErg = Table[Null, {jqy, 1, nqqy}, {jqx, 1, nqqx}]; 
      hVec = Table[Null, {jqy, 1, nqqy}, {jqx, 1, nqqx}]; 
      Table[kx = qqx[[jqx]]; ky = qqx[[jqy]]; hamilt = hk @@ {kx, ky}; 
        {erg, vec} = Eigensystem[N[hamilt]]; where = DeleteDuplicates[
          Flatten[(Position[erg, #1] & ) /@ Sort[erg, Less]]]; 
        erg = erg[[where]]; vec = vec[[where]]; hErg[[jqy,jqx]] = erg; 
        hVec[[jqy,jqx]] = vec; , {jqy, nqqy}, {jqx, nqqx}]; {hErg, hVec}]
 
calcChern[hVec_] := Module[{cnlist = Table[0, Dimensions[hVec][[3]]], u12, 
      u23, u34, u41, t12, t23, t34, t41, tplaquet, 
      nqqy = Dimensions[hVec][[1]], nqqx = Dimensions[hVec][[2]], Null}, 
     Table[u12 = Conjugate[hVec[[jqy,jqx + 1]]] . Transpose[hVec[[jqy,jqx]]]; 
        u23 = Conjugate[hVec[[jqy + 1,jqx + 1]]] . Transpose[
           hVec[[jqy,jqx + 1]]]; u34 = Conjugate[hVec[[jqy + 1,jqx]]] . 
          Transpose[hVec[[jqy + 1,jqx + 1]]]; 
        u41 = Conjugate[hVec[[jqy,jqx]]] . Transpose[hVec[[jqy + 1,jqx]]]; 
        t12 = DiagonalMatrix[Diagonal[u12]]; 
        t23 = DiagonalMatrix[Diagonal[u23]]; 
        t34 = DiagonalMatrix[Diagonal[u34]]; 
        t41 = DiagonalMatrix[Diagonal[u41]]; tplaquet = 
         t41 . t34 . t23 . t12; cnlist += Arg[Diagonal[tplaquet]]; , 
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; cnlist]
