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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk[{0., 0.}]]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[FullSimplify[hamkxky - ConjugateTranspose[
             hamkxky]]]] === ConstantArray[0, {ndim, ndim}], 
       "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[FullSimplify[hamkxky - Transpose[hamkxky]]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
Attributes[Assert] = {HoldAllComplete}
 
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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk[{0., 0.}]]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; {hErg, hVec} = 
       diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_] := Module[{nQq = ndiscretized, 
      hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
hk := 2*t2*Cos[\[Phi]]*(Cos[k . b1] + Cos[k . b2] + Cos[k . b3])*
      PauliMatrix[0] + t1*((Cos[k . a1] + Cos[k . a2] + Cos[k . a3])*
        PauliMatrix[1] + (Sin[k . a1] + Sin[k . a2] + Sin[k . a3])*
        PauliMatrix[2]) + (m - 2*t2*Sin[\[Phi]]*(Sin[k . b1] + Sin[k . b2] + 
         Sin[k . b3]))*PauliMatrix[3]
 
t2 = 1
 
\[Phi] = Pi/2
 
k := {kx/Sqrt[3], ky*(2/3)}
 
b1 := a2 - a3
 
a2 = {0, -1}
 
a3 = {-Sqrt[3]/2, 1/2}
 
b2 := a3 - a1
 
a1 = {Sqrt[3]/2, 1/2}
 
b3 := a1 - a2
 
t1 = 1
 
m = 0.5
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
$AssertFunction /: $AssertFunction::usage = 
     "$AssertFunction specifies a function to apply to assertions that fail."
 
hamkxky = {{0.5 - 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky]), 
      Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] - 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3])}, 
     {Cos[kx/2 - ky/3] + Cos[kx/2 + ky/3] + Cos[(2*ky)/3] + 
       I*(-Sin[kx/2 - ky/3] + Sin[kx/2 + ky/3] - Sin[(2*ky)/3]), 
      -0.5 + 2*(-Sin[kx] + Sin[kx/2 - ky] + Sin[kx/2 + ky])}}
 
ndim = 2
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
chern[hamiltkxky_, ndiscretized_:16] := 
    Module[{nQq = ndiscretized, hk = If[funcQ[hamiltkxky], hamiltkxky, 
        hamiltkxky /. {kx -> #1, ky -> #2} & ], qqx, qqy, hErg, hVec, 
      cnlist}, {qqx, qqy} = disCretize[nQq]; If[check[hk] === 0, 
       cnlist = ConstantArray[0, Length[hk @@ {0., 0.}]], 
       {hErg, hVec} = diagn[hk, qqx, qqy]; cnlist = calcChern[hVec]; ]; 
      cnlist]
 
funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True
 
funcQ[f_Symbol] := DownValues[f] =!= {} || MemberQ[Attributes[f], 
      NumericFunction]
 
funcQ[_] := False
 
disCretize[nQq_] := Module[{nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, 
      qqy}, nqqx = nQx + 1; dqx = 2.*(Pi/nQx); 
      qqx = N[Table[(jqx - 1.)*dqx - Pi, {jqx, 1, nqqx}]]; nqqy = nQy + 1; 
      dqy = 2.*(Pi/nQy); qqy = N[Table[(jqy - 1.)*dqy - Pi, {jqy, 1, nqqy}]]; 
      {qqx, qqy}]
 
check[hk_] := Block[{$AssertFunction = Throw[{##1}] & , kx, ky, hamkxky, 
      ndim}, hamkxky = hk @@ {kx, ky}; ndim = Length[hamkxky]; 
      Assert[Dimensions[Dimensions[hamkxky]] === {2}, "input not matrix"]; 
      Assert[Dimensions[hamkxky][[2]] === ndim, "input not square matrix"]; 
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"]; 
      Assert[Chop[ComplexExpand[hamkxky - ConjugateTranspose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], "Not Hermittian matrix"]; 
      If[Chop[ComplexExpand[hamkxky - Transpose[hamkxky]]] === 
        ConstantArray[0, {ndim, ndim}], 0, 1]]
 
Attributes[Assert] = {HoldAllComplete}
 
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
       {jqy, 1, nqqy - 1}, {jqx, 1, nqqx - 1}]; cnlist /= 2*Pi; Chop[cnlist]]
