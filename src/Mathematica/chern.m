(* Mathematica Package *)

(* :Title: Chern Calculator for 2D materials *)
(* :Context: ChernCalculator` *)
(* :Author: Ning Sun *)
(* :Date: 2019-11-01 *)

(* :Package Version: 0.3 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["ChernCalculator`"];
(* Exported symbols added here with SymbolName::usage *)
chern::usage =
    "Chern number calculator for 2D materials.
    Input hamiltonian should be either
    1) bivariate function that generates a Hermitian matrix, like
      hk[kx_, ky_] := ({{1 - Cos[kx] - Cos[ky], Sin[kx] - I Sin[ky]},
          {Sin[kx] + I Sin[ky], -1 + Cos[kx] + Cos[ky]}})
      or
    2) some set-delayed expression that contains 'kx' and 'ky' (who are
      actually Global`kx, Global`ky), like
      hk := ({{1 - Cos[kx] - Cos[ky], Sin[kx] - I Sin[ky]},
          {Sin[kx] + I Sin[ky], -1 + Cos[kx] + Cos[ky]}})
    The former is preferred.
    ";

Begin["`Private`"];

chern[hamiltonian_, ndiscretized_: 16] :=
    Module[
      {nQq = ndiscretized,
        hk = If[funcQ[hamiltonian], hamiltonian,
          (hamiltonian/.{Global`kx->#1, Global`ky->#2})&],
        qqx, qqy, hErg, hVec, cn
      },
      {qqx, qqy} = disCretize[nQq];
      If[
        check[hk] === 0,
        cn = ConstantArray[0, Length[hk @@ {0., 0.}]],
        {hErg, hVec} = diagn[hk, qqx, qqy]; cn = calcChern[hVec];
      ];
      cn
    ];

check[hamiltonian_] :=
    Block[
      {$AssertFunction = Throw[{##}]&, kx, ky, hk, ndim},
      hk = hamiltonian @@ {kx, ky};
      ndim = Length[hk];
      Assert[Dimensions[Dimensions[hk]] === {2}, "input not matrix"];
      Assert[Dimensions[hk][[2]] === ndim, "input not square matrix"];
      Assert[AllTrue[hk @@ {0, 0}, NumericQ, 2], "NumericQ fails at {0,0}"];
      Assert[Chop[ComplexExpand[hk- ConjugateTranspose[hk]]] ===
          ConstantArray[0, {ndim, ndim}], "Not Hermitian matrix"];
      If[Chop[ComplexExpand[hk - Transpose[hk]]] ===
          ConstantArray[0, {ndim, ndim}], 0, 1]
    ];

funcQ[_Function | _IntepolatingFunction | _CompiledFunction] := True;
funcQ[f_Symbol] :=
    Or[DownValues[f] =!= {}, MemberQ[Attributes[f], NumericFunction]];
funcQ[_] := False;

disCretize[nQq_] :=
    Module[
      {nQx = nQq, nQy = nQq, nqqx, dqx, qqx, nqqy, dqy, qqy},
      nqqx = nQx + 1;
      dqx = 2. * \[Pi] / nQx;
      qqx = N[Table[(jqx - 1.) * dqx - \[Pi], {jqx, 1, nqqx}]];
      nqqy = nQy + 1;
      dqy = 2. * \[Pi] / nQy;
      qqy = N[Table[(jqy - 1.) * dqy - \[Pi], {jqy, 1, nqqy}]];
      {qqx, qqy}
   ];

diagn[hk_, qqx_, qqy_] :=
    Module[
      {nqqx = Length[qqx], nqqy = Length[qqy],
        kx, ky, hErg, hVec, hamilt, erg, vec, where},
      hErg = Table[Null, {jqy, 1, nqqy}, {jqx, 1, nqqx}];
      hVec = Table[Null, {jqy, 1, nqqy}, {jqx, 1, nqqx}];
      Table[
        kx = qqx[[jqx]];
        ky = qqx[[jqy]];
        hamilt = hk @@ {kx, ky};
        {erg, vec} = Eigensystem[N[hamilt]];
        where = DeleteDuplicates[Flatten[Map[Position[erg, #]&,
          Sort[erg, Less]]]];
        erg = erg[[where]];
        vec = vec[[where]];
        Part[hErg, jqy, jqx] = erg;
        Part[hVec, jqy, jqx] = vec;
        , {jqy, nqqy}, {jqx, nqqx}];
      {hErg, hVec}
   ];

calcChern[hVec_] :=
 Module[
  {cnlist = Table[0, Dimensions[hVec][[3]]],
    u12, u23, u34, u41, t12, t23, t34, t41, tplaquet,
    nqqy = Dimensions[hVec][[1]], nqqx = Dimensions[hVec][[2]]},
   Table[
     u12 = Conjugate[hVec[[jqy, jqx + 1]]].Transpose[hVec[[jqy, jqx]]];
     u23 = Conjugate[hVec[[jqy + 1, jqx + 1]]].Transpose[hVec[[jqy, jqx + 1]]];
     u34 = Conjugate[hVec[[jqy + 1, jqx]]].Transpose[hVec[[jqy + 1, jqx + 1]]];
     u41 = Conjugate[hVec[[jqy, jqx]]].Transpose[hVec[[jqy + 1, jqx]]];
     t12 = DiagonalMatrix[Diagonal[u12]];
     t23 = DiagonalMatrix[Diagonal[u23]];
     t34 = DiagonalMatrix[Diagonal[u34]];
     t41 = DiagonalMatrix[Diagonal[u41]];
     tplaquet = t41.t34.t23.t12;
     cnlist += Arg[Diagonal[tplaquet]];
     , {jqy, 1, nqqy-1}, {jqx, 1, nqqx-1}];
   cnlist /= 2 \[Pi];
   Chop[cnlist]
 ];

End[]; (* `Private` *)

EndPackage[]
