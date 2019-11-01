(* Mathematica Source File *)
(* Chern Calculator tutorial *)
(* :Author: Ning Sun *)
(* :Date: 2019-11-01 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 atom-sun *)

Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
<< "../src/Mathematica/chern.m";


(* show context path *)
$ContextPath


(* show namespace *)
Names["ChernCalculator`*"]


(* show usage or docstring *)
?chern


(* Preferred input form *)
hk[kx_, ky_] := ( {
   {1 - Cos[kx] - Cos[ky], Sin[kx] - I Sin[ky]},
   {Sin[kx] + I Sin[ky], -1 + Cos[kx] + Cos[ky]}
  } );

chern[hk, 64]


(* Not preferred input form *)
hk1 := ( {
  {1 - Cos[kx] - Cos[ky], Sin[kx] - I Sin[ky]},
  {Sin[kx] + I Sin[ky], -1 + Cos[kx] + Cos[ky]}
} );

chern[hk1, 64]
