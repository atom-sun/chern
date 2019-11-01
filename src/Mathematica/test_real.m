(* Mathematica Source File *)
(* :Author: Ning Sun *)
(* :Date: 2019-11-01 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 atom-sun *)

Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
<< "chern.m";


(* All real Hamiltonian should have zero Chern number *)
hhh[kx_, ky_] := {{3 - Cos[kx] - Cos[ky], Sin[kx], Sin[ky]}, {Sin[kx],
   Cos[kx]*Cos[ky], Sin[kx]*Sin[ky]}, {Sin[ky], Sin[kx]*Sin[ky], 1}};

chern[hhh, 64]
