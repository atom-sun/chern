(* Mathematica Source File *)
(* :Author: Ning Sun *)
(* :Date: 2019-11-01 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 atom-sun *)

Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
<< "../src/Mathematica/chern.m";

(* haldane honeycomb *)
k := {kx/Sqrt[3], ky*2/3};
t1 = t2 = 1;
a1 = {Sqrt[3]/2, 1/2};
a2 = {0, -1};
a3 = {-Sqrt[3]/2, 1/2};
b1 := a2 - a3;
b2 := a3 - a1;
b3 := a1 - a2;
hk := 2*t2*Cos[\[Phi]]*(Cos[k.b1] + Cos[k.b2] + Cos[k.b3])*
   PauliMatrix[0] +
  t1*((Cos[k.a1] + Cos[k.a2] + Cos[k.a3])*
      PauliMatrix[1] + (Sin[k.a1] + Sin[k.a2] + Sin[k.a3])*
      PauliMatrix[2]) + (m -
     2*t2*Sin[\[Phi]]*(Sin[k.b1] + Sin[k.b2] + Sin[k.b3]))*
   PauliMatrix[3];

m = 0.5;
\[Phi] = \[Pi]/2;
hk // MatrixForm;

(* default discretization 16x16 *)
chern[hk]

(* discretized to 64x64 *)
chern[hk, 64]
