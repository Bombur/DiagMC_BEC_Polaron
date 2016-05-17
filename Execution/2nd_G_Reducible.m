(* ::Package:: *)

(* ::Title:: *)
(*2nd Order reduciable Diagram in \[Tau] p=0*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<<params.mx;
ResetDirectory[];
p = 0;
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha] = Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ];
taus = Import["data/zero_order", "table"][[All,1]];
FP = Import["DiagMC_BEC.json", {"Data", "Froehlich_Polaron"} ];


(* ::Chapter:: *)
(*t Integral*)


IntegFP[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_,\[Tau]3_,\[Tau]4_] := vq2FP[q1]*vq2FP[q2]*g0FP[0,0,\[Tau]1]*g0FP[-q1,\[Tau]1,\[Tau]2]*dqFP[\[Tau]1,\[Tau]2]*g0FP[0,\[Tau]2,\[Tau]3]*g0FP[-q2,\[Tau]3,\[Tau]4]*dqFP[\[Tau]3,\[Tau]4]*g0FP[0,\[Tau]4,\[Tau]];
IntegBEC[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_,\[Tau]3_,\[Tau]4_] := vq2[q1]*vq2[q2]*g0[0,0,\[Tau]1]*g0[-q1,\[Tau]1,\[Tau]2]*dq[q1,\[Tau]1,\[Tau]2]*g0[0,\[Tau]2,\[Tau]3]*g0[-q2,\[Tau]3,\[Tau]4]*dq[q2,\[Tau]3,\[Tau]4]*g0[0,\[Tau]4,\[Tau]];
f[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_,\[Tau]3_,\[Tau]4_]:=(4*Pi)^2/(2*Pi)^6 *q1^2*q2^2*If[FP,IntegFP[q1,q2,\[Tau],\[Tau]1, \[Tau]2,\[Tau]3,\[Tau]4], IntegBEC[q1,q2,\[Tau],\[Tau]1, \[Tau]2,\[Tau]3,\[Tau]4]];
g[q1_,q2_, \[Tau]_,\[Tau]2_,\[Tau]3_,\[Tau]4_] := Integrate[f[q1,q2,\[Tau],\[Tau]1,\[Tau]2,\[Tau]3,\[Tau]4],  {\[Tau]1,0,\[Tau]2}];
h[q1_,q2_, \[Tau]_,\[Tau]3_,\[Tau]4_] := Integrate[g[q1,q2,\[Tau],\[Tau]2,\[Tau]3,\[Tau]4],  {\[Tau]2,0,\[Tau]3}];
i[q1_,q2_, \[Tau]_,\[Tau]4_] := Integrate[h[q1,q2,\[Tau],\[Tau]3,\[Tau]4],  {\[Tau]3,0,\[Tau]4}];
j[q1_,q2_, \[Tau]_] := Integrate[i[q1,q2,\[Tau],\[Tau]4],  {\[Tau]4,0,\[Tau]}];
(*LogPlot[ {j[1, q2, 1], j[q2, 1, 1]}, {q2, 0,qc}, PlotPoints \[Rule] 10, MaxRecursion\[Rule]2, PlotRange -> Full]*)


(* ::Chapter:: *)
(*q Integral*)


g2[\[Tau]_?NumericQ]:=NIntegrate[j[q1, q2, \[Tau]],{q1,0, qc}, {q2,0, qc}, AccuracyGoal->100];
(*g2[1]
Plot[g2[t],{t, 0,1}, PlotPoints \[Rule] 10, MaxRecursion\[Rule]2, PlotRange -> Full]*)
g2data = Table[{t, g2[t]}, {t, taus}];


(* ::Chapter:: *)
(*Export*)


Export["data/mat_2nd_G_reducible", g2data, "Table"];
