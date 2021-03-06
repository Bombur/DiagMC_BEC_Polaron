(* ::Package:: *)

(* ::Title:: *)
(*2nd Order Rainbow G0SE in \[Tau] p=0*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<<params.mx;
ResetDirectory[];
FP = Import["DiagMC_BEC.json", {"Data", "Froehlich_Polaron"} ];
p = 0;
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha] = Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = If[FP, 100,Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ]];
taus = Import["data/zero_order", "table"][[All,1]];


(* ::Chapter:: *)
(*Rainbow*)


IntegFP[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_,\[Tau]3_, ctheta_] := vq2FP[q1]*vq2FP[q2]*g0FP[0,0,\[Tau]1]*g0FP[q1,\[Tau]1,\[Tau]2]*dqFP[\[Tau]1,\[Tau]]*g0cosFP[q1,q2,\[Tau]2,\[Tau]3, ctheta]*dqFP[\[Tau]2,\[Tau]3]*g0FP[q1,\[Tau]3,\[Tau]];
IntegBEC[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_, \[Tau]2_,\[Tau]3_, ctheta_] := vq2[q1]*vq2[q2]*g0[0,0,\[Tau]1]*g0[q1,\[Tau]1,\[Tau]2]*dq[q1,\[Tau]1,\[Tau]]*g0cos[q1,q2,\[Tau]2,\[Tau]3, ctheta]*dq[q2,\[Tau]2,\[Tau]3]*g0[q1,\[Tau]3,\[Tau]];
f[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_,\[Tau]3_, ctheta_]:=4*Pi*2*Pi/(2*Pi)^6 *q1^2*q2^2*If[FP,IntegFP[q1,q2,\[Tau],\[Tau]1, \[Tau]2, \[Tau]3, ctheta], IntegBEC[q1,q2,\[Tau],\[Tau]1, \[Tau]2, \[Tau]3,ctheta]];
g[q1_,q2_, \[Tau]_,\[Tau]2_,\[Tau]3_,ctheta_] := Integrate[f[q1,q2,\[Tau],\[Tau]1,\[Tau]2,\[Tau]3,ctheta],  {\[Tau]1,0,\[Tau]2}];
h[q1_,q2_, \[Tau]_,\[Tau]3_, ctheta_] := Integrate[g[q1,q2,\[Tau],\[Tau]2,\[Tau]3, ctheta],  {\[Tau]2,0,\[Tau]3}];
i[q1_,q2_, \[Tau]_, ctheta_] := Integrate[h[q1,q2,\[Tau],\[Tau]3, ctheta],  {\[Tau]3,0,\[Tau]}];
(*Integrate[h[q1,q2,\[Tau],\[Tau]3, ctheta],  {\[Tau]3,0,\[Tau]}]*)
(*LogPlot[ {h[1, q2, 1, 0], h[q2, 1, 1, 0]}, {q2, 0,qc}, PlotPoints \[Rule] 10, MaxRecursion\[Rule]2, PlotRange -> Full]*)


g0se2rainb[\[Tau]_?NumericQ]:=NIntegrate[i[q1, q2, \[Tau], ctheta],{q1,0, qc}, {q2,0, qc}, {ctheta, -1,1}, AccuracyGoal->100];
(*LogPlot[se2rainb[t],{t, 0,10}, PlotPoints -> 10, MaxRecursion\[Rule]3, PlotRange -> Full]*)
(*g0se2rainb[1]*)


(* ::Chapter:: *)
(*Export*)


g0se2data = Table[{t, g0se2rainb[t]}, {t, taus}];


Export["data/mat_2nd_g0se", g0se2data, "Table"];



