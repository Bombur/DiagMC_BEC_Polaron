(* ::Package:: *)

(* ::Title:: *)
(*2nd Order Crossed G0SE Diagram in \[Tau] p=0*)


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
(*Other*)


Integ2FP[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_, \[Tau]3_,  ctheta_] := vq2FP[q1]*vq2FP[q2]*g0FP[0,0,\[Tau]1]*g0FP[q1,\[Tau]1, \[Tau]2]*dqFP[\[Tau]1,\[Tau]3]*g0cosFP[q1,q2,\[Tau]2,\[Tau]3, ctheta]*dqFP[\[Tau]2,\[Tau]]*g0FP[q2,\[Tau]3,\[Tau]];
Integ2BEC[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_, \[Tau]3_, ctheta_] := vq2[q1]*vq2[q2]*g0[0,0,\[Tau]1]*g0[q1,\[Tau]1,\[Tau]2]*dq[q1,\[Tau]1,\[Tau]3]*g0cos[q1,q2,\[Tau]2,\[Tau]3, ctheta]*dq[q2,\[Tau]2,\[Tau]]*g0[q2,\[Tau]3,\[Tau]];
f2[q1_,q2_,\[Tau]_,\[Tau]1_, \[Tau]2_, \[Tau]3_, ctheta_]:=4*Pi*2*Pi/(2*Pi)^6 *q1^2*q2^2*If[FP,Integ2FP[q1,q2,\[Tau],\[Tau]1, \[Tau]2, \[Tau]3, ctheta], Integ2BEC[q1,q2,\[Tau],\[Tau]1, \[Tau]2, \[Tau]3, ctheta]];
g2[q1_,q2_, \[Tau]_,\[Tau]2_,\[Tau]3_, ctheta_] := Integrate[f2[q1,q2,\[Tau],\[Tau]1,\[Tau]2,\[Tau]3,ctheta],  {\[Tau]1,0,\[Tau]2}];
h2[q1_,q2_, \[Tau]_,\[Tau]3_, ctheta_] := Integrate[g2[q1,q2,\[Tau],\[Tau]2,\[Tau]3,ctheta],  {\[Tau]2,0,\[Tau]3}];
i2[q1_,q2_, \[Tau]_, ctheta_] := Integrate[h2[q1,q2,\[Tau],\[Tau]3, ctheta],  {\[Tau]3,0,\[Tau]}];


g0se2oth[\[Tau]_?NumericQ]:=NIntegrate[i2[q1, q2, \[Tau], ctheta],{q1,0, qc}, {q2,0, qc}, {ctheta, -1,1}, AccuracyGoal->100];
(*LogPlot[se2oth[t],{t, 0,10}, PlotPoints -> 10, MaxRecursion->3, PlotRange -> Full]*)


(* ::Chapter:: *)
(*Export*)


g0se2data = Table[{t, g0se2oth[t]}, {t, taus}]


Export["data/mat_2nd_g0se_crossed", g0se2data, "Table"];



