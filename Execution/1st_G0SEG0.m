(* ::Package:: *)

(* ::Title:: *)
(*First Order G0SEG0 in \[Tau] p=0*)


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


f[q_,\[Tau]_,\[Tau]1_, \[Tau]2_]:=(4*Pi)/(2*Pi)^3 *q^2*If[FP,vq2FP[q]*g0FP[0,0,\[Tau]1]*g0FP[-q,\[Tau]1,\[Tau]2]*dqFP[\[Tau]1,\[Tau]2]*g0FP[0,\[Tau]2,\[Tau]], vq2[q]*g0[0,0,\[Tau]1]*g0[-q,\[Tau]1,\[Tau]2]*dq[q,\[Tau]1,\[Tau]2]*g0[0,\[Tau]2,\[Tau]]];
g[q_,\[Tau]_,\[Tau]2_]:=Integrate[f[q,\[Tau],\[Tau]1,\[Tau]2],  {\[Tau]1,0,\[Tau]2}];
h[q_,\[Tau]_]:=Integrate[g[q,\[Tau],\[Tau]2],{\[Tau]2,0,\[Tau]}];


(* ::Chapter:: *)
(*q Integral*)


g0seg0[\[Tau]_]:=NIntegrate[h[q, \[Tau]],{q,0, qc}, AccuracyGoal->100];
(* Plot[\[CapitalSigma]1[t],{t, 0,1}, PlotPoints \[Rule] 10, MaxRecursion->5, PlotRange -> Full]*)
g0seg0data = Table[{t, g0seg0[t]}, {t, taus}];


(* ::Chapter:: *)
(*Export*)


Export["data/mat_1st_g0seg0", g0seg0data, "Table"];
