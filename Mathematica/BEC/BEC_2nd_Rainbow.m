(* ::Package:: *)

(* ::Title:: *)
(*Second Order in \[Tau] Rainbow *)


(* ::Chapter:: *)
(*\[CapitalSigma] *)


(* ::Section:: *)
(*Import*)


SetDirectory["/home/h/H.Guertner/theorie/BEC_Polaron/Mathematica/"];
<<params.mx;
ResetDirectory[];
p = Import["DiagMC_BEC.json", {"Data", "Momentum"} ];
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha] = Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ];
taus = Import["data/zero_order", "table"][[100;;200,1]];


(* ::Section:: *)
(*\[Tau] Integral*)


f[q1_,q2_,Theta_,\[Tau]_,\[Tau]1_, \[Tau]2_]:=(8*Pi^2)/(2*Pi)^6*q1^2*q2^2*vq2[q1]*vq2[q2]*g0[q1,0,\[Tau]1]*g0cos[q1,q2,\[Tau]1,\[Tau]2,Theta]*g0[q1,\[Tau]2,\[Tau]]*dq[q1,0,\[Tau]]*dq[q2,\[Tau]1,\[Tau]2];
g[q1_,q2_,Theta_,\[Tau]_,\[Tau]2_]:=Integrate[f[q1,q2,Theta,\[Tau],\[Tau]1,\[Tau]2],  {\[Tau]1,0,\[Tau]2}];
h[q1_,q2_,Theta_,\[Tau]_]:=Integrate[g[q1,q2,Theta,\[Tau],\[Tau]2],{\[Tau]2,0,\[Tau]}];


(*Plot[h[q,1],{q, 0,qc}, PlotPoints \[Rule] 10, MaxRecursion->5, PlotRange -> Full]*)


(* ::Section:: *)
(*q Integral*)


se[\[Tau]_]:=NIntegrate[h[q1,q2,Theta, \[Tau]],{q1,0, qc}, {q2,0,qc},{Theta,-1,1}];
(* Plot[\[CapitalSigma]1[t],{t, 0,1}, PlotPoints \[Rule] 10, MaxRecursion->5, PlotRange -> Full]*)
sedata = Table[{t, se[t]}, {t, taus}]


(* ::Chapter:: *)
(*Export*)


Export["data/mat_2nd_rainbow_se", sedata, "Table"];
