(* ::Package:: *)

(* ::Title:: *)
(*Polaron Energy in 1.Order p=0*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<<params.mx;
ResetDirectory[];
p = 0;
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha]= Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ];
wlist = Import["DiagMC_BEC.json", {"Data", "Ws_for_Epol"} ];
FP = Import["DiagMC_BEC.json", {"Data", "Froehlich_Polaron"} ];


(* ::Chapter:: *)
(*Integrand for q Integral*)


g[q_,\[Omega]_, \[CurlyEpsilon]_ ]:=(4*Pi)/(2*Pi)^3 *q^2* If[FP,vq2FP[q]/(\[Omega]-q^2/2-1+I*\[CurlyEpsilon]), vq2[q]/(\[Omega]-q^2/2/mr*Sqrt[2]-wq[q]+I*\[CurlyEpsilon])];
(*Plot[-g[q,\[Mu] + 1,0],{q,0,qc}]*)




(* ::Chapter:: *)
(*\[CapitalSigma]1 in omega*)


\[CapitalSigma]1[\[Omega]_]:=NIntegrate[g[q, \[Omega], 0],{q,0, qc}, AccuracyGoal->100];


wgerade[w_] := w;
Epol= w /. FindRoot[\[CapitalSigma]1[w]-w, {w,0}];
Epol2 =Epol + eRenqc[\[Alpha],mr,qc];


(* ::Chapter:: *)
(*\[CapitalSigma]1 Transform from tau*)


f2[q_,\[Tau]_,\[Omega]_]:=(4*Pi)/(2*Pi)^3 *q^2* Exp[(\[Omega]- \[Mu])*\[Tau]]*If[FP, vq2FP[q]*g0FP[-q,0,\[Tau]]*dqFP[0,\[Tau]], vq2[q]*g0[-q,0,\[Tau]]*dq[q,0,\[Tau]]];
(*Plot[f2[q,1, \[Mu] + 2],{q,0,2}]*)
g2[q_,\[Omega]_]= Integrate[f2[q,t,\[Omega]], {t,0,10}];
(*Plot[g2[q, \[Mu] + 1],{q,0,qc}]*)
\[CapitalSigma]1trans[\[Omega]_]:= NIntegrate[g2[q, \[Omega]],{q,0, qc}, AccuracyGoal->1000];
Epol2= w /. FindRoot[\[CapitalSigma]1trans[-w]-w, {w,0}]

Epplot = Show[Plot[{-\[CapitalSigma]1[-w], \[CapitalSigma]1trans[-w], wgerade[w]},{w,-wlist[[-1]],-wlist[[1]]},PlotPoints->10,MaxRecursion ->3,PlotRange->{-Epol-1,-Epol+1},AxesLabel->{"-\[Omega]", "-\!\(\*SubscriptBox[\(E\), \(p\)]\)[c/\[Xi]]"}], ListPlot[{{-Epol, -Epol}}, PlotStyle-> {Red}, PlotLegends-> Placed[Epol2, {Right, Center}]]];


(* ::Chapter:: *)
(*Export*)


Epdata= Table[{w,-\[CapitalSigma]1[w]}, {w, wlist} ];
Export["data/Ep/mat_Ep1st", Epdata, "Table"];
Export["ana/Ep/1st_Epvsws.pdf", Epplot];



