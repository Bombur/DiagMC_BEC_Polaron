(* ::Package:: *)

(* ::Title:: *)
(*First Order SE in \[Tau] p=0*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<<params.mx;
ResetDirectory[];
p = 0;
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha] = Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ];
taus = Import["data/zero_order", "table"][[All,1]];
FP= Import["DiagMC_BEC.json", {"Data", "Froehlich_Polaron"} ];




(* ::Chapter:: *)
(*SE in iw analytical*)


f2[q_,w_, w1_]:=-(4*Pi)/(2*Pi)^4 *q^2*If[FP,vq2FP[q]*g0pwFPmu[-q,I w - I w1,0]*dqwFP[I w1,0], vq2[q]*g0pwmu[-q,I w - I w1,0]*dqw[q,I w1,0]];
g2[q_,w_] := Integrate[f2[q,w,w1], {w1,-Infinity, Infinity}]
h2[q_,w_]:=2* Pi* I*Residue[f2[q,w, w1], {w1, w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]}];
h2[q,w];
(*f2[q,w,w1]
Simplify[f2[q,w,w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]]]
Integrate[f2[q,w,w1], {w1,-Infinity, Infinity}]
(*Integrate[g2[q,w]/2/Pi*Exp[-I w t], {w,-Infinity, Infinity}]*)
2* Pi* I*Residue[f2[q,w,w1], {w1, w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]}]
2* Pi *I* Residue[h2[q,w], {w, - I q^2/2 + I \[Mu] - I wq[q]}]*)


(*Plot[Re[h2[1,w]/2/Pi], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->3]
Plot[Im[h2[1,w]/2/Pi], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->3]*)


wc=10^7;
sew[t_]:=2*Re[NIntegrate[h2[q,w]/2/Pi*Exp[-I w *t],{q,0, qc},{w,0, wc},AccuracyGoal->5]];


(*sewtest[t_, wc_]:=2*Re[NIntegrate[h2[q,w]/2/Pi*Exp[-I w *t],{q,0, qc},{w,0, wc},AccuracyGoal->5]];
LogLinearPlot[sewtest[1,wc], {wc,10,100000000}, MaxRecursion->2, PlotPoints->10]*)


(* ::Chapter:: *)
(*t Integral*)


f[q_,\[Tau]_]:=(4*Pi)/(2*Pi)^3 *q^2* If[FP,vq2FP[q]*g0FP[-q,0,\[Tau]]*dqFP[0,\[Tau]],vq2[q]*g0[-q,0,\[Tau]]*dq[q,0,\[Tau]]];


(*g0mu[pabs_, t1_, t2_,\[Mu]_ ]:= Exp[-(pabs^2/(2*mr)*Sqrt[2]-\[Mu])(t2-t1)];
test[q_,t_,mu_]:= (4*Pi)/(2*Pi)^3 *q^2*vq2[q]*g0mu[-q,0,t,mu]*dq[q,0,t];
setest[t_, mu_]:= NIntegrate[test[q, t, mu],{q,0, qc},AccuracyGoal->1000];
data= Import["data/SE/mat_1st_se_mu-1", "Table"];
data2= Import["data/SE/mat_1st_se_mu-2", "Table"];
data3= Import["data/SE/mat_1st_se_mu-2.5", "Table"];
ListLogPlot[{data, data2,data3}]
Show[LogPlot[{setest[t, -1], setest[t, -2],setest[t, -2.5]}, {t, 0.1, 6}, MaxRecursion\[Rule]1], ListLogPlot[{data, data2,data3}]]*)



(* ::Chapter:: *)
(*q Integral*)


se[\[Tau]_]:=NIntegrate[f[q, \[Tau]],{q,0, qc},AccuracyGoal->1000];


(* ::Section:: *)
(*Comparison to Analytical Result*)


(*se[1]
sew[1]
Plot[{-sew[t], se[t]}, {t,0,3}, MaxRecursion->1, PlotPoints->10]*)


(* ::Chapter:: *)
(*Export*)


sedata = Table[{t, se[t]}, {t, taus}];


Export["data/SE/mat_1st_se", sedata, "Table"];
