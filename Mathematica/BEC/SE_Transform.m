(* ::Package:: *)

(* ::Title:: *)
(*Transformation of SE to G using Matsubara for p=0*)


(* ::Chapter:: *)
(*Import*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<< params.mx;
ResetDirectory[];
p = 0;
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
\[Alpha]= Import["DiagMC_BEC.json", {"Data", "Alpha"} ];
mr = Import["DiagMC_BEC.json", {"Data", "Impurity_Mass"} ];
qc = Import["DiagMC_BEC.json", {"Data", "Q_Cutoff"} ];
Bins = Import["DiagMC_BEC.json", {"Data", "Bins"}];
FP = Import["DiagMC_BEC.json", {"Data", "Froehlich_Polaron"} ];
data=Import["data/SE/SE_>1", "Table"];
lastbin = If[SequenceCount[data[[All,2]], {0}] == Bins[[-1]], Bins[[-1]],  If[SequenceCount[data[[All,2]], {0}] == 0, Length[data]-10, SequencePosition[data[[All,2]],{0}][[1,1]]]];
data1= Import["data/SE/mat_1st_se", "Table"];
data[[All, 2]] += data1[[All, 2]];
tmax = data[[lastbin,1]];
wintmax = 1/(data[[3,1]]-data[[2,1]]);



(* ::Chapter:: *)
(*SE1 in iw analytical*)


f2[q_,w_, w1_]:=-(4*Pi)/(2*Pi)^4 *q^2*If[FP,vq2FP[q]*g0pwFPmu[-q,I w - I w1,0]*dqwFP[I w1,0], vq2[q]*g0pwmu[-q,I w - I w1,0]*dqw[q,I w1,0]];
g2[q_,w_] := Integrate[f2[q,w,w1], {w1,-Infinity, Infinity}]
h2[q_,w_]:=2* Pi* I*Residue[f2[q,w, w1], {w1, w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]}];
(*h2[q,w]
f2[q,w,w1]
f2[q,w,w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]]
Integrate[f2[q,w,w1], {w1,-Infinity, Infinity}]
Integrate[g2[q,w]/2/Pi*Exp[-I w t], {w,-Infinity, Infinity}]
2* Pi* I*Residue[f2[q,w,w1], {w1, w+ If[FP,I q^2/2, I q^2/2/mr*Sqrt[2]] - I \[Mu]}]
2* Pi *I* Residue[h2[q,w], {w, - I q^2/2 + I \[Mu] - I wq[q]}]*)


(*Plot[Re[h2[1,w]/2/Pi], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->3]
Plot[Im[h2[1,w]/2/Pi], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->3]*)


sewana[w_?NumericQ] := NIntegrate[h2[q,w], {q,0,qc}, AccuracyGoal->10];
(*sewana[0.01]*)


wc=100000;
se[t_]:=2*Re[NIntegrate[h2[q,w]/2/Pi*Exp[-I w *t],{q,0, qc},{w,0, wc},AccuracyGoal->10]];
(*se[0.01]*)
(*NIntegrate[h[q,w]/2/Pi*Exp[-I w *0.1],{q,0, qc},{w,-wc, wc},AccuracyGoal->1000]*)
(*Show[Plot[-se[t],{t,0,5},PlotPoints->10,MaxRecursion ->2], ListPlot[data1[[All,1;;2]]]]*)


(* ::Chapter:: *)
(*All Orders*)


(* ::Section:: *)
(*Interpolation*)


interpolse = Interpolation[data[[All,1;;2]]];
interpolsew[w_?NumericQ]:= -NIntegrate[interpolse[t]*Exp[I*w*t], {t, 0, tmax}, AccuracyGoal->10];
(*interpolsew[0.01]*)


(* ::Section:: *)
(*Simpson*)


h= (data[[3, 1]]-data[[2, 1]])/3;
lpf[itau_, w_] := data[[itau, 2]]*Exp[I*w*data[[itau, 1]]];
sew[\[Omega]_] := -h(lpf[1, \[Omega]]+lpf[lastbin, \[Omega]]+ Sum[4*lpf[2*itau, \[Omega]]+2*lpf[2*itau-1, \[Omega]], {itau, 1, lastbin/2}]); 
(*sew[0.01]
Plot[{Re[sew[w]], Re[interpolsew[w]]}, {w,0,wintmax}, PlotPoints->10, MaxRecursion\[Rule] 1]*)


(* ::Section:: *)
(*Comparison to analytical Result*)


(*Plot[{Re[sew[w]], Re[sewana[w]]}, {w, -10, 10}, PlotPoints->10, MaxRecursion->1]*)


(* ::Chapter:: *)
(*G-G0 in w*)


interpolgminusg0w[w_?NumericQ]:= (1/If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]] - interpolsew[w])^(-1) - If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]];
gminusg0w[w_]:= (1/If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]] - sew[w])^(-1) - If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]];
gminusg0wana[w_]:= (1/If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]] - sewana[w])^(-1) - If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]];
test[w_] := (1/If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]] - 1 )^(-1) - If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]];


(*Plot[{Re[gminusg0w[w]], Re[gminusg0wana[w]]}, {w,-wintmax+10,wintmax-10}, PlotPoints->10, MaxRecursion->2, PlotRange->Full]
Plot[Im[gminusg0w[w]], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->5, PlotRange->Full]*)


(*LogPlot[{Re[test[w]],Re[interpolgminusg0w[w]],Re[gminusg0w[w]], Re[gminusg0wana[w]]}, {w,0,1000}, PlotPoints->10, MaxRecursion->2, PlotRange->Full]*)


(* ::Chapter:: *)
(*ReFT of G-G0*)


interpolgminusg0[t_?NumericQ]:= NIntegrate[1/2/Pi*interpolgminusg0w[w]*Exp[-(I *w*t)], {w, -wintmax,wintmax}, AccuracyGoal->10];
gminusg0[t_?NumericQ]:= 2*Re[NIntegrate[1/2/Pi*gminusg0w[w]*Exp[-(I *w*t)], {w, 0,wintmax}, AccuracyGoal->10]];
gminusg0ana[t_?NumericQ]:= 2*Re[NIntegrate[1/2/Pi*gminusg0wana[w]*Exp[-(I *w*t)], {w, 0,wintmax}, AccuracyGoal->10]];
(*gminusg0[1]
interpolgminusg0[1]*)


gminusg0prove[t_?NumericQ]:= NIntegrate[1/2/Pi*gminusg0w[w]*Exp[-(I *w*t)], {w, -wintmax,wintmax}, AccuracyGoal->10];
gminusg0proveana[t_?NumericQ]:= NIntegrate[1/2/Pi*gminusg0wana[w]*Exp[-(I *w*t)], {w, -wintmax,wintmax}, AccuracyGoal-> 10];
(*gminusg0prove[1]
gminusg0proveana[1]
Plot[gminusg0prove[t], {t,0,3}, MaxRecursion->1, PlotPoints->10]*)


(* ::Chapter:: *)
(*complete G*)


interpolghans[t_]:=Abs[interpolgminusg0[t]]+If[FP,g0FP[p,0,t],g0[p,0,t]];
ghans[t_]:= -gminusg0[t]+If[FP,g0FP[p,0,t],g0[p,0,t]];
ghansana[t_]:= -gminusg0ana[t]+If[FP,g0FP[p,0,t],g0[p,0,t]];
(*ghans[1]
ghansana[1]
Plot[ghans[t], {t,0,10}, PlotPoints->10, MaxRecursion->3, PlotRange->Full]*)


(* ::Chapter:: *)
(*Export*)


(*Export Green function of second order*)
(*output = Table[{t, ghansana[t]}, {t, 0.1, 10, 0.1}];
interpoloutput = Table[{data[[itau,1]], interpolghans[data[[itau,1]]]}, {itau, 1, Length[data]}]
Export["data/Trans_1st_se_all_orders", output, "Table"];
Export["data/FP_control/mat_FP_all_transform_interpol", output, "Table"];*)


(*Export Green function of second order*)
output2 = Table[{t, ghans[t]}, {t, 0.1, 10, 0.1}];
(*interpoloutput = Table[{data[[itau,1]], interpolghans[data[[itau,1]]]}, {itau, 1, Length[data]}]*)
Export["data/Trans_se_all_orders", output2, "Table"];
(*Export["data/FP_control/mat_FP_all_transform_interpol", output, "Table"];*)
