(* ::Package:: *)

(* ::Title:: *)
(*Transformation of G0SE to G using Matsubara for p=0*)


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
fog0seg0 = Import["DiagMC_BEC.json", {"Data", "FOG0SEG0"} ];
data=Import["data/all_orders", "Table"];
lastbin = If[SequenceCount[data[[All,2]], {0}] == Bins[[-1]], Bins[[-1]],  If[SequenceCount[data[[All,2]], {0}] == 0, Length[data]-10, SequencePosition[data[[All,2]],{0}][[1,1]]]]
data1= Import["data/mat_1st_g0se", "Table"];
data[[All, 2]] += data1[[All, 2]];
tmax = data[[lastbin,1]];
wintmax = 1/(data[[3,1]]-data[[2,1]]);
(*data3=Import["data/all_orders", "Table"];
ListPlot[{data[[All,2]],data1[[All,2]],data3[[All,2]],data1[[All,2]]+data3[[All,2]] }, PlotRange\[Rule]All]*)


(* ::Chapter:: *)
(*All Orders*)


(* ::Section:: *)
(*Interpolation*)


interpolg0se = Interpolation[data[[All,1;;2]]];
interpolg0sew[w_?NumericQ]:= -NIntegrate[interpolg0se[t]*Exp[I*w*t], {t, 0, tmax}];
(*interpolg0sew[0.01]*)


(* ::Section:: *)
(*Simpson*)


h= (data[[3, 1]]-data[[2, 1]])/3;
lpf[itau_, \[Omega]_] := data[[itau, 2]]*Exp[I*\[Omega]*data[[itau, 1]]];
g0sew[\[Omega]_] := -h(lpf[1, \[Omega]]+lpf[lastbin, \[Omega]]+ Sum[4*lpf[2*itau, \[Omega]]+2*lpf[2*itau-1, \[Omega]], {itau, 1, lastbin/2}]); 
(*g0sew[0.01]*)


(* ::Chapter:: *)
(*G-G0 in w*)


interpolgminusg0w[w_?NumericQ]:= If[FP, g0pwFPmu[p, I w, 0], g0pwmu[p, I w, 0]]*(interpolg0sew[w]/(1-interpolg0sew[w]));
gminusg0w[\[Omega]_]:= If[FP, g0pwFPmu[p, I \[Omega], 0], g0pwmu[p, I \[Omega], 0]]*(g0sew[\[Omega]]/(1-g0sew[\[Omega]]));
(*Plot[{Re[gminusg0w[w]], Re[interpolgminusg0w[w]]}, {w,-wintmax+10,wintmax-10}, PlotPoints->10, MaxRecursion->3, PlotRange->Full]*)


(*Plot[Im[gminusg0w[w]], {w,-wintmax,wintmax}, PlotPoints->10, MaxRecursion->5, PlotRange->Full]*)


(* ::Chapter:: *)
(*ReFT of G-G0*)


interpolgminusg0[t_?NumericQ]:= NIntegrate[1/2/Pi*interpolgminusg0w[w]*Exp[-(I *w*t)], {w, -wintmax,wintmax}];
gminusg0[t_?NumericQ]:= 2*Re[NIntegrate[1/2/Pi*gminusg0w[w]*Exp[-(I *w*t)], {w, 0,wintmax}]];
gminusg0[0.01]


gminusg0prove[t_?NumericQ]:= NIntegrate[1/2/Pi*gminusg0w[w]*Exp[-(I *w*t)], {w, -wintmax,wintmax}];
gminusg0prove[0.01]


(* ::Chapter:: *)
(*complete G*)


interpolghans[t_]:=Abs[interpolgminusg0[t]]+If[FP,g0FP[p,0,t],g0[p,0,t]];
ghans[t_]:= -gminusg0[t]+If[FP,g0FP[p,0,t],g0[p,0,t]];
(*ghans[0.01]
Plot[ghans[t], {t,0,10}, PlotPoints->10, MaxRecursion->3, PlotRange->Full]*)


(* ::Chapter:: *)
(*Export*)


(*Export Green function of second order*)
output = Table[{t, ghans[t]}, {t, 0.1, tmax, 0.1}]
(*interpoloutput = Table[{data[[itau,1]], interpolghans[data[[itau,1]]]}, {itau, 1, Length[data]}]*)
Export["data/Trans_g0se_all_orders", output, "Table"];
(*Export["data/FP_control/mat_FP_all_transform_interpol", output, "Table"];*)



