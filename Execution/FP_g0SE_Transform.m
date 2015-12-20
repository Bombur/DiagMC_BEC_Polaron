(* ::Package:: *)

(* ::Title:: *)
(*FP Transformation of G0SE*)


(* ::Chapter:: *)
(*Import*)


SetDirectory["/project/theorie/h/H.Guertner/BEC_Polaron/Mathematica/"];
<< params.mx;
ResetDirectory[];
p =Import["DiagMC_BEC.json", {"Data", "Momentum"} ];
\[Mu] = Import["DiagMC_BEC.json", {"Data", "Chemical_Potential"} ];
peter=Import["/home/h/H.Guertner/theorie/BEC_Polaron/peter/FP_all", "Table"];
peter = Insert [peter, {0,1,0}, 1];
peter1=Import["/home/h/H.Guertner/theorie/BEC_Polaron/peter/FP_1", "Table"];
peter1 = Insert [peter1, {0,0,0}, 1];
data=Import["data/all_orders", "Table"];
data = Insert [data, {0,1,0}, 1];
data1= Import["data/first_order", "Table"];
data1 = Insert [data1, {0,0,0}, 1];


(* ::Chapter:: *)
(*First Order*)


(* ::Section:: *)
(*Interpolation*)


interpolg0se1 = Interpolation[data1[[All,1;;2]]];
interpolg0se1w[w_]:= NIntegrate[interpolg0se1[t]*Exp[I*w*t], {t, 0, 5}];


(* ::Section:: *)
(*Simpson*)


h1= (data1[[3, 1]]-data1[[2, 1]])/3;
lpf1[itau_, \[Omega]_] := data1[[itau, 2]]*Exp[I*\[Omega]*data1[[itau, 1]]];
g0se1w[\[Omega]_] := h1*(lpf1[1, \[Omega]]+lpf1[-1, \[Omega]]+ Sum[4*lpf1[itau, \[Omega]], {itau, 2, Length[data1]-1, 2}] + Sum[2*lpf1[itau, \[Omega]], {itau, 3, Length[data1]-2, 2}]); 


(* ::Chapter:: *)
(*Orders >= 2*)


(* ::Section:: *)
(*Interpolation*)


interpolg0sebig2 = Interpolation[data[[All,1;;2]]];
interpolg0sebig2w[w_]:= NIntegrate[interpolg0sebig2[t]*Exp[I*w*t], {t, 0, 4.99}];
interpolg0sew[w_] :=interpolg0sebig2w[w]+ interpolg0se1w[w]; 


(* ::Section:: *)
(*Simpson*)


h= (data[[3, 1]]-data[[2, 1]])/3;
lpf[itau_, \[Omega]_] := data[[itau, 2]]*Exp[I*\[Omega]*data[[itau, 1]]];
g0sebig2w[\[Omega]_] := h(lpf[1, \[Omega]]+lpf[-1, \[Omega]]+ Sum[4*lpf[itau, \[Omega]], {itau, 2, Length[data]-1, 2}] + Sum[2*lpf[itau, \[Omega]], {itau, 3, Length[data]-2, 2}]); 
g0sew[\[Omega]_] := g0se1w[\[Omega]]+g0sebig2w[\[Omega]];


(* ::Chapter:: *)
(*G-G0 in w*)


interpolgminusg0w[w_]:= g0pw[p, I w, 0]*(interpolg0sew[w]/(1-interpolg0sew[w]));
gminusg0w[\[Omega]_]:= g0pw[p, I \[Omega], 0]*(g0sew[\[Omega]]/(1-g0sew[\[Omega]]));


(* ::Chapter:: *)
(*ReFT of G-G0*)


interpolgminusg0[t_]:= NIntegrate[1/2/Pi*interpolgminusg0w[w]*Exp[-(I *w*t)], {w, -10,10}];
gminusg0[t_]:= NIntegrate[1/2/Pi*gminusg0w[w]*Exp[-(I *w*t)], {w, -10,10}];


(* ::Chapter:: *)
(*complete G*)


interpolghans[t_]:=Abs[interpolgminusg0[t]]+g0[p,0,t];
ghans[t_]:= Abs[gminusg0[t]]+g0[p,0,t];


(* ::Chapter:: *)
(*Export*)


(*Export Green function of second order*)
output = Table[{data[[itau,1]], ghans[data[[itau,1]]]}, {itau, 1, Length[data]}]
(*interpoloutput = Table[{data[[itau,1]], interpolghans[data[[itau,1]]]}, {itau, 1, Length[data]}]*)
Export["data/FP_control/mat_FP_all_transform", output, "Table"];
(*Export["data/FP_control/mat_FP_all_transform_g0se_interpol", output, "Table"];*)
