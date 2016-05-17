(* ::Package:: *)

(* ::Title:: *)
(*G0SE in Matsubara in 1.Order p=0*)


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
wbin = Import["DiagMC_BEC.json", {"Data", "Omega_bin"} ];
wmax = Import["DiagMC_BEC.json", {"Data", "Omega_max"} ];


(* ::Chapter:: *)
(*Integrand for q Integral*)


f[q_,\[Omega]_]:=(4*Pi)/(2*Pi)^3 *q^2* If[FP,vq2FP[q]/(I \[Omega]-q^2/2-1+\[Mu])/(I \[Omega] +\[Mu]) , vq2[q]/(I \[Omega]-q^2/2/mr*Sqrt[2]-wq[q]+\[Mu])/(I \[Omega] +\[Mu])];
(*Plot[-g[q,\[Mu] + 1,0],{q,0,qc}]*)




(* ::Chapter:: *)
(*G0SE1 in omega*)


g0seiw[\[Omega]_]:=NIntegrate[f[q, \[Omega]],{q,0, qc}, AccuracyGoal->100];


g0seiwdata = Table[{N[w],g0seiw[w]}, {w, 0,wmax-wmax/wbin, wmax/wbin} ]


(* ::Chapter:: *)
(*Export*)


writeComplex[file_String, data_List] := 
 Module[{str = OpenWrite[file, FormatType -> OutputForm]}, 
  Scan[Write[str, #] &, data /. {Complex[x_, y_] :> ToString[x] <> "," <> ToString[y] <> "i"}]; 
  Close[str];
]
writeComplex["data/G0SEiw/mat_1st_g0seiw", g0seiwdata[[All,2]]]
