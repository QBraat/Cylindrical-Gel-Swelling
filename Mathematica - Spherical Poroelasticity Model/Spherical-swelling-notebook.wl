(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]


(* ::Text:: *)
(*Keeping track of time for large computations*)


(* ::Input::Initialization:: *)
start2 = AbsoluteTime[];


(* ::Text:: *)
(*Set the total number of iterations (titer) and indicate when you want a notification each 'check' number of iterations*)


(* ::Input::Initialization:: *)
titer =15000;
check =1000;


(* ::Text:: *)
(*Allow for continuation after termination. The data must be in the same folder as the notebook itself in order for the notebook to continue on the cluster. *)


(* ::Input::Initialization:: *)
NEW = "Y";  "N" 


(* ::Subsection::Initialization::Closed:: *)
(*(*Parameters*)*)


(* ::Text:: *)
(*Model parameters*)


(* ::Input::Initialization:: *)
Ns =15; 
\[CapitalOmega]=1.28*10^-4;
\[Chi]=.4;
\[Alpha]=250.;
\[Beta]=3/2;


(* ::Text:: *)
(*Environment*)


(* ::Input::Initialization:: *)
\[Mu]fs =-2.*10^3;
\[Mu]f0 =-5*10^3;


(* ::Text:: *)
(*Time discretization*)


(* ::Input::Initialization:: *)
dt=5.*10^-8;
t0=0; 


(* ::Subsection::Initialization:: *)
(*(*Functions*)*)


(* ::Subsubsection::Initialization:: *)
(*(*Functions*)*)


(* ::Text:: *)
(*Functions for all quantities such as stresses etcetera*)


(* ::Input::Initialization:: *)
ur[r_]:=r-(r^3-3Integrate[\[Rho]^2 \[Phi]f[\[Rho]],{\[Rho],0,r}])^(1/3);
\[Lambda]r[r_]:=1/(r^2 ((1-\[Phi]f[r])/(r-ur[r])^2));
\[Lambda]\[Theta][r_]:=1/(1-ur[r]/r);
J[r_]:=1/(1-\[Phi]f[r]);
\[Sigma]rp[r_]:=(\[Lambda]r[r]^2-1)/J[r];
\[Sigma]\[Theta]p[r_]:=(\[Lambda]\[Theta][r]^2-1)/J[r];
\[CapitalPi][r_]:= -(1/\[CapitalOmega])*(1/J[r]+Log[1-1/J[r]]-1/(\[Alpha] J[r])+\[Chi]/J[r]^2);


(* ::Subsubsection::Initialization:: *)
(*(*Import previous data - continue*)*)


(* ::Text:: *)
(*Import old data and continue to write data to these datasets*)


(* ::Input::Initialization:: *)
ImporterenOudeData[]:=Module[{},
dtbefore = Import[NotebookDirectory[]<>"dt-before.mat"][[1]];
dtupdate =  Import[NotebookDirectory[]<>"dt-update.mat"][[1]];

\[Mu]import = Import[NotebookDirectory[]<>"data\[Mu]-excludingdt.mat"];
\[Phi]fupdateimport = Import[NotebookDirectory[]<>"\[Phi]fupdate-excludingdt.mat"];
\[Sigma]rleftimport =Import[NotebookDirectory[]<>"\[Sigma]rleft-excludingdt.mat"];
\[Sigma]\[Theta]leftimport =Import[NotebookDirectory[]<>"\[Sigma]\[Theta]left-excludingdt.mat"];
\[CapitalPi]leftimport =Import[NotebookDirectory[]<>"\[CapitalPi]left-excludingdt.mat"];
urimport = Import[NotebookDirectory[]<>"ur-excludingdt.mat"];
a0import = Import[NotebookDirectory[]<>"a0.mat"];

data\[Mu] = Transpose[Join[{Transpose[dtbefore][[1]],\[Mu]import}]];
data\[Phi]farray = Transpose[Join[{Transpose[dtupdate][[1]],\[Phi]fupdateimport}]];
data\[Sigma]r = Transpose[Join[{Transpose[dtbefore][[1]],\[Sigma]rleftimport}]];
data\[Sigma]\[Theta] = Transpose[Join[{Transpose[dtbefore][[1]],\[Sigma]\[Theta]leftimport}]];
data\[CapitalPi]=Transpose[Join[{Transpose[dtbefore][[1]],\[CapitalPi]leftimport}]];
dataur = Transpose[Join[{Transpose[dtbefore][[1]],urimport}]];
dataa0 = a0import[[1]];

datad\[Mu]={};
data\[Delta]a={}; 

lengthN = Length[dtbefore];
]


(* ::Subsubsection::Initialization:: *)
(*(*Working modules*)*)


(* ::Text:: *)
(*Create empty lists for storage*)


(* ::Input::Initialization:: *)
EmptyArray[]:=Module[{},
data\[Sigma]r = {};
data\[Sigma]\[Theta] = {};
data\[CapitalPi]={};
data\[Mu]={};
data\[Phi]farray = {};
data\[Delta]a={};
dataa0 = {{0,a0}};
dataur={};
]


(* ::Text:: *)
(*Initialize the set of parameters *)


(* ::Input::Initialization:: *)
ParametersLoop[]:=Module[{},
If[NEW == "N" && ind == lengthN+1, anew = dataa0[[-1,2]]]; (*check if you need to continue with previous set or not*)
If[ind==1, a0=a0, a0=anew];(* initialize with \[Mu]f0 parameter is ind = 1, and continue with the a(t) when ind > 1*)

dr = a0/Ns; set gridsize
If[NEW=="N"&&ind==lengthN+1,dadtnow= (dataa0[[-1,2]]-dataa0[[-2,2]])*dt]; (*set estimate for da/dt when you continue from previous set*)

If[ind==1,dadtguess=4000., dadtguess=dadtnow]; (*set estimate for da/dt*)
]


(* ::Text:: *)
(*Calculation of \[Phi]f(t) and the corresponding deformation*)


(* ::Input::Initialization:: *)
Initialize\[Phi]fur[]:=Module[{},
If[NEW=="N" && ind == lengthN+1,\[Phi]farr = data\[Phi]farray[[-1,2]]];

If[ind == 1, \[Phi]f[x_]:=(1-\[Phi]feq)*HeavisideTheta[x-a0]+\[Phi]feq,
\[Phi]f[x_]:=(Sum[\[Phi]farr[[i,2]]*(HeavisideTheta[x-(i-1)*dr]-HeavisideTheta[x-i*dr]),{i,1,Ns}]+HeavisideTheta[x-a0])]; (*set function \[Phi]f (t). Since Subscript[\[Phi]f, i] is discrete, it is a combination of Heaviside (step) functions*)

urlist = Table[{n/Ns*a0,ur[n/Ns*a0]},{n,0,Ns}]; (*calculate the deformation of the polymeric network in the radial direction*)
]


(* ::Text:: *)
(*Calculate the effective stresses and osmotic pressure; for a given \[Phi]f(t), one can use the functions defined earlier *)


(* ::Input::Initialization:: *)
GenerateStressPressure[]:=Module[{},
rcenter = Table[Subscript[r, i]=(i-1/2)*dr,{i,1,Ns}];
list\[Sigma]r = Table[0,{i,1,Ns}];
list\[Sigma]\[Theta]= Table[0,{i,1,Ns}];
list\[CapitalPi] = Table[0,{i,1,Ns}];

For[i=1,i<= Ns,i++,{
list\[Sigma]r[[i]]={Subscript[r, i],\[Sigma]rp[Subscript[r, i]]},
list\[Sigma]\[Theta][[i]]={Subscript[r, i],\[Sigma]\[Theta]p[Subscript[r, i]]},
list\[CapitalPi][[i]]={Subscript[r, i],\[CapitalPi][Subscript[r, i]]}
}];
]


(* ::Text:: *)
(*Calculation of the pore pressure and the chemical potential. We have used a different approach than in Bertrand et al. It is outlined in the report, and we refer to the different steps here. *)


(* ::Input::Initialization:: *)
GenerateChemicalPotential[]:=Module[{},
dpdr = Table[{i*dr,(list\[Sigma]r[[i+1,2]]-list\[Sigma]r[[i,2]])/dr+2*1/(i*dr) (1/2 (list\[Sigma]r[[i+1,2]]+list\[Sigma]r[[i,2]])-1/2 (list\[Sigma]\[Theta][[i+1,2]]+list\[Sigma]\[Theta][[i,2]]))},{i,1,Ns-1}]; (*calculation of the dp/dr from the mechanical equilibrium definition*)

pa = list\[Sigma]r[[-1]]; (*set p(a) = \[Sigma]r'(a) due to the stress free condition*)

p0 = pa[[2]] - Sum[dpdr[[i,2]]*dr,{i,1,Ns-1}]; (*calculate p(0) via a repeated application of the Taylor expansion over the different boundaries*)

\[Mu]f0 = p0-list\[CapitalPi][[1,2]]; (*use definition of \[Mu] = p - Pi to calculate \[Mu](0)*)

d\[Mu] = Table[{i*dr,(list\[Sigma]r[[i+1,2]]-list\[Sigma]r[[i,2]])/dr+2*1/(i*dr) (1/2 (list\[Sigma]r[[i+1,2]]+list\[Sigma]r[[i,2]])-1/2 (list\[Sigma]\[Theta][[i+1,2]]+list\[Sigma]\[Theta][[i,2]]))-(list\[CapitalPi][[i+1,2]]-list\[CapitalPi][[i,2]])/dr},{i,1,Ns-1}]; (*calculate d\[Mu]/dr from mechanical equilibrium
*)
\[Mu]=Table[{(j+1/2)*dr,\[Mu]f0+Sum[d\[Mu][[i,2]]*dr,{i,1,j}]},{j,1,Ns-1}]; (*calculate \[Mu](r) via a repeated application of the Taylor expansion over the different boundaries*)
PrependTo[\[Mu],{1/2*dr,\[Mu]f0}]; (*prepend \[Mu]0*)
AppendTo[\[Mu],{dr*(Ns+1/2),\[Mu]fs}]; (*append environment*)

d\[Mu]edge = (\[Mu][[-1,2]]-\[Mu][[-2,2]])/dr; (*calculate gradient over the outer boundary by (\[Mu](r>a) - \[Mu](r<a))/dr *)
AppendTo[d\[Mu],{dr*Ns,d\[Mu]edge}];(* append to d\[Mu]/dr*)
]


(* ::Text:: *)
(*Update the quantities  \[Phi]f(t+dt) and a(t+dt) once they are determined (below)*)


(* ::Input::Initialization:: *)
Update\[Phi]f[]:=Module[{},
anew = a0+Re[dadt]*dt/.solution;
\[Delta]a = dadt*dt/.solution;
dadtnow = dadt/.solution;

drnew = anew/Ns;
\[Phi]farr = Table[{(i-1/2)*drnew,Re[Subscript[\[Phi]fnew, i]]},{i,1,Ns}]/.solution;

]


(* ::Text:: *)
(*Save data*)


(* ::Input::Initialization:: *)
SaveData[]:=Module[{},
AppendTo[dataur,{(ind-1)*dt,urlist}];
AppendTo[data\[Sigma]r,{(ind-1)*dt,list\[Sigma]r}];
AppendTo[data\[Sigma]\[Theta],{(ind-1)*dt, list\[Sigma]\[Theta]}];
AppendTo[data\[CapitalPi],{(ind-1)*dt, list\[CapitalPi]}];
AppendTo[data\[Mu],{(ind-1)*dt,\[Mu]}];
AppendTo[data\[Phi]farray,{ind*dt,\[Phi]farr}]; 
AppendTo[data\[Delta]a,{ind*dt,\[Delta]a}];
AppendTo[dataa0,{ind*dt,anew}];
]


(* ::Text:: *)
(*Setting up the set of discrete numerical equations, and solve these equations numerically *)


(* ::Input::Initialization:: *)
SolveSystem[]:=Module[{},
rarr=Table[Subscript[r, i]=(i-1/2)*dr,{i,1,Ns}];
rleft= Table[Subscript[r, i-1/2]=(i-1)*dr,{i,1,Ns}];
rright = Table[Subscript[r, i+1/2]=i*dr,{i,1,Ns}];
\[Phi]fold = Table[Subscript[\[Phi]f, i]=\[Phi]f[Subscript[r, i]],{i,1,Ns}];
Subscript[\[Phi]f, Ns+1]=1;

eqs= {4/3 \[Pi](Subscript[r, 3/2]^3-Subscript[r, 1/2]^3)*((Subscript[\[Phi]fnew, 1]-Subscript[\[Phi]f, 1])+((3Subscript[\[Phi]f, 1])/a0)dadt*dt)-4\[Pi]*(((Subscript[r, 3/2]^3*1/2*(Subscript[\[Phi]f, 2]+Subscript[\[Phi]f, 1]))/a0)dadt*dt+Subscript[r, 3/2]^2*((1/2 (Subscript[\[Phi]f, 2]+Subscript[\[Phi]f, 1]))/(1-1/2 (Subscript[\[Phi]f, 2]+Subscript[\[Phi]f, 1]))^(\[Beta]-1))d\[Mu][[1,2]]*dt)==0};
For[i=2,i<= Ns,i++, AppendTo[eqs,
4/3 \[Pi](Subscript[r, i+1/2]^3-Subscript[r, i-1/2]^3)*((Subscript[\[Phi]fnew, i]-Subscript[\[Phi]f, i])+((3Subscript[\[Phi]f, i])/a0)dadt*dt)-4\[Pi]*(((Subscript[r, i+1/2]^3*1/2*(Subscript[\[Phi]f, i+1]+Subscript[\[Phi]f, i]))/a0)dadt*dt+Subscript[r, i+1/2]^2*((1/2 (Subscript[\[Phi]f, i+1]+Subscript[\[Phi]f, i]))/(1-1/2 (Subscript[\[Phi]f, i+1]+Subscript[\[Phi]f, i]))^(\[Beta]-1))d\[Mu][[i,2]]*dt)+4\[Pi](((Subscript[r, i-1/2]^3*1/2*(Subscript[\[Phi]f, i]+Subscript[\[Phi]f, i-1]))/a0)dadt*dt+Subscript[r, i-1/2]^2*((1/2 (Subscript[\[Phi]f, i-1]+Subscript[\[Phi]f, i]))/(1-1/2 (Subscript[\[Phi]f, i-1]+Subscript[\[Phi]f, i]))^(\[Beta]-1))d\[Mu][[i-1,2]]*dt)==0]];

(*above, the N discretized conservation laws are listed*)

AppendTo[eqs, (a0+dadt*dt)^3*(1-Sum[Subscript[\[Phi]fnew, i]*((i/Ns)^3-((i-1)/Ns)^3),{i,1,Ns}])==1]; (*discrete material boundary condition*)
(*
unknowns including the initial guesses*)
list\[Phi]new = Table[{Subscript[\[Phi]fnew, i],\[Phi]f[Subscript[r, i]],0.,1.},{i,1,Ns}];
AppendTo[list\[Phi]new,{dadt,dadtguess}];
(*
solution to the discrete equations*)
solution =FindRoot[eqs,list\[Phi]new,WorkingPrecision->100,PrecisionGoal->100,MaxIterations->1000];
]


(* ::Text:: *)
(*Check the deformation boundary condition ur(a,t) = a(t)-1. Error when not satisfied (within a numerical error)*)


(* ::Input::Initialization:: *)
CheckCondition[]:=Module[{},
If[ind==1,
\[Phi]f[x_]:=(1-\[Phi]feq)*HeavisideTheta[x-a0]+\[Phi]feq,
\[Phi]f[x_]:=(Sum[\[Phi]farr[[i,2]]*(HeavisideTheta[x-(i-1)*drnew]-HeavisideTheta[x-i*drnew]),{i,1,Ns}]+HeavisideTheta[x-anew])
];

If[ind==1,anew =a0];
If[Limit[ur[x],x-> anew,Direction->"FromBelow"]-(anew-1) >10^-10,Print["error in Condition:", Limit[ur[x],x-> anew,Direction->"FromBelow"]-(anew-1)]];
] 


(* ::Subsection::Initialization:: *)
(*(*Iteration loop*)*)


(* ::Subsubsection::Initialization:: *)
(*(*Initialize*)*)


(* ::Text:: *)
(*Set the equilibrium solution for the hydrogel that is initially in equilibrium with \[Mu]f0*)


(* ::Input::Initialization:: *)
resaeq = FindRoot[(aeq^2-1)/aeq^3==-1/\[CapitalOmega] (1/aeq^3+Log[1-1/aeq^3]-1/(\[Alpha] aeq^3)+\[Chi]/aeq^6)+\[Mu]f0,{aeq,1.1}];
a0 = aeq/.resaeq;
\[Phi]feq = 1-1/a0^3;


(* ::Text:: *)
(*Check if you need to continue from a previous dataset*)


(* ::Input::Initialization:: *)
If[NEW=="Y",  EmptyArray[],ImporterenOudeData[]]
If[NEW == "Y", ind = 1, ind = lengthN+1]


(* ::Subsubsection::Initialization:: *)
(*(*Loop*)*)


(* ::Text:: *)
(*Perform the computations from the flow chart in the correct order. The flow chart is provided in the graduation thesis*)


(* ::Input::Initialization:: *)
Do[{startloop = AbsoluteTime[],
ParametersLoop[],
Initialize\[Phi]fur[],
If[Mod[ind,check]==0,Print["ur(a0)= ",ur[a0]," , a0-1 = ", a0-1.]],
GenerateStressPressure[],
GenerateChemicalPotential[],
SolveSystem[],
Update\[Phi]f[],
SaveData[],
ind=ind+1;
If[Mod[ind,check]==0,{Print["Loop ", ind-1," is finished. It took ",AbsoluteTime[]-startloop , " seconds."]}]
},titer]


(* ::Subsection::Initialization:: *)
(*(*Exporteren*)*)


(* ::Text:: *)
(*Uncomment if you want to export the plot*)


(* ::Input::Initialization:: *)
(*Export["data\[Mu]-excludingdt.mat",Transpose[data\[Mu]][[2]]];
Export["dt-before.mat",Transpose[data\[Mu]][[1]]];

Export["\[Phi]fupdate-excludingdt.mat",Transpose[data\[Phi]farray][[2]]];
Export["dt-update.mat",Transpose[data\[Phi]farray][[1]]];

Export["\[Sigma]rleft-excludingdt.mat",Transpose[data\[Sigma]r][[2]]];
Export["\[Sigma]\[Theta]left-excludingdt.mat",Transpose[data\[Sigma]\[Theta]][[2]]];
Export["\[CapitalPi]left-excludingdt.mat",Transpose[data\[CapitalPi]][[2]]];

Export["ur-excludingdt.mat",Transpose[dataur][[2]]];
Export["a0.mat",dataa0];*)
