(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*Parameters*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
titer =35000;
check =100;


(* ::Input::Initialization:: *)
Ns=10;
dt=1.*10^-8;


(* ::Input::Initialization:: *)
\[Alpha]=250.;
\[CapitalOmega]=1.28*10^-4;
\[Chi]=0.4;
\[Gamma]=1.5;


(* ::Input::Initialization:: *)
RHstart = 0.5;
RHend = 0.9; 


(* ::Input::Initialization:: *)
\[Mu]fs = 1/\[CapitalOmega] Log[RHend];
\[Mu]fs0=\[Mu]fs;
\[Mu]f0 = 1/\[CapitalOmega] Log[RHstart];


(* ::Input::Initialization:: *)
\[Mu]sh = 5*10^-4;
Kl=1;
Kt=1.0;


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*Functions*)*)*)*)*)*)*)


(* ::Subsubsection:: *)
(*Chemische potential and Kt transition*)


(* ::Text:: *)
(*Set varying RH by 1% *)


(* ::Input::Initialization:: *)
steps = 1;  (*either in 1 step or in (0.9-0.5)/100 = 40 steps*)
RHsteps = Table[RHstart + (RHend-RHstart)i/steps,{i,1,steps}] (*height of each step*)


(* ::Text:: *)
(*Function for RH(t) and \[Mu](t)*)


(* ::Input::Initialization:: *)
RH[t_]:=Sum[RHsteps[[j]]*(HeavisideTheta[t-(j-1)/steps-10^-10]-HeavisideTheta[t-(j)/steps-10^-10]),{j,1,steps}] (*RH is discrete so we defined it as a function of repeated step functions. small shift since RH is discontinuous at the exact transition point*)
\[Mu]f[t_]:=1/\[CapitalOmega] Log[RH[t]]; (*conversion from RH to Subscript[\[Mu], f]*)


(* ::Text:: *)
(*Initialize the parameters and varying Kt*)


(* ::Input::Initialization:: *)
ParametersLoop[]:=Module[{rgrid},
If[ind==1,{at=a0,\[Beta]t=\[Beta]0,dadtguess=1000.},{at=an,\[Beta]t=\[Beta]n}]; (*set parameters at time t*)

(*Uncomment line below if one wants to implement manual transition*)
(*If[at>1.03*a0,Kt=0.5];*)


\[Mu]fs = \[Mu]f[ind/titer]; (* set the environment *)

\[CapitalDelta]r = at/Ns;
rgrid = Table[Subscript[r, i]=(i-1/2)*\[CapitalDelta]r,{i,1/2,Ns+1/2,1/2}];
]


(* ::Subsubsection::Initialization:: *)
(*(*(*(*(*(*(*Functions*)*)*)*)*)*)*)


(* ::Text:: *)
(*Standard functions*)


(* ::Input::Initialization:: *)
ur[r_,\[Beta]_]:= r-(r^2-2Integrate[\[Phi]f1[\[Rho]]\[Rho],{\[Rho],0,r}])^(1/2);
\[Lambda]r[r_,\[Beta]_]:=(1-ur[x,b]/x)/(1-\[Phi]f1[x])/.{x-> r,b-> \[Beta]};
\[Lambda]\[Theta][r_,\[Beta]_]:=1/(1-ur[x,b]/x)/.{x-> r,b-> \[Beta]};
\[Lambda]z[r_,\[Beta]_]:=\[Beta];
J[r_,\[Beta]_]:=b/(1-\[Phi]f1[x])/.{x-> r,b-> \[Beta]};
\[CapitalPi][r_,\[Beta]_]:=-1/\[CapitalOmega] (1/J[x,b]+Log[1-1/J[x,b]]-1/(\[Alpha] J[x,b])+\[Chi]/J[x,b]^2)/.{x-> r,b-> \[Beta]};
k[\[Phi]_,\[Beta]_]:=(1-1/b+\[Phi]/b)/(1/b-\[Phi]/b)^\[Gamma]/.{b-> \[Beta]};


(* ::Text:: *)
(*Anisotropic effective stress functions*)


(* ::Input::Initialization:: *)
\[Sigma]rp[r_,\[Beta]_]:=(\[Lambda]r[x,b](3Kt(\[Lambda]r[x,b]-1)+2(2\[Lambda]r[x,b]-\[Lambda]\[Theta][x,b]-\[Lambda]z[x,b])\[Mu]sh))/(3 J[x,b])/.{x-> r,b-> \[Beta]};

\[Sigma]\[Theta]p[r_,\[Beta]_]:= (\[Lambda]\[Theta][x,b](3Kt(\[Lambda]\[Theta][x,b]-1)-2(\[Lambda]r[x,b]-2\[Lambda]\[Theta][x,b]+\[Lambda]z[x,b])\[Mu]sh))/(3J[x,b])/.{x-> r,b-> \[Beta]};

\[Sigma]zp[r_,\[Beta]_]:=(\[Lambda]z[x,b](3Kl(\[Lambda]z[x,b]-1)-2(\[Lambda]r[x,b]+\[Lambda]\[Theta][x,b]-2\[Lambda]z[x,b])\[Mu]sh))/(3J[x,b])/.{x-> r,b-> \[Beta]};


(* ::Subsubsection::Initialization:: *)
(*(*(*(*(*(*(*Modules*)*)*)*)*)*)*)


(* ::Text:: *)
(*Data storage*)


(* ::Input::Initialization:: *)
EmptyArray[]:=Module[{},
dataur={};
data\[Sigma]r={};
data\[Sigma]\[Theta]={};
data\[Sigma]z={};
data\[CapitalPi]={};
datap={};
dataJ={};
data\[Mu]={};
data\[Phi]real={};
data\[Phi]1={};
data\[Mu]env = {};
dataa={{0,a0}};
data\[Beta]={{0,\[Beta]0}};
data\[Beta]arr={}; 
dataVfluid ={{0,\[Pi] a0^2 \[Beta]0 \[Phi]f0}};
dataVtot={{0,\[Pi] a0^2 \[Beta]0}};
dataVpol={{0,\[Pi] a0^2 \[Beta]0 (1-\[Phi]f0)}};
]


(* ::Input::Initialization:: *)
SaveData[]:=Module[{},
AppendTo[dataur,{(ind-1)*dt,listur}];
AppendTo[data\[Sigma]r,{(ind-1)*dt,list\[Sigma]r}];
AppendTo[data\[Sigma]\[Theta],{(ind-1)*dt,list\[Sigma]\[Theta]}];
AppendTo[data\[Sigma]z,{(ind-1)*dt,list\[Sigma]z}];
AppendTo[data\[CapitalPi],{(ind-1)*dt,list\[CapitalPi]}];
AppendTo[dataJ,{(ind-1)*dt,listJ}];
AppendTo[datap,{(ind-1)*dt,listp}];
AppendTo[data\[Mu],{(ind-1)*dt,list\[Mu]}];

AppendTo[data\[Mu]env,{ind*dt,\[Mu]fs}];

AppendTo[data\[Phi]real,{ind*dt,\[Phi]frealarr}];
AppendTo[data\[Phi]1,{ind*dt,\[Phi]f1arr}];
AppendTo[dataa,{ind*dt,an}];
AppendTo[data\[Beta],{ind*dt,\[Beta]n}];
AppendTo[data\[Beta]arr,{ind*dt,\[Beta]newarr}];
AppendTo[dataVfluid,{ind*dt,Vfluid}];
AppendTo[dataVtot,{ind*dt,Vtot}];
AppendTo[dataVpol,{ind*dt,Vpol}];
]


(* ::Text:: *)
(*Initialization of \[Phi]f(t) and calculation of Subscript[u, r](r)*)


(* ::Input::Initialization:: *)
Initialize\[Phi]ur[]:=Module[{},
If[ind==1,{\[Phi]f1[x_]:=(1-\[Phi]f10)HeavisideTheta[x-at]+\[Phi]f10}, {\[Phi]f1[x_]:=(Sum[\[Phi]f1arr[[i,2]]*(HeavisideTheta[x-(i-1)*\[CapitalDelta]r]-HeavisideTheta[x-i*\[CapitalDelta]r]),{i,1,Ns}]+HeavisideTheta[x-at])}];
listur=Table[{(n/Ns)at,ur[(n/Ns)*at,\[Beta]t]},{n,1,Ns}];
PrependTo[listur,{0,0}];
AppendTo[listur,{at,at-1}];
]


(* ::Text:: *)
(*Calculation of the effective stresses and osmotic pressure. Note that it is done manually because it decrease the computation time per loop*)


(* ::Input::Initialization:: *)
GenerateStressOsmotic[]:=Module[{rgrid},
list\[Lambda]r=Table[{Subscript[r, i],\[Lambda]r[Subscript[r, i],\[Beta]t]},{i,1,Ns}];
list\[Lambda]\[Theta]=Table[{Subscript[r, i],\[Lambda]\[Theta][Subscript[r, i],\[Beta]t]},{i,1,Ns}];
list\[Lambda]z=Table[{Subscript[r, i],\[Lambda]z[Subscript[r, i],\[Beta]t]},{i,1,Ns}];
listJ=Table[{Subscript[r, i],J[Subscript[r, i],\[Beta]t]},{i,1,Ns}];

list\[Sigma]r=Table[{Subscript[r, i],(list\[Lambda]r[[i,2]]*(3Kt(list\[Lambda]r[[i,2]]-1)+2\[Mu]sh(2list\[Lambda]r[[i,2]]-list\[Lambda]\[Theta][[i,2]]-list\[Lambda]z[[i,2]])))/(3 listJ[[i,2]])},{i,1,Ns}];
list\[Sigma]\[Theta]=Table[{Subscript[r, i], (list\[Lambda]\[Theta][[i,2]](3Kt(list\[Lambda]\[Theta][[i,2]]-1)-2(list\[Lambda]r[[i,2]]-2list\[Lambda]\[Theta][[i,2]]+list\[Lambda]z[[i,2]])\[Mu]sh))/(3listJ[[i,2]])},{i,1,Ns}];
list\[Sigma]z=Table[{Subscript[r, i],(list\[Lambda]z[[i,2]](3Kl(list\[Lambda]z[[i,2]]-1)-2(list\[Lambda]r[[i,2]]+list\[Lambda]\[Theta][[i,2]]-2list\[Lambda]z[[i,2]])\[Mu]sh))/(3listJ[[i,2]])},{i,1,Ns}];
list\[CapitalPi]=Table[{Subscript[r, i],-1/\[CapitalOmega] (1/listJ[[i,2]]+Log[1-1/listJ[[i,2]]]-1/(\[Alpha] listJ[[i,2]])+\[Chi]/listJ[[i,2]]^2)},{i,1,Ns}];
]


(* ::Text:: *)
(*Calculation of the pore pressure, chemical potential and Subscript[\[Beta], new] via the methods outlined in the report*)


(* ::Input::Initialization:: *)
GenerateChemicalPotential[]:=Module[{dpdr,pedge,pcenter,\[Mu]center},
dpdr=Table[{i*\[CapitalDelta]r,(list\[Sigma]r[[i+1,2]]-list\[Sigma]r[[i,2]])/\[CapitalDelta]r+(1/2 (list\[Sigma]r[[i,2]]+list\[Sigma]r[[i+1,2]])-1/2 (list\[Sigma]\[Theta][[i,2]]+list\[Sigma]\[Theta][[i+1,2]]))/Subscript[r, i+1/2]},{i,1,Ns-1}];
pedge = list\[Sigma]r[[-1]]; (*\[Sigma]r = 0 at edge. edge approx equal to Subscript[r, Ns] (center of last segment). \[Sigma]r=0 \[Rule] p=\[Sigma]r'*)
pcenter = pedge[[2]] - Sum[dpdr[[i,2]]*\[CapitalDelta]r,{i,1,Ns-1}]; (*use Taylor approx p(r+\[CapitalDelta]r) = p(r)+dpdr|r \[CapitalDelta]r repeated*)
pr = Table[{(j+1/2)*\[CapitalDelta]r,pcenter + Sum[dpdr[[i,2]]*\[CapitalDelta]r,{i,1,j}]},{j,1,Ns-1}]; (*j is the index of the element. found p(Subscript[r, 1]) as pcenter. so we only need to calculate elements 2 until Ns*)
PrependTo[pr,{1/2*\[CapitalDelta]r,pcenter}];
listp = pr;

\[Beta]newarr = Table[{(i-1/2)*\[CapitalDelta]r,\[Beta]new/.FindRoot[\[Sigma]zp[Subscript[r, i],\[Beta]new]==pr[[i,2]],{\[Beta]new,\[Beta]t},WorkingPrecision->75,MaxIterations->1000,AccuracyGoal->5,PrecisionGoal->6]},{i,1,Ns}];(*calculation of \[Beta](r)*)

\[Beta]new = Mean[\[Beta]newarr][[2]]; (*average*)

\[Beta]new = 1/2 (\[Beta]new+\[Beta]t); (*suppress fluctuations*)

\[Mu]center = pcenter-list\[CapitalPi][[1,2]]; (*condition \[Mu]=p-\[CapitalPi]*)
d\[Mu]dr=Table[{i*\[CapitalDelta]r,dpdr[[i,2]]-(list\[CapitalPi][[i+1,2]]-list\[CapitalPi][[i,2]])/\[CapitalDelta]r},{i,1,Ns-1}]; (*mechanical equilibrium equation*)
\[Mu]=Table[{(j+1/2)*\[CapitalDelta]r,\[Mu]center + Sum[d\[Mu]dr[[i,2]]*\[CapitalDelta]r,{i,1,j}]},{j,1,Ns-1}];(*use Taylor approx repeated*)
PrependTo[\[Mu],{1/2 \[CapitalDelta]r,\[Mu]center}];
AppendTo[\[Mu],{(Ns+1/2)*\[CapitalDelta]r,\[Mu]fs}];(*\[Mu] in environment*)
AppendTo[d\[Mu]dr,{\[CapitalDelta]r*Ns,(\[Mu][[-1,2]]-\[Mu][[-2,2]])/\[CapitalDelta]r}]; (*d\[Mu]dr over the outerboundary *)
list\[Mu]=\[Mu];
]


(* ::Text:: *)
(*Generation of the solution to the numerical equations*)


(* ::Input::Initialization:: *)
SolveSystem[]:=Module[{},
Table[Subscript[\[Phi]f1old, i]=\[Phi]f1[Subscript[r, i]],{i,1,Ns}]; Subscript[\[Phi]f1old, Ns+1]=1; (*environment is 1*)

(*here, the equations are definied*)
EQS = {\[Pi](Subscript[r, 3/2]^2-Subscript[r, 1/2]^2)*((Subscript[\[Phi]f1new, 1]-Subscript[\[Phi]f1old, 1])+(1-Subscript[\[Phi]f1old, 1])/\[Beta]t*(\[Beta]new-\[Beta]t)+(2Subscript[\[Phi]f1old, 1])/at*(anew-at))-2\[Pi]((Subscript[r, 3/2]^2*1/2 (Subscript[\[Phi]f1old, 1]+Subscript[\[Phi]f1old, 2]))/at*(anew-at)+Subscript[r, 3/2]*(1-1/2 (Subscript[\[Phi]f1old, 2]+Subscript[\[Phi]f1old, 1]))k[1/2 (Subscript[\[Phi]f1old, 2]+Subscript[\[Phi]f1old, 1]),\[Beta]t]*d\[Mu]dr[[1,2]]*dt)==0};
For[i=2,i<= Ns,i++,
AppendTo[EQS, 
\[Pi](Subscript[r, i+1/2]^2-Subscript[r, i-1/2]^2)*((Subscript[\[Phi]f1new, i]-Subscript[\[Phi]f1old, i])+(1-Subscript[\[Phi]f1old, i])/\[Beta]t*(\[Beta]new-\[Beta]t)+(2Subscript[\[Phi]f1old, i])/at*(anew-at))-2\[Pi]((Subscript[r, i+1/2]^2*1/2 (Subscript[\[Phi]f1old, i]+Subscript[\[Phi]f1old, i+1]))/at*(anew-at)+Subscript[r, i+1/2]*(1-1/2 (Subscript[\[Phi]f1old, i]+Subscript[\[Phi]f1old, i+1]))k[1/2 (Subscript[\[Phi]f1old, i+1]+Subscript[\[Phi]f1old, i]),\[Beta]t]*d\[Mu]dr[[i,2]]*dt)+2\[Pi]((Subscript[r, i-1/2]^2*1/2 (Subscript[\[Phi]f1old, i]+Subscript[\[Phi]f1old, i-1]))/at*(anew-at)+Subscript[r, i-1/2]*(1-1/2 (Subscript[\[Phi]f1old, i]+Subscript[\[Phi]f1old, i-1]))k[1/2 (Subscript[\[Phi]f1old, i-1]+Subscript[\[Phi]f1old, i]),\[Beta]t]*d\[Mu]dr[[i-1,2]]*dt)==0]];
COND = {(anew)^2 (1-Sum[Subscript[\[Phi]f1new, i]((i/Ns)^2-((i-1)/Ns)^2),{i,1,Ns}])-1==0};
set = Join[EQS,COND];

(*here, the unknowns and initial guesses (previous results) are defined and the equations are solved*)
unknown\[Phi]= Table[{Subscript[\[Phi]f1new, i],\[Phi]f10q },{i,1,Ns}];
unknowna = {{anew,aq}};
listunknown = Join[unknowna,unknown\[Phi]];

sol = FindRoot[set,listunknown,MaxIterations->2500,WorkingPrecision->100];
solution=sol;
]


(* ::Text:: *)
(*Update all quantities at the new moment in time*)


(* ::Input::Initialization:: *)
UpdateAll[]:=Module[{},
an = anew/.solution; 
\[Beta]n=\[Beta]new;
\[CapitalDelta]rn=an/Ns;

\[Phi]f1arr = Table[{(i-1/2)\[CapitalDelta]rn,Subscript[\[Phi]f1new, i]},{i,1,Ns}]/.solution;
\[Phi]frealarr = Table[{(i-1/2)*\[CapitalDelta]rn,1-1/\[Beta]n+1/\[Beta]n Subscript[\[Phi]f1new, i]},{i,1,Ns}]/.solution;
\[Phi]freal[x_]:=Sum[\[Phi]frealarr[[i,2]]*(HeavisideTheta[x-(i-1)*\[CapitalDelta]rn]-HeavisideTheta[x-i*\[CapitalDelta]rn]),{i,1,Ns}]+HeavisideTheta[x-an];

(*calculate volumes at time t+dt*)
Vfluid = Integrate[\[Phi]freal[x]*2\[Pi] x,{x,0,an},{z,0,\[Beta]n}];
Vtot = Integrate[2\[Pi] x,{x,0,an},{z,0,\[Beta]n}];
Vpol = Vtot-Vfluid;

]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*(*(*Solving*)*)*)*)*)*)*)


(* ::Input::Initialization:: *)
\[Beta]d=1;ad=1;ind=1;EmptyArray[]


(* ::Subsubsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*Initialize*)*)*)*)*)*)*)


(* ::Text:: *)
(*Solution for energy minimum for \[Mu]0*)


(* ::Input::Initialization:: *)
reseq0 =FindRoot[{(3Kt(aeq-1)+2(aeq-\[Beta]eq)\[Mu]sh)/(3aeq \[Beta]eq)==-1/\[CapitalOmega] (1/(aeq^2 \[Beta]eq)+Log[1-1/(aeq^2 \[Beta]eq)]-1/(\[Alpha] aeq^2 \[Beta]eq)+\[Chi]/(aeq^2 \[Beta]eq)^2)+\[Mu]f0,(3Kl(\[Beta]eq-1)-4(aeq-\[Beta]eq)\[Mu]sh)/(3aeq^2)==-1/\[CapitalOmega] (1/(aeq^2 \[Beta]eq)+Log[1-1/(aeq^2 \[Beta]eq)]-1/(\[Alpha] aeq^2 \[Beta]eq)+\[Chi]/(aeq^2 \[Beta]eq)^2)+\[Mu]f0},{{aeq,1.1},{\[Beta]eq,1.1}},MaxIterations->1000]

a0=aeq/.reseq0;\[Beta]0 =\[Beta]eq/.reseq0; 
\[Phi]f0 = 1-(ad^2 \[Beta]d)/(a0^2 \[Beta]0); \[Phi]f10 = 1-\[Beta]0+\[Beta]0 \[Phi]f0;
Print["a0: ", a0,". \[Beta]0: ",\[Beta]0,". \[Phi]f0: ",\[Phi]f0,"."]


(* ::Subsubsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*Equilibrium*)*)*)*)*)*)*)


(* ::Text:: *)
(*Solution for energy minimum for \[Mu]s*)


(* ::Input::Initialization:: *)
equilibrium= FindRoot[{((3Kt(aeq-1)+2(aeq-\[Beta]eq)\[Mu]sh)/(3aeq \[Beta]eq))==-1/\[CapitalOmega] (1/(aeq^2 \[Beta]eq)+Log[1-1/(aeq^2 \[Beta]eq)]-1/(\[Alpha] aeq^2 \[Beta]eq)+\[Chi]/(aeq^2 \[Beta]eq)^2)+\[Mu]fs,(3Kl(\[Beta]eq-1)-4(aeq-\[Beta]eq)\[Mu]sh)/(3aeq^2)==-1/\[CapitalOmega] (1/(aeq^2 \[Beta]eq)+Log[1-1/(aeq^2 \[Beta]eq)]-1/(\[Alpha] aeq^2 \[Beta]eq)+\[Chi]/(aeq^2 \[Beta]eq)^2)+\[Mu]fs},{{aeq,1.1},{\[Beta]eq,1.1}},MaxIterations->1000]

aq = aeq/.equilibrium; \[Beta]q=\[Beta]eq/.equilibrium;
\[Phi]f0q = 1-(ad^2 \[Beta]d)/(aq^2 \[Beta]q); \[Phi]f10q = 1-\[Beta]q+\[Beta]q \[Phi]f0q;


(* ::Subsubsection::Initialization:: *)
(*(*(*(*(*(*(*Iterations*)*)*)*)*)*)*)


(* ::Text:: *)
(*Iterative process - see flowchart*)


(* ::Input::Initialization:: *)
Do[{startloop=AbsoluteTime[],
ParametersLoop[],
Initialize\[Phi]ur[],
If[Mod[ind,check]==0,Print["ur(a0)= ",DecimalForm[ur[at,\[Beta]t],6],", a-1 = ", at-1.,". Absolute difference: ",ur[at,\[Beta]t]-(at-1.)]],
GenerateStressOsmotic[],
GenerateChemicalPotential[],
SolveSystem[],
UpdateAll[],
SaveData[],
ind=ind+1,
If[Mod[ind,check]==0,{Print["Loop ", ind-1," is finished. It took ",AbsoluteTime[]-startloop , " seconds."]}]
},titer]


(* ::Subsection::Initialization::Closed:: *)
(*(*(*(*(*(*(*Exporteren Data*)*)*)*)*)*)*)


(* ::Text:: *)
(*Exportation of the data. Note that all elements that are combinations of {timestep,array} are split, and must be recombined in the reading-file*)


(* ::Input::Initialization:: *)
Export["data\[Mu]-excludingdt.mat",Transpose[data\[Mu]][[2]]];
Export["dt-before.mat",Transpose[data\[Mu]][[1]]];

Export["\[Phi]freal-excludingdt.mat",Transpose[data\[Phi]real][[2]]];
Export["\[Phi]f1-excludingdt.mat",Transpose[data\[Phi]1][[2]]];
Export["dt-update.mat",Transpose[data\[Phi]real][[1]]];

Export["\[Sigma]rleft-excludingdt.mat",Transpose[data\[Sigma]r][[2]]];
Export["\[Sigma]\[Theta]left-excludingdt.mat",Transpose[data\[Sigma]\[Theta]][[2]]];
Export["\[Sigma]zleft-excludingdt.mat",Transpose[data\[Sigma]z][[2]]];
Export["\[CapitalPi]left-excludingdt.mat",Transpose[data\[CapitalPi]][[2]]];

Export["ur-excludingdt.mat",Transpose[dataur][[2]]];
Export["a0.mat",dataa];
Export["bz.mat", data\[Beta]];
Export["\[Mu]env.mat",data\[Mu]env];
Export["bzarr.mat",Transpose[data\[Beta]arr][[2]]];

Export["Vtot.mat",Transpose[dataVtot][[2]]];
Export["Vfluid.mat",Transpose[dataVfluid][[2]]];
Export["Vpol.mat",Transpose[dataVpol][[2]]];


(* ::Input::Initialization:: *)
(*Quit[]*)
