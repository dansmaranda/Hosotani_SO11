(* Use MaTeX to generate TeX axis for the plots *)
(* << MaTeX` *)
(* Set up the machine 0*)
machineZero = 10^(-$MachinePrecision);
(*Number of Families, and bound where we consider non-perturbative behaviour of the couplings*)
nbOfFamilies = 1;
(* nonPertBound = 1; *)
timeOutSolve = 30;

(*Initial beta coefficients for the QED & the weinberg Angle RGE running*)
initBetaEM = 22;
initBetaSin2ThW = -1.0864;

(*plotHandle*)
makePlots = False;

SetSystemOptions[
  "MachineRealPrintPrecision" -> Round[$MachinePrecision]];
Clear[zL, k, \[Mu]11, \[Mu]1, \[Mu]11Prime, \[Mu]2Tilde, c0Prime, c0,
  c1, c2, \[Theta]HRule, \[Theta]Hmin, sin2\[Theta]W, \[Xi]Gauge, M,
  mB, \[Alpha]EM, fH];


(* Print["-------------------------"] *)
JsonNb=$ScriptCommandLine[[2]];
(* JsonNb = "ThreadNb-0" *)
jsonName="dataIn"<>JsonNb<>".json";
jsonNameOut="weinbergAngleOut"<>JsonNb<>".json";


(* Model parameters *)
MZ = 91.1876
sin2\[Theta]W = 0.2312;
\[Xi]Gauge = 0;
M = -10^7;
mB = 1.145*10^12;
\[Alpha]EM = 1/127.96;

dataRule = Import[jsonName];
k = "k" /. dataRule;
zL = "zL" /. dataRule;
c0 = "c0" /. dataRule;
c1 = "c1" /. dataRule;
c2 = "c2" /. dataRule;
c0Prime = "c0Prime" /. dataRule;
\[Mu]1 = "Mu1" /. dataRule;
\[Mu]2Tilde = "Mu2Tilde" /. dataRule;
\[Mu]11 = "Mu11" /. dataRule;
\[Mu]11Prime = "Mu11Prime" /. dataRule;
\[Theta]Hmin = "ThetaHiggs"/. dataRule;


(* k = 89130;
zL = 35;
c0 = 0.3325;
c1 = 0.0;
c2 = -0.7;
c0Prime = 0.5224;
\[Theta]Hmin = 0.149481;

\[Mu]1 = 11.18;
\[Mu]2Tilde = 0.7091;
\[Mu]11 = 0.108;
\[Mu]11Prime = 0.108;

mKK5 = \[Pi] k /(zL - 1);
L5 = N[Log[zL]/(k)]; *)

(* Print["Loaded in parameters, starting analysis."] *)
Print["------------ Loaded Params; Starting 4D Analysis ---------------"]

(*------------------- DEFINITIONS --------------------*)
(*Bessel basis functions Subscript[FHat, \
\[Alpha]\[Beta]],Subscript[F, \[Alpha]\[Beta]] *)
Fab[\[Alpha]_, \[Beta]_, u_, v_] :=
  BesselJ[\[Alpha], u] BesselY[\[Beta], v] -
   BesselY[\[Alpha], u] BesselJ[\[Beta], v];
FHat[\[Alpha]_, \[Beta]_, u_, v_] :=
  BesselI[\[Alpha], u] * BesselK[\[Beta], v] -
   Exp[-I (\[Alpha] - \[Beta] ) * Pi] * BesselK[\[Alpha], u] *
    BesselI[\[Beta], v];
FHatcPM[q_, c_, pm1_, pm2_, zL_] :=
  FHat[ c + pm1*1/2, c + pm2*1/2, q*(zL^(-1)), q  ];

(*   Fermions  *)
SL[\[Lambda]_, z_, zL_, c_] := -\[Pi]/
   2 \[Lambda] Sqrt[z zL] Fab[c + 1/2,
    c + 1/2, \[Lambda] z, \[Lambda] zL];
SR[\[Lambda]_, z_, zL_, c_] := +\[Pi]/
   2 \[Lambda] Sqrt[z zL] Fab[c - 1/2,
    c - 1/2, \[Lambda] z, \[Lambda] zL];
CR[\[Lambda]_, z_, zL_, c_] := -\[Pi]/
   2 \[Lambda] Sqrt[z zL] Fab[c - 1/2,
    c + 1/2, \[Lambda] z, \[Lambda] zL];
CL[\[Lambda]_, z_, zL_, c_] := +\[Pi]/
   2 \[Lambda] Sqrt[z zL] Fab[c + 1/2,
    c - 1/2, \[Lambda] z, \[Lambda] zL];
(*   Bosons   *)
Cbasis[\[Lambda]_, z_, zL_] := \[Pi]/
   2 \[Lambda] (z^(3/2)) (zL^(1/2)) Fab[3/2,
    1/2, \[Lambda] z, \[Lambda] zL];
Sbasis[\[Lambda]_, z_, zL_] := -\[Pi]/
   2 \[Lambda] (z^(3/2)) (zL^(1/2)) Fab[3/2,
    3/2, \[Lambda] z , \[Lambda] zL];
Cprime[\[Lambda]_, z_, zL_] := \[Pi]/
   2 (\[Lambda]^2) (z^(3/2)) (zL^(1/2)) Fab[1/2,
    1/2, \[Lambda] z, \[Lambda] zL];
Sprime[\[Lambda]_, z_, zL_] := -\[Pi]/
   2 (\[Lambda]^2) (z^(3/2)) (zL^(1/2)) Fab[1/2,
    3/2, \[Lambda] z, \[Lambda] zL];

(*Mass Approximations*)
mWboson[k_, zL_, \[Theta]H_] :=
  Sqrt[3/2] * k *((zL)^(-3/2)) * Sin[\[Theta]H];
mZboson[mWboson_, sin2\[Theta]W_] := mWboson / Sqrt[1 - sin2\[Theta]W];
muTypeQuark[c0_, k_, zL_, \[Theta]H_ ] := If[c0 < 0.5,
    (k/zL )*Sqrt[1 - 4 (c0)^2] Sin[\[Theta]H/2],
    k *((zL)^-1/2 - c0)*Sqrt[4 (c0)^2 - 1] Sin[\[Theta]H/2]
   ];
mdTypeQuark[c0_, c1_, \[Mu]11_, \[Mu]1_, k_, zL_, \[Theta]H_] :=
  If[c0 < 0.5,
   k * ((zL)^(c0 - c1 - 1)) * \[Mu]11/\[Mu]1*
    Sqrt[(1 - 2 c0) (1 + 2 c1)]*Sin[\[Theta]H/2],
   k * ((zL)^(-c1 - 1/2)) * \[Mu]11/\[Mu]1*
    Sqrt[(2 c0 - 1) (2 c1 + 1)]*Sin[\[Theta]H/2]
   ];
meTypeLepton[c2_, c0_, \[Mu]11Prime_, \[Mu]2Tilde_, k_,
   zL_, \[Theta]H_] := If[ c2 < 0.5,
   k * ((zL)^(-1 + c2 - c0)) * \[Mu]11Prime/\[Mu]2Tilde Sqrt[(2 c0 +
        1) (1 - 2 c2)] Sin[\[Theta]H/2],
   k *  ((zL)^(-1/2 - c0)) * \[Mu]11Prime/\[Mu]2Tilde Sqrt[(2 c0 +
        1) (2 c2 - 1)] Sin[\[Theta]H/2]
   ];
mPsiDarkType[c0Prime_, k_, zL_, \[Theta]H_] := If[c0Prime < 0.5,
    k *(zL^(-1)) Sqrt[1 - 4 (c0Prime)^2] Cos[\[Theta]H/2],
    k *(zL^(-1/2 - c0Prime)) Sqrt[4 (c0Prime)^2 - 1] Cos[\[Theta]H/2]
   ];
m\[Nu]Type [mu_, M_, mB_, c0_, zL_] := -((
   2 (mu)^2 *M* ((zL)^(2 c0 + 1)))/((2 c0 + 1) *(mB^2)));




   replacementRules = {z -> 1, c0var -> c0,
      zLvar -> zL, \[Mu]1var -> \[Mu]1, \[Mu]11var -> \[Mu]11,
      c1var ->
       c1, \[Mu]2TildeVar -> \[Mu]2Tilde, \[Mu]11PrimeVar -> \
   \[Mu]11Prime, c2var -> c2, c0PrimeVar -> c0Prime};

   (* Fermion basis functions for c0, c0Prime, c1, c2 *)

   SL0 = SL[\[Lambda], z, zLvar, c0var];
   SR0 = SR[\[Lambda], z, zLvar, c0var];
   CR0 = CR[\[Lambda], z, zLvar, c0var];
   CL0 =  CL[\[Lambda], z, zLvar, c0var];
   SL1 = SL[\[Lambda], z, zLvar, c1var];
   CR1 = CR[\[Lambda], z, zLvar, c1var];
   SR2 = SR[\[Lambda], z, zLvar, c2var];
   CL2 = CL[\[Lambda], z, zLvar, c2var];
   SL0prime = SL[\[Lambda], z, zLvar, c0PrimeVar];
   SR0prime = SR[\[Lambda], z, zLvar, c0PrimeVar];


   (* Standard Model Fermion and Boson towers*)

   uQuarkSol = SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 /. replacementRules;
   bQuarkSol =
     SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 + ((\[Mu]1var^2) SR0* CR0 *
          SL1 *CR1)/(\[Mu]11var^2 (CR1)^2 - (SL1)^2) /.
      replacementRules;
   tauLeptonSol =
     SL0 * SR0  + (Sin[\[Theta]Hmin/
           2])^2 + (\[Mu]2TildeVar^2 SL0 CL0 SR2 CL2) / \
   (\[Mu]11PrimeVar^2 (CL2)^2  - (SR2)^2) /. replacementRules;

   WbosonSol = (2* Cprime[\[Lambda], z, zLvar]*
         Sbasis[\[Lambda], z,
          zLvar] + \[Lambda] (Sin[\[Theta]Hmin])^2) /. replacementRules;
   ZbosonSol = (2* Cprime[\[Lambda], z, zLvar]*
         Sbasis[\[Lambda], z, zLvar] +
        1/(cos2\[Theta]W) \[Lambda] (Sin[\[Theta]Hmin])^2) /.
      replacementRules;



   (* Exotic Fields: Fermion and Boson towers*)
   (*The EXACT soln for the \
   \[Gamma]Tower is \[Pi]/(zL-1), don't need to solve it , and Sbasis \
   has the 2nd solution after Subscript[\[Lambda], \[Gamma]]*)

   psiDarkQuarkSol =
     SL0prime * SR0prime + (Cos[\[Theta]Hmin/2])^2 /. replacementRules;

   WRSol = Cbasis[\[Lambda], z, zLvar] /. replacementRules;
   \[Gamma]Tower = Cprime[\[Lambda], z, zLvar] /. replacementRules;
   HiggsSolTower = Sbasis[\[Lambda], z, zLvar] /. replacementRules;

   NeutrinoSol = -(k * \[Lambda] + M)/(mB) (
         SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2) - (mB/(2 k))*SR0 CR0 /.
      replacementRules;


      cutOffReg = \[Pi] /(zL - 1);

      (*Solving the tau tower.*)

      tauGuess =
        meTypeLepton[c2, c0, \[Mu]11Prime, \[Mu]2Tilde, k,
          zL, \[Theta]Hmin] / k;
      \[Lambda]ListTauGuess =
        TimeConstrained[
         NSolve[tauLeptonSol == 0 &&
           tauGuess >= \[Lambda] >= machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];
      \[Lambda]ListTauFullRange  =
        TimeConstrained[
         NSolve[tauLeptonSol == 0 &&
           cutOffReg >= \[Lambda] >= tauGuess, \[Lambda] \[Element] Reals],
         timeOutSolve];
      \[Lambda]ListTau =
        DeleteDuplicates[
         Join[\[Lambda]ListTauGuess, \[Lambda]ListTauFullRange]];

      (*Solving the Bottom tower.*)

      bottomGuess =
        mdTypeQuark[c0, c1, \[Mu]11, \[Mu]1, k, zL, \[Theta]Hmin]/k;
      \[Lambda]ListBottomGuess =
        TimeConstrained[
         NSolve[bQuarkSol == 0 &&
           bottomGuess >= \[Lambda] >= machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];
      \[Lambda]ListBottomFullRange  =
        TimeConstrained[
         NSolve[bQuarkSol == 0 &&
           cutOffReg >= \[Lambda] >= bottomGuess, \[Lambda] \[Element]
           Reals], timeOutSolve];
      \[Lambda]ListBottom =
        DeleteDuplicates[
         Join[\[Lambda]ListBottomGuess, \[Lambda]ListBottomFullRange]];


      (*Solving Top Tower*)
      \[Lambda]ListTop =
        TimeConstrained[
         NSolve[uQuarkSol == 0 &&
           cutOffReg >= \[Lambda] > machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];

      (*Solving Psi Dark Tower*)
      \[Lambda]ListPsiDark =
        TimeConstrained[
         NSolve[psiDarkQuarkSol == 0 &&
           cutOffReg >= \[Lambda] > machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];


      (*Solving W, Z Boson Tower *)
      \[Lambda]ListW =
        TimeConstrained[
         NSolve[WbosonSol == 0 &&
           cutOffReg >= \[Lambda] > machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];
      \[Lambda]ListZ =
        TimeConstrained[
         NSolve[ZbosonSol == 0 &&
           cutOffReg >= \[Lambda] > machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];

      (*Solving WR Boson Tower (Note that it has the same equations as X \
      gluons, and ZR)*)
      \[Lambda]ListWR =
        TimeConstrained[
         NSolve[WRSol == 0 &&
           cutOffReg >= \[Lambda] > machineZero, \[Lambda] \[Element]
           Reals], timeOutSolve];

      (*Solving the Neutrino Tower*)
      \[Lambda]List\[Nu] =
        TimeConstrained[
         NSolve[NeutrinoSol == 0 &&
           2*mKK5/k >= \[Lambda] > machineZero, \[Lambda] \[Element] Reals],
          timeOutSolve];



          makeSUNassocDict[n_] := Module[
             {dRDict, casKey, dynKey, descKey},
             casKey = "QuadCasimir";
             dynKey = "DynkinIdx";
             descKey = "IrrepDesc";



             dRDict = <|
               "1" -> <|casKey -> 0 ,  dynKey -> 0,
                 descKey -> "Scalar Irrep (0000...0000)"|> |>;

             If[n >= 2,
              AppendTo[dRDict,
               ToString[n] -> <|casKey -> (n^2 - 1)/2 n ,   dynKey -> 1/2,
                 descKey -> "Fundamental Irrep (1000...0000)"|>];

              AppendTo[dRDict,
               ToString[n * (n + 1)/2] ->   <|casKey -> (n - 1) (n + 2)/2,
                 dynKey -> (n + 2)/2,
                 descKey -> "Irrep (2000...0000)"|>];

              AppendTo[dRDict,
               ToString[n^2 - 1] ->             <|casKey -> n, dynKey -> n,
                 descKey -> "Adjoint Irrep (1000...0001)"|>];


              ];

             If[n >= 3,
              AppendTo[dRDict,
                ToString[n * (n - 1)/2] ->   <|casKey -> (n + 1) (n - 2)/n ,
                  dynKey -> (n - 2)/2, descKey -> "Irrep (0100...0000)"|>];
              ];
             If[n >= 5,
              AppendTo[dRDict,
               ToString[n * (n - 1) (n - 2)/6] ->   <|
                 casKey -> (n - 2) (n - 3) (n - 4)/12,
                 dynKey -> (n - 3) (n - 2)/4, descKey -> "Irrep (0010...0000)"|>];
              AppendTo[dRDict,
               ToString[n * (n - 1) (n + 1)/3] ->   <|
                 casKey -> 3 (n^2 - 3)/(2 n), dynKey -> (n^2 - 3)/2,
                 descKey -> "Irrep (1100...0000)"|>];
              AppendTo[dRDict,
               ToString[n * (n + 1) (n + 2)/6] ->   <|
                 casKey -> (n + 2) (n + 3) (n + 4)/12,
                 dynKey -> (n + 2) (n + 3)/4, descKey -> "Irrep (3000...0000)"|>];
              AppendTo[dRDict,
               ToString[n^2 * (n + 1) (n - 1)/12] ->   <|
                 casKey -> n (n^2 - 16)/3, dynKey -> n (n - 2) (n + 2)/6,
                 descKey -> "Irrep (0200...0000)"|>];
              ];


             Return[dRDict]
             ];


             (* Makes a dictionary with all the names of the towers and their \
             masses*)
             makeDictTowers[\[Lambda]List_, \[Lambda]Name_] := Module[
                {i, unpackedDitcList},

                unpackedDitcList = <||>;
                For[i = 1, i <= Length[\[Lambda]List], i++,
                 unpackedDitcList[\[Lambda]Name <> ToString[i]] = \[Lambda]List[[
                    i]];];

                Return[unpackedDitcList];];


             (* Adds a contribution to the beta coefficient depending on weather \
             we're adding a fermion or a boson*)

             addToBetaCoeff[currBetaCoef_,  dofToAdd_, SUN_, SUNcharge_] := Block[
                {DynkinIndex , toAddBetaCoeff, QuadraticCasimir, irrepDictSUN,
                 Ta\[Psi], strippedSUN, C2SUN},
                 (* Printing here. *)
                 (* Print["SU3C, adding ", dofToAdd]; *)

                If [SUN != "U1EM",
                 strippedSUN =
                  ToExpression[
                   StringCases[SUN, RegularExpression["\\d+"]][[1]] ];
                 irrepDictSUN = makeSUNassocDict[strippedSUN];

                 DynkinIndex = irrepDictSUN[[ ToString[SUNcharge] ]][["DynkinIdx"]];
                 QuadraticCasimir =
                  irrepDictSUN[[ ToString[SUNcharge] ]][["QuadCasimir"]];
                 (*Print[dofToAdd, " has a Dynkin index of ", Ta\[Psi],
                 " and a Quadratic Casimir ",C2SUN , " under ", SUN];*)


                 (*CASIMIR HERE !!!*)
                 (*DynkinIndex=1/2;
                 QuadraticCasimir= 1;*)


                 If[dofToAdd == "Fermion",
                  toAddBetaCoeff = 4/3 * 1/2*DynkinIndex*nbOfFamilies;
                  (* Printing here. *)
                  (* Print["SU3C:", currBetaCoef + toAddBetaCoeff]; *)
                  Return[currBetaCoef + toAddBetaCoeff],


                  If[ dofToAdd == "Boson",
                    toAddBetaCoeff = -(11/3) * QuadraticCasimir ;
                    (* Printing here. *)
                    (* Print["SU3C:", currBetaCoef + toAddBetaCoeff]; *)
                    Return[currBetaCoef + toAddBetaCoeff]];]


                 ,
                 Print["Need to run separate U1EM Routine!"];

                 ]];


             addToBetaCoeffU1EM[currBetaCoef_,  dofToAdd_, U1EMDict_, dofName_,
                colorFactor_] := Block[
                {toAddBetaCoeff,  gammaFactor, chargeFactor},

                If[dofToAdd == "Fermion",
                  (* Print["Adding Fermion"]; *)
                  (*All the fremions are chiral*)

                  gammaFactor = 8;
                  chargeFactor = U1EMDict[["QEM"]]^2;

                  (* Printing here. *)
                  (* Print[chargeFactor, "U1EM: ", dofName, "  ",
                   currBetaCoef +
                    colorFactor * gammaFactor * chargeFactor * nbOfFamilies]; *)
                  Return[
                   currBetaCoef +
                    colorFactor * gammaFactor * chargeFactor * nbOfFamilies];

                  ,


                  If[ dofToAdd == "Boson",
                    (* Print["Adding Boson"]; *)

                    gammaFactor = -22;
                    chargeFactor = U1EMDict[[ "QEM"]]^2;
                    (* Printing here. *)
                    (* Print["U1EM:", dofName, "  ",currBetaCoef +
                      colorFactor * gammaFactor * chargeFactor]; *)
                    Return[
                     currBetaCoef + colorFactor * gammaFactor * chargeFactor]
                    ];
                  ];

                ];

             makeOrderedTowerDict[classDict_] := Block[
               {dictMerged, allPartNames, partName, i, newDict},
               dictMerged = <||>;
               allPartNames = Keys[classDict];

               For [i = 1, i <= Length[allPartNames], i++,
                partName = allPartNames[[i]];
                newDict =
                 makeDictTowers[classDict[[partName]][["Sols"]], partName ];
                AppendTo[dictMerged, newDict];
                ];
               Return[ SortBy[dictMerged, # &]];
               ]



               makePiecewiseRGE[classDict_, gFct_, startBetaCoeff_,
                  SUN_,  \[Mu]Renorm_, k_, startReg_: 80] := Block[
                  { \[Beta]FctNoCoeff, startBetaCoeffTemp, newBetaCoeff, regList,
                   dictKeys, \[Beta]CoeffList, pieceWiseList, matterType, partIndex,
                   dictMerged, allPartNames, i, partName, newDict, orderedDict,
                   SUNcharge, strippedSUN},

                  (* Convert the full dictionary into a ordered dictionary,
                  sorted by the mass value of each state *)
                  dictMerged = <||>;
                  allPartNames = Keys[classDict];

                  For [i = 1, i <= Length[allPartNames], i++,
                   partName = allPartNames[[i]];
                   newDict =
                    makeDictTowers[classDict[[partName]][["Sols"]], partName ];
                   AppendTo[dictMerged, newDict];
                   ];
                  orderedDict = SortBy[dictMerged, # &];

                  (* Building the g's*)
                  If[ SUN == "U1EM",
                   \[Beta]FctNoCoeff = 1/(4 \[Pi])^2 (gFct[\[Mu]Renorm]^3) /6
                   ,
                   \[Beta]FctNoCoeff = gFct[\[Mu]Renorm]^3/(4 \[Pi])^2;
                   ];
                  startBetaCoeffTemp = startBetaCoeff;


                  dictKeys = Keys[orderedDict];
                  regList = {startReg};
                  \[Beta]CoeffList = {startBetaCoeff};



                  For[partIndex = 1, partIndex <= Length[ dictKeys], partIndex++,

                   partName =
                    StringTrim[dictKeys[[partIndex]], RegularExpression["\\d+$"]];
                   SUNcharge = classDict[[ partName]][["GSMCharges"]][[SUN]]  ;

                   (*Making the \[Beta] fct coefficients for the various regions*)

                     matterType = classDict[[ partName ]][["Type"]];
                   (*newBetaCoeff = addToBetaCoeff[startBetaCoeffTemp, matterType ,
                   SUN, SUNcharge];*)

                   If[
                     SUN == "U1EM",
                     U1EMDict = classDict[[partName]][["GSMCharges"]][["U1EM"]] ;
                     newBetaCoeff =
                      addToBetaCoeffU1EM[startBetaCoeffTemp, matterType, U1EMDict,
                       partName , classDict[[partName]][["GSMCharges"]][["SU3C"]] ];
                     ,
                     (*Printing here*)
                     (* Print[partName, "  ", SUNcharge, "  ",SUN]; *)
                     newBetaCoeff =
                       addToBetaCoeff[startBetaCoeffTemp, matterType , SUN,
                        SUNcharge];
                     ]


                    (*Condition is here to assure we're not introducing \
               noncontributing factors which mess up the regions*)

                    If [newBetaCoeff != startBetaCoeffTemp,
                     AppendTo[\[Beta]CoeffList, newBetaCoeff];
                     startBetaCoeffTemp = newBetaCoeff;

                     (*Defining the \[Mu] region thresholds*)

                     AppendTo[regList,
                      k*  orderedDict[[  dictKeys[[partIndex]] ]][[1, 2]]   ];
                     ]
                   ];
                  (*Adding the last mKK5 threshold*)
                  AppendTo[regList, N[ mKK5]];

                  pieceWiseList = {};
                  For[regIdx = 1, regIdx <= Length[regList] - 1, regIdx++,
                   AppendTo[pieceWiseList,
                     {\[Beta]FctNoCoeff* \[Beta]CoeffList[[regIdx]],
                      regList[[regIdx]] <= \[Mu]Renorm <= regList[[regIdx + 1]] }
                     		];
                   ];

                  Return[Piecewise[pieceWiseList]];
                  ];


                  addToBetaCoeffSinThW[currBetaCoef_,  GSMChargeDict_,
                     dofToAdd_, \[Mu]Fake_] := Block[
                     {sin2\[Theta]WFct},
                     chargeFactor = GSMChargeDict[["U1EM"]][[ "QEM"]];
                     colorFactor = GSMChargeDict[["SU3C"]];
                     isospinFactor =  GSMChargeDict[["U1EM"]][[ "T3"]];

                     If[dofToAdd == "Fermion",
                      (*All the fremions are chiral*)
                      gammaFactor = 8;
                      toAdd =
                       colorFactor * gammaFactor *
                        chargeFactor *(isospinFactor -
                          chargeFactor*sin2\[Theta]WFct[\[Mu]Fake])*nbOfFamilies;
                      ,

                      If[ dofToAdd == "Boson",
                        gammaFactor = -22;
                        toAdd =
                         colorFactor * gammaFactor *
                          chargeFactor *(isospinFactor -
                            chargeFactor*sin2\[Theta]WFct[\[Mu]Fake]);

                        ];
                      ];


                     (* Printing here. *)
                     (* Print["sin2ThW:",currBetaCoef + toAdd]; *)
                     Return[currBetaCoef + toAdd];

                     ];


                  getWeinbergSolvedRGE[solved\[Alpha]EM_, classDict_ , initSinThW_,
                    initBetaCoeff_, \[Mu]Renorm_] := Block[
                    {orderedDict, dictKeys, i, branchNb, startReg, endReg, KKNumber,
                     KKName, KKNameStripped, KKMass, \[Mu]0},
                    orderedDict = makeOrderedTowerDict[classDict];

                    dictKeys = Keys[orderedDict];
                    sinThWList = {initSinThW};
                    betaCoeff = initBetaCoeff;

                    pieceList = {};
                    sinThWList = {initSinThW};


                    For[branchNb = 1, branchNb <= Length[solved\[Alpha]EM[[1]]],
                     branchNb++,
                     (* Select the \[Mu] start and end caps,
                     and solve for the currentBC *)

                     startReg = solved\[Alpha]EM[[1, branchNb, 2, 1]] /k;
                     endReg = solved\[Alpha]EM[[1, branchNb, 2, 3]] / k;


                     For[KKNumber = 1, KKNumber <= Length[dictKeys], KKNumber++,

                      (* Get the particle Names and masses *)

                      KKName = dictKeys[[KKNumber]];
                      KKNameStripped =
                       StringTrim[dictKeys[[KKNumber]], RegularExpression["\\d+$"]];
                      KKMass = \[Lambda] /. orderedDict[[KKName]];


                      If[
                       startReg <= KKMass < endReg,

                       GSMChargeDict = classDict[[KKNameStripped]][["GSMCharges"]];
                       dofType =  classDict[[KKNameStripped]][["Type"]];

                       newBetaCoeff =
                        addToBetaCoeffSinThW[0, GSMChargeDict, dofType, \[Mu]Renorm];
                       betaCoeff = Simplify[betaCoeff + newBetaCoeff];
                       (*Print[betaCoeff];
                       Print[KKName]*)
                       (*Print[
                       solved\[Alpha]EM/.{\[Mu]Renorm\[Rule] KKMass * k} ];*)
                       ]

                      ];

                     \[Mu]0 = startReg * k;
                     branchBeta = betaCoeff /. {\[Mu]Renorm -> \[Mu]0};
                     (*Print[branchBeta];*)
                     sin2Val = sinThWList[[branchNb]];
                     newBranch =
                      sin2\[Theta]WFct[\[Mu]0]*(1 +
                          1/(24 \[Pi]) solved\[Alpha]EM[[1, branchNb, 1]] * 1 /
                           sin2\[Theta]WFct[\[Mu]0] * Log[\[Mu]0 / \[Mu]Renorm]*
                           branchBeta ) /. {sin2\[Theta]WFct[\[Mu]0] -> sin2Val};

                     AppendTo[
                      pieceList, {newBranch, startReg*k <= \[Mu]Renorm <= endReg*k}];
                     sin2Piece = Piecewise[pieceList];

                     (*Print[sin2Piece];*)
                     (*Print[startReg*k,"  ", endReg*k];
                     Print[newBranch];*)
                     (*Print[newBranch/.{\[Mu]\[Rule] endReg*
                     k}];*)
                     AppendTo[sinThWList, newBranch /. {\[Mu] -> endReg*k}];

                     ];



                    Return[sin2Piece];
                    ]
                    solvePieceWiseRGE[pieceWiseRGE_, startBoundaryCond_, gFct_, \[Mu]_] :=
                      Module[
                      {startReg, endReg, branchSol, branchNb, BCtoSolve, solList},

                      BCtoSolve = startBoundaryCond;
                      solList = {};

                      For[branchNb = 1, branchNb <= Length[pieceWiseRGE[[1]]],
                       branchNb++,

                       (*Select the \[Mu] start and end caps, and solve for the currentBC*)

                          startReg = pieceWiseRGE[[1, branchNb, 2, 1]];
                       endReg = pieceWiseRGE[[1, branchNb, 2, 3]];
                       branchSol = NDSolve[{\[Mu] \!\(
                    \*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\(gFct[\[Mu]]\)\) ==
                           pieceWiseRGE[[1, branchNb, 1]],
                          gFct[startReg] == BCtoSolve}, {gFct}, {\[Mu], startReg, endReg}];
                       (*The new BC corresponds to the value of the current solution at \
                    the en point of the interval*)

                       BCtoSolve = gFct[\[Mu]] /. branchSol[[1]] /. {\[Mu] -> endReg};

                       AppendTo[
                        solList, {gFct[\[Mu]] /. branchSol[[1, 1]],
                         startReg <= \[Mu] <= endReg}];

                       ];


                      Return[Piecewise[solList]];
                      ]

                    plotRGEsol[solvedRGEList_, orderedTowerDict_, \[Mu]var_, startReg_,
                      plotLegend_, plotRange_: Automatic, sin2WEnable_: False,
                      plotStyle_: {}] := Module[
                      {epilogList, lineWithText, newReg, towerKeys, tickList, i,
                       frameLabel},


                      epilogList = {};
                      tickList = {};
                      towerKeys = Keys[orderedTowerDict];

                      For[i = 1, i <= Length[towerKeys], i++,

                       newReg = k*\[Lambda] /. orderedTowerDict[[ towerKeys[[i]] ]];
                       lineWithText =  InfiniteLine[{{newReg, 0.0}, {newReg, 0.5}}];
                       AppendTo[epilogList, lineWithText];
                       AppendTo[tickList, {newReg, towerKeys[[i]], {0.5, 0}} ]

                       ];

                      If[sin2WEnable == True,
                       frameLabel =
                         MaTeX[{"\\mu (GeV)", " \\sin \\theta_W"}, FontSize -> 24];,
                       frameLabel = MaTeX[{"\\mu (GeV)", " \\alpha^{-1}"}, FontSize -> 24];
                       ];

                      texStyle = {FontFamily -> "Latin Modern Roman", FontSize -> 13};
                      Plot[solvedRGEList,  {\[Mu]var, startReg, mKK5},

                       PlotRange -> plotRange,
                       Epilog -> {Dashed, Thickness[0.001], epilogList},
                       (*AxesLabel\[Rule]{"\[Mu] (GeV)", "(4\[Pi])/\[Alpha]"},*)

                       Frame -> True, FrameStyle -> BlackFrame,
                       FrameLabel -> frameLabel,

                       BaseStyle -> texStyle,
                       ImageSize -> Large,
                       PlotStyle -> plotStyle,

                       PlotLegends -> MaTeX[plotLegend, FontSize -> 18]
                       (*{"\\frac{4 \\pi}{g^2_{3C}}", "\\frac{4 \\pi}{g^2_{2L}}",
                       "\\frac{4 \\pi}{g^2_{3C}}"}*)
                       (*,PlotRange\[Rule]{Automatic,{0,
                       100}}*)
                       (*,Ticks\[Rule]{tickList, Automatic}*)

                       ]



                      ]
(*------------------- Analysis run --------------------*)

(*Put all the solutions in a list , and make a oredered dictionary \
with their names and their \[Lambda] values*)
(*fullList = Sort[Join[\
\[Lambda]ListBottom,\[Lambda]ListTau, \[Lambda]ListPsiDark]];*)

(*********************************************************************************************************************)
classDict = <|
   "Tau" ->       <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListTau, 1] ,

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 2,
       "U1EM" ->  <|"QEM" -> -1, "T3" -> -1/2|> |>
     				     |> ,

   "Bottom" -> <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListBottom, 1],

     "GSMCharges" -> <|"SU3C" -> 3, "SU2L" -> 2,
       "U1EM" ->  <|"QEM" -> -1/3, "T3" -> -1/2|>  |>
     				      |>,

   "Top" ->         <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListTop, 1] ,

     "GSMCharges" -> <|"SU3C" -> 3, "SU2L" -> 1,
       "U1EM" ->  <|"QEM" -> 2/3, "T3" -> +1/2|>  |>
     				       |>,
   (*******************NEED TO CHECK THE CHARGE ASSIGNMENT/
   S FOR MPSIDARK*************)

   "mPsiDarkUno" ->         <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListPsiDark, 1] ,

     "GSMCharges" -> <|"SU3C" -> 3, "SU2L" -> 1,
       "U1EM" ->  <|"QEM" -> 2/3, "T3" -> +1/2|>  |>
     				       |>,
   "mPsiDarkDos" ->         <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListPsiDark, 1] ,

     "GSMCharges" -> <|"SU3C" -> 3, "SU2L" -> 1,
       "U1EM" ->  <|"QEM" -> -2/3, "T3" -> -1/2|>  |>
     				       |>,
   "mPsiDarkTres" ->         <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListPsiDark, 1] ,

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 1,
       "U1EM" ->  <|"QEM" -> -1, "T3" -> -1/2|>  |>
     				       |>,
   "mPsiDarkQuatro" ->         <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]ListPsiDark, 1] ,

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 1,
       "U1EM" ->  <|"QEM" -> 0, "T3" -> +1/2|>  |>
     				       |>,
   (*******************NEED TO CHECK THE CHARGE ASSIGNMENT/
   S FOR MPSIDARK*************)


   		     (*"W-"\[Rule]        <|"Type"\[Rule]"Boson",
   "Sols"\[Rule]Delete[\[Lambda]ListW,1],
   					 "GSMCharges"\[Rule]<|"SU3C"\[Rule]1, "U1EM"\[Rule]  -1 |>
   					|>,*)

   "Neutrino" -> <|"Type" -> "Fermion",
     "Sols" -> Delete[\[Lambda]List\[Nu], 1],

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 2,
       "U1EM" -> <|"QEM" -> 0,  "T3" -> +1/2|> |>
     					|>,
   		    "W-" ->        <|"Type" -> "Boson",
     "Sols" -> Delete[\[Lambda]ListW, 1],

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 3,
       "U1EM" ->  <|"QEM" -> -1, "T3" -> -1|>  |>
     					|>,
   "W+" ->        <|"Type" -> "Boson",
     "Sols" -> Delete[\[Lambda]ListW, 1],

     "GSMCharges" -> <|"SU3C" -> 1, "SU2L" -> 3,
       "U1EM" ->  <|"QEM" -> +1, "T3" -> +1|>  |>
     					|>
   |>;
(*********************************************************************************************************************)


dictMerged = <||>;

allPartNames = Keys[classDict];

For [i = 1, i <= Length[allPartNames], i++,
  partName = allPartNames[[i]];
  newDict = makeDictTowers[classDict[[partName]][["Sols"]], partName ];
  AppendTo[dictMerged, newDict];
  ];


orderedTowerDict = SortBy[dictMerged, # &];


couplingDict = <|
   "SU2L" -> <|"Coupling" -> g2L, "InitBeta" -> -19/16,
     "InitBC" -> Sqrt[4 \[Pi] /127.91]/Sqrt[sin2\[Theta]W]|>,
   			   "SU3C" -> <|"Coupling" -> g3C, "InitBeta" -> -7,
     "InitBC" -> Sqrt[4 \[Pi] 0.11822]|>
   |>;

solvedRGEdict = <||>;
solvedRGEList = {};
gaugeGroupList = Keys[couplingDict];

For[iGauge = 1, iGauge <= Length[gaugeGroupList], iGauge++,
 gaugeGroup = gaugeGroupList[[iGauge]];
 gaugeCoupling = couplingDict[[gaugeGroup]][["Coupling"]];
 gaugeInitBeta = couplingDict[[gaugeGroup]][["InitBeta"]];
 gaugeInitBC = couplingDict[[gaugeGroup]][["InitBC"]];

 rgeSystem =
  makePiecewiseRGE[classDict, gaugeCoupling, gaugeInitBeta,
   gaugeGroup,  \[Mu], k];
 solvedRGE =
  solvePieceWiseRGE[rgeSystem, gaugeInitBC, gaugeCoupling, \[Mu]];

 solvedRGEdict[[ gaugeGroup]] = solvedRGE;
 AppendTo[solvedRGEList, solvedRGE];

 ]

(************************************** RGE SOLVING **********************************************)
rgeSystem =  makePiecewiseRGE[classDict, \[Alpha]1EM, initBetaEM, "U1EM",  \[Mu], k];
solvedRGE\[Alpha]EM =  solvePieceWiseRGE[rgeSystem,   Sqrt[4 \[Pi] /127.91] , \[Alpha]1EM, \[Mu]];
sin2WSolved =  getWeinbergSolvedRGE[solvedRGE\[Alpha]EM, classDict, sin2\[Theta]W, initBetaSin2ThW, \[Mu]];

(************************************************************************************)

\[Alpha]1YatMZinv = 3/5 (1 - sin2WSolved) (4 \[Pi] / solvedRGE\[Alpha]EM^2) /. {\[Mu] ->     MZ};
\[Alpha]2LatMZinv = sin2WSolved * (4 \[Pi] / solvedRGE\[Alpha]EM^2) /. {\[Mu] -> MZ};
\[Alpha]3CatMZinv = 4 \[Pi] / solvedRGEdict[["SU3C"]]^2 /. {\[Mu] -> MZ};

Print["-------------------------------------"];
Print["At MZ, (\[Alpha]1Y)^-1 :" ,\[Alpha]1YatMZinv];
Print["At MZ, (\[Alpha]2L)^-1 :" ,\[Alpha]2LatMZinv];
Print["At MZ, (\[Alpha]3C)^-1 :" ,\[Alpha]3CatMZinv];
Print["-------------------------------------"];


(**************************************1**********************************************)

\[Alpha]1YatMKK5inv = 3/5 (1 - sin2WSolved) (4 \[Pi] / solvedRGE\[Alpha]EM^2) /. {\[Mu] ->     mKK5};
\[Alpha]2LatMKK5inv = sin2WSolved * (4 \[Pi] / solvedRGE\[Alpha]EM^2) /. {\[Mu] -> mKK5};
\[Alpha]3CatMKK5inv = 4 \[Pi] / solvedRGEdict[["SU3C"]]^2 /. {\[Mu] -> mKK5};

\[Alpha]4inv = \[Alpha]3CatMKK5inv + 1/(12 \[Pi])
\[Alpha]2Linv = \[Alpha]2LatMKK5inv
\[Alpha]2Rinv = 5/3 \[Alpha]1YatMKK5inv - 2/3 \[Alpha]3CatMKK5inv + 8/(45 \[Pi])

(* Print["Found the boundary conditions."] *)
Print["-------------------------------------"];
Print["At MKK5, (\[Alpha]1Y)^-1 :" ,\[Alpha]1YatMKK5inv];
Print["At MKK5, (\[Alpha]2L)^-1 :" ,\[Alpha]2LatMKK5inv];
Print["At MKK5, (\[Alpha]3C)^-1 :" ,\[Alpha]3CatMKK5inv];
Print["-------------------------------------"];
Print["-------------------------------------"];
Print["At MKK5, (\[Alpha]4C)^-1 :" ,\[Alpha]4inv];
Print["At MKK5, (\[Alpha]2L)^-1 :" ,\[Alpha]2Linv];
Print["At MKK5, (\[Alpha]2R)^-1 :" ,\[Alpha]2Rinv];
(* Print[\[Alpha]4inv, "  ", \[Alpha]2Linv, " ", \[Alpha]2Rinv]; *)
Print["-------------------------------------"];
If[makePlots == True,
  Print["--> Plotting 4D  RGES"];

  pltsin2ThW = plotRGEsol[sin2WSolved, orderedTowerDict, \[Mu], 80 , \
                          {"\\sin\\theta_W"}, Automatic, True, ColorData[97, "ColorList"][[8]]];

  pltg3C = plotRGEsol[ 4 \[Pi] /    solvedRGEdict[["SU3C"]]^2, orderedTowerDict, \[Mu], 80 , {"\\frac{4 \
                    \\pi}{g^2_{3C}}"}, Automatic, False, {ColorData[97, "ColorList"][[3]]}
                   ];

  pltgEM = plotRGEsol[{4 \[Pi] /   solvedRGE\[Alpha]EM^2}, orderedTowerDict, \[Mu], 80 , {"\\frac{4 \
                    \\pi}{g^2_{EM}}"} , Automatic, False, ColorData[97, "ColorList"][[7]]
                   ];

  plt2L = plotRGEsol[ sin2WSolved * (4 \[Pi] /   solvedRGE\[Alpha]EM^2), orderedTowerDict, \[Mu], 80 , {"\\frac{4 \
                    \\pi}{g^2_{2L}}"} , Automatic, False, {Orange}
                   ];
  plt1Y = plotRGEsol[ 3/5 (1 - sin2WSolved) (4 \[Pi] /  solvedRGE\[Alpha]EM^2), orderedTowerDict, \[Mu], 80 , {"\\frac{4 \
                    \\pi}{g^2_{1Y}}"}
                   ];
  pltGSM = plotRGEsol[{3/5 (1 - sin2WSolved) (4 \[Pi] / solvedRGE\[Alpha]EM^2),
                      sin2WSolved * (4 \[Pi] / solvedRGE\[Alpha]EM^2),
                      4 \[Pi] /  solvedRGEdict[["SU3C"]]^2}, orderedTowerDict, \[Mu], 80 , {"\\frac{4 \
                              \\pi}{g^2_{1Y}}", "\\frac{4 \\pi}{g^2_{2L}}", "\\frac{4 \\pi}{g^2_{3C}}" } , Automatic
                   ];


  Export["WeinbergPlots/sin2ThW_4D.pdf", pltsin2ThW];
  Export["WeinbergPlots/g3C_4D.pdf", pltg3C];
  Export["WeinbergPlots/gEM_4D.pdf", pltgEM];
  Export["WeinbergPlots/g2L_4D.pdf", plt2L];
  Export["WeinbergPlots/g1Y_4D.pdf", plt1Y];
  Export["WeinbergPlots/gGSM_4D.pdf", pltGSM];

]

Print["------------ Starting 5D Analysis ---------------"]
(*************============================================================================================*******)
(*************============================================================================================*******)
(*************-------------------------5D ANALYSIS--------------------------------------------------------*******)
(*************============================================================================================*******)
(*************============================================================================================*******)

(* Subscript[N, ++] *)

Npp[\[Mu]_, s\[Xi]_, r0\[Xi]_, r\[Pi]\[Xi]_, \[Alpha]\[Xi]_] :=
  Module[ {\[Mu]Aux},
   Return[
     (( s\[Xi]/2 - r0\[Xi] ) BesselJ[\[Alpha]\[Xi], \[Mu] /
             k ]  + \[Mu]/
           k D[BesselJ[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/
           k}) * (  (
            s\[Xi]/2 -
             r\[Pi]\[Xi] ) BesselY[\[Alpha]\[Xi], \[Mu] / (k zL) ]  + \
\[Mu]/(k zL)
            D[BesselY[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/(
           k zL)}) + (  (
            s\[Xi]/2 -
             r\[Pi]\[Xi] ) BesselJ[\[Alpha]\[Xi], \[Mu] / (k zL) ]  + \
\[Mu]/(k zL)
            D[BesselJ[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/(
           k zL)}) *(  (
            s\[Xi]/2 - r0\[Xi] ) BesselY[\[Alpha]\[Xi], \[Mu] /
             k ]  + \[Mu]/
           k D[BesselY[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/k}) ];
   ];

(*Subscript[N, +-]*)

Npm[\[Mu]_, s\[Xi]_, r0\[Xi]_, r\[Pi]\[Xi]_, \[Alpha]\[Xi]_] :=
  Module[{\[Mu]Aux},
   Return[
     -  BesselY[\[Alpha]\[Xi], \[Mu] / (k zL) ] *(  (

          s\[Xi]/2 - r0\[Xi] ) BesselJ[\[Alpha]\[Xi], \[Mu] /
             k ]  + \[Mu]/
           k D[BesselJ[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/k})  +
      BesselJ[\[Alpha]\[Xi], \[Mu] / (k zL) ] * (  (
            s\[Xi]/2 - r0\[Xi] ) BesselY[\[Alpha]\[Xi], \[Mu] /
             k ]  + \[Mu]/
           k D[BesselY[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/k})];
   ];

(*Subscript[N, -+]*)

Nmp[\[Mu]_, s\[Xi]_, r0\[Xi]_, r\[Pi]\[Xi]_, \[Alpha]\[Xi]_] :=
  Module[{\[Mu]Aux},
   Return[
     BesselJ[\[Alpha]\[Xi], \[Mu] /
         k ] * (  (
            s\[Xi]/2 -
             r\[Pi]\[Xi] ) BesselY[\[Alpha]\[Xi], \[Mu] / (k zL) ]  + \
\[Mu]/(k zL)
            D[BesselY[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/(k zL)})  -
       BesselY[\[Alpha]\[Xi], \[Mu] /
         k ]*(  (
            s\[Xi]/2 -
             r\[Pi]\[Xi] ) BesselJ[\[Alpha]\[Xi], \[Mu] / (k zL) ]  + \
\[Mu]/(k zL)
            D[BesselJ[\[Alpha]\[Xi], \[Mu]Aux], \[Mu]Aux] /. \
{\[Mu]Aux -> \[Mu]/(k zL)}) ];
   ];


(*Subscript[N, --]*)

Nmm[\[Mu]_, s\[Xi]_, r0\[Xi]_, r\[Pi]\[Xi]_, \[Alpha]\[Xi]_] :=
  Module[{\[Mu]Aux},
   Return[
     BesselJ[\[Alpha]\[Xi], \[Mu] / k ] *
       BesselY[\[Alpha]\[Xi], \[Mu] / (k zL) ]  -
      BesselJ[\[Alpha]\[Xi], \[Mu] / (k zL) ]  *
       BesselY[\[Alpha]\[Xi], \[Mu] / k ]
      ];
   ];

(*N functions that take into account the \[PlusMinus]\[PlusMinus] \
assignment and the field type to assign the right r0, r\[Pi], \
\[Alpha], s, *)

NFor\[Xi][\[Mu]Var_, \[Xi]Type_, M\[Xi]_, sign1_, sign2_] := Module[
   {s\[Xi], \[Alpha]\[Xi], r0\[Xi], r\[Pi]\[Xi]},

   If[\[Xi]Type == "Boson",
    If[sign1 == 1 && sign2 == 1,
     s\[Xi] = 4; \[Alpha]\[Xi] = 0; r0\[Xi] = r\[Pi]\[Xi] = 2;
     Return[
      Nmm[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
     ,
     s\[Xi] = 2; \[Alpha]\[Xi] = 1; r0\[Xi] = r\[Pi]\[Xi] = 0;]

    ,
    If[\[Xi]Type == "FermionL",
     If[sign1 == 1 && sign2 == 1,
      s\[Xi] = 1;
      r0\[Xi] = r\[Pi]\[Xi] = M\[Xi]; \[Alpha]\[Xi] =
       Abs[M\[Xi] - 1/2];
      Return[
       Nmm[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
      ,
      s\[Xi] = 1;
      r0\[Xi] = r\[Pi]\[Xi] = -M\[Xi]; \[Alpha]\[Xi] =
       Abs[M\[Xi] + 1/2];]

     ,
     If[\[Xi]Type == "FermionR",
      If[sign1 == 1 && sign2 == 1,
       s\[Xi] = 1;
       r0\[Xi] = r\[Pi]\[Xi] = -M\[Xi]; \[Alpha]\[Xi] =
        Abs[M\[Xi] + 1/2];
       Return[
        Nmm[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
       ,
       s\[Xi] = 1;
       r0\[Xi] = r\[Pi]\[Xi] = M\[Xi]; \[Alpha]\[Xi] =
        Abs[M\[Xi] - 1/2];]

      ,
      If[\[Xi]Type == "Scalar",
       s\[Xi] = 4; r0\[Xi] = r\[Pi]\[Xi] = 0; \[Alpha]\[Xi] = 2;
       If[sign1 == 1 && sign2 == 1,
        Return[
          Npp[\[Mu]Var, s\[Xi], r0\[Xi],
           r\[Pi]\[Xi], \[Alpha]\[Xi]]];  ]]
      ]]];


   If[sign1 == 1 && sign2 == -1,
    Return[Npm[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
    ,
    If[sign1 == -1 && sign2 == 1,
     Return[
       Nmp[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
     ,
     If[ sign1 == -1 && sign2 == -1,
      Return[
        Nmm[\[Mu]Var, s\[Xi], r0\[Xi], r\[Pi]\[Xi], \[Alpha]\[Xi]]];
      ]]];
   ];
   makeDeltaForSUNscalar[n_, sign1_, sign2_, \[Mu]_, M\[Xi]_,
     dR_, \[CapitalLambda]Cutoff_] := Module[
     {Nfunc, \[CapitalDelta]a\[Phi], irrepDictSUN, Ta\[Phi], Nvar,
      nIntFct, numIntRes, u},



     (* Performing the N function numerical integration *)
     (*
     NOTE! You should only keep the real part of the integration for \
   (-,-) since the imaginary part is super oscillatory and should be 0 \
   in the theory. All other cases have Im =
     0 See commented plotting code for clarification  *)

     Nvar = (I u)/2* Sqrt[\[Mu]^2];


     Nfunc = NFor\[Xi][Nvar, "Scalar", M\[Xi], sign1, sign2];
     nIntFct = u (1 - u^2)^(1/2)* Log[Nfunc];
     (*numIntRes = NIntegrate[Re @ nIntFct,{u, 0.0, 1.0}];*)


     numIntRes =
      Quiet[ NIntegrate[Re @ nIntFct, {u, machineZero, 1.0}]];
     (*Plot[{
     Re @ N[Log[NFor\[Xi][\[ImaginaryI] x, "Scalar",0,pm1,
     pm2]/.{x\[Rule]50*uval}]] * uval(1-uval^2)^(1/2),
     Im @ N[Log[NFor\[Xi][\[ImaginaryI] x, "Scalar",0, pm1,
     pm2]/.{x\[Rule]50*uval}]] * uval(1-uval^2)^(1/2)},
     {uval, 0.0, 1.0}, PlotLegends\[Rule]{"Real", "Imaginary"}]*)


     (*Make the gauge structure dictionary*)

     irrepDictSUN = makeSUNassocDict[n];
     Ta\[Phi] = irrepDictSUN[[ ToString[dR] ]][["DynkinIdx"]];

     (*Construct and return the \[CapitalDelta] corrections*)

     If[sign1 == 1 && sign2 == 1,
      \[CapitalDelta]a\[Phi] =
        1/12 * Ta\[Phi] * (Log[\[CapitalLambda]Cutoff/k] -
           3*numIntRes  )  ;
      ,
      If[(sign1 == 1 && sign2 == -1) || (sign1 == -1 && sign2 == +1),
       \[CapitalDelta]a\[Phi] = 1/12*((-3)* Ta\[Phi] * numIntRes);
       ,
       If[sign1 == -1 && sign2 == -1,
        \[CapitalDelta]a\[Phi] = -(1/12) *
           Ta\[Phi] * (Log[\[CapitalLambda]Cutoff/k] +
             3*numIntRes  )  ;

        ]]];


     Return[\[CapitalDelta]a\[Phi]];]

   (******************  Fermions  *********************)

   makeDeltaForSUNFermion[n_, sign1_, sign2_, \[Mu]_, M\[Xi]_, dR_, L5_,
      FermionStr_] := Module[
      {Nfunc, \[CapitalDelta]a\[Psi], irrepDictSUN, Ta\[Psi], Nvar, u,
       nIntFct, numIntRes},


      (* Performing the N function numerical integration *)

      Nvar = (I u)/2* Sqrt[\[Mu]^2];

      If[sign1 == 1 && sign2 == 1 &&  FermionStr == "FermionL",
       Nfunc = NFor\[Xi][Nvar, "FermionR", M\[Xi], sign1, sign2];,

       If[sign1 == 1 && sign2 == 1 &&  FermionStr == "FermionR",
        Nfunc = NFor\[Xi][Nvar, "FermionL", M\[Xi], sign1, sign2];,
        Nfunc = NFor\[Xi][Nvar, FermionStr, M\[Xi], sign1, sign2];

        ]];



      nIntFct = (  u (1 - u^2)^(1/2) - u (1 - u^2)^(-1/2) )* Log[Nfunc];
      numIntRes = Quiet[ NIntegrate[Re @ nIntFct, {u, machineZero, 1.0}]];

      (*Make the gauge structure dictionary*)

      irrepDictSUN = makeSUNassocDict[n];
      Ta\[Psi] = irrepDictSUN[[ ToString[dR] ]][["DynkinIdx"]];


      If[ (sign1 == 1 && sign2 == 1) || (sign1 == -1 && sign2 == -1),
       \[CapitalDelta]a\[Psi] =
         1/3 * Ta\[Psi] * (2 Log[k/\[Mu] ] - k L5 +  3*numIntRes  ) ;
       ,

       If[sign1 == 1 && sign2 == -1,
        \[CapitalDelta]a\[Psi] = 1/3 * Ta\[Psi] *(-k L5 + 3 numIntRes);
        ,

        If[sign1 == -1 && sign2 == 1,
         \[CapitalDelta]a\[Psi] =
           1/3 * Ta\[Psi] * (k L5 + 3 numIntRes);]]];

      Return[\[CapitalDelta]a\[Psi]];];


   (******************  Vector bosons  *********************)

   makeDeltaForSUNVector[n_, sign1_, sign2_, \[Mu]_, M\[Xi]_, dR_,
      L5_, \[CapitalLambda]_] := Module[
      {Nfunc, \[CapitalDelta]aA, irrepDictSUN, TaA, Nvar, u, nIntFct,
       numIntRes},


      (* Performing the N function numerical integration *)

      Nvar = (I u)/2* Sqrt[\[Mu]^2];
      Nfunc = NFor\[Xi][Nvar, "Boson", M\[Xi], sign1, sign2];
      nIntFct = (  -9 u (1 - u^2)^(1/2) + 24 u (1 - u^2)^(-1/2) )*
        Log[Nfunc];
      numIntRes = Quiet[ NIntegrate[Re @ nIntFct, {u, machineZero, 1.0}]];

      (*Make the gauge structure dictionary*)

      irrepDictSUN = makeSUNassocDict[n];
      TaA = irrepDictSUN[[ ToString[dR] ]][["DynkinIdx"]];


      If[ sign1 == 1 && sign2 == 1,
       \[CapitalDelta]aA =
         1/12*TaA*(23 Log[\[Mu]/ \[CapitalLambda]] + 21 Log[\[Mu] /k] +
            22 k L5 + numIntRes);
       ,

       If[ sign1 == +1 && sign2 == -1,
        \[CapitalDelta]aA = 1/12* TaA * (-k L5 + numIntRes);
        ,

        If[sign1 == -1 && sign2 == +1,
         \[CapitalDelta]aA = 1/12* TaA * (k L5 + numIntRes);
         ,

         If[ sign1 == -1 && sign2 == -1,
          \[CapitalDelta]aA =
            1/12*TaA*(23 Log[\[Mu]/ \[CapitalLambda]] +
               2 Log[\[Mu] /k] - k L5 + numIntRes);
          ]


         ]]];
      Return[\[CapitalDelta]aA];
      ];
      (*Clear[]*)

      SU4Corr[\[Mu]Val_, \[CapitalLambda]Cut_, L5_, c0_, c0Prime_, c1_,
        c2_] := Module[
        {nbOfQuarkGens, SUN, SU4ScalarCorr, SU4FermionCorr, SU4VectorCorr,
         fullSU4Corr},

        SUN = 4;

        (* Scalar contribution*)

        SU4ScalarCorr =
         2*makeDeltaForSUNscalar[SUN, -1, -1, \[Mu]Val, 0,
           15, \[CapitalLambda]Cut];

        (*Fermion Contributions*)
        nbOfQuarkGens = nbOfFamilies;

        SU4FermionCorr =
         nbOfQuarkGens*(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c0, 4,
               L5, "FermionL"] +
             makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c0, 4, L5,
              "FermionR"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c0, 4, L5,
              "FermionL"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c0, 4, L5,
              "FermionR"]) +
          nbOfQuarkGens*(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c1,
              6, L5, "FermionL"] +
             makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c1, 6, L5,
              "FermionR"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c1, 6, L5,
              "FermionL"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c1, 6, L5,
              "FermionR"]) +
          makeDeltaForSUNFermion[SUN, +1, -1, \[Mu]Val, c0Prime, 4, L5,
           "FermionL"] +
          makeDeltaForSUNFermion[SUN, -1, +1, \[Mu]Val, c0Prime, 4, L5,
           "FermionR"] +
          makeDeltaForSUNFermion[SUN, -1, +1, \[Mu]Val, c0Prime, 4, L5,
           "FermionL"] +
          makeDeltaForSUNFermion[SUN, +1, -1, \[Mu]Val, c0Prime, 4, L5,
           "FermionR"];

        (*Note that 3bar is identical to 3 when our concern is RGE runnings \
      since it only depends on the Dynkin index which is the same*)

        SU4VectorCorr =
         makeDeltaForSUNVector[SUN - 1, +1, +1, \[Mu]Val, 0, 8,
           L5, \[CapitalLambda]Cut] +
          makeDeltaForSUNVector[SUN - 1, -1, +1, \[Mu]Val, 0, 3,
           L5, \[CapitalLambda]Cut] +
          makeDeltaForSUNVector[SUN - 1, -1, +1, \[Mu]Val, 0, 3,
           L5, \[CapitalLambda]Cut];

        fullSU4Corr = SU4ScalarCorr + SU4FermionCorr + SU4VectorCorr;
        Return[1/(8 \[Pi]^2) fullSU4Corr];
        ]

      SU2LCorr[\[Mu]Val_, \[CapitalLambda]Cut_, L5_, c0_, c0Prime_, c1_,
        c2_] := Module[
        {nbOfQuarkGens, SUN, SU2LScalarCorr, SU2LFermionCorr,
         SU2LVectorCorr, fullSU2LCorr},

        SUN = 2;

        (* Scalar contribution*)

        SU2LScalarCorr =
         makeDeltaForSUNscalar[SUN, +1, +1, \[Mu]Val, 0,
           2, \[CapitalLambda]Cut] +
          2* makeDeltaForSUNscalar[SUN, -1, -1, \[Mu]Val, 0,
            3, \[CapitalLambda]Cut] ;

        (*Fermion Contributions*)
        nbOfQuarkGens = nbOfFamilies;
        SU2LFermionCorr =
         nbOfQuarkGens*(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c0, 2,
               L5, "FermionL"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c0, 2, L5,
              "FermionR"]) +
          nbOfQuarkGens *(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c2,
              2, L5, "FermionL"] +

             makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c2, 2, L5,
              "FermionR"] +

             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c2, 2, L5,
              "FermionL"] +

             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c2, 2, L5,
              "FermionR"]
            ) +
          	makeDeltaForSUNFermion[SUN, +1, -1, \[Mu]Val, c0Prime, 2, L5,
           "FermionL"] +
          makeDeltaForSUNFermion[SUN, -1, +1, \[Mu]Val, c0Prime, 2, L5,
           "FermionR"];

        SU2LVectorCorr =
         makeDeltaForSUNVector[SUN, +1, +1, \[Mu]Val, 0, 3,
           L5, \[CapitalLambda]Cut] +

          makeDeltaForSUNVector[SUN, -1, -1, \[Mu]Val, 0, 2,
           L5, \[CapitalLambda]Cut]
        ;

        fullSU2LCorr = SU2LScalarCorr + SU2LFermionCorr + SU2LVectorCorr;
        Return[1/(8 \[Pi]^2) fullSU2LCorr];

        ]


      SU2RCorr[\[Mu]Val_, \[CapitalLambda]Cut_, L5_, c0_, c0Prime_, c1_,
        c2_] := Module[
        {nbOfQuarkGens, SUN, SU2RScalarCorr, SU2RFermionCorr,
         SU2RVectorCorr, fullSU2RCorr},

        SUN = 2;

        (* Scalar contribution*)

        SU2RScalarCorr =
         makeDeltaForSUNscalar[SUN, +1, +1, \[Mu]Val, 0,
           2, \[CapitalLambda]Cut] +
          2* makeDeltaForSUNscalar[SUN, -1, -1, \[Mu]Val, 0,
            3, \[CapitalLambda]Cut] ;

        (*Fermion Contributions*)
        nbOfQuarkGens = nbOfFamilies;
        SU2RFermionCorr =
         nbOfQuarkGens*(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c0, 2,
               L5, "FermionL"] +
             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c0, 2, L5,
              "FermionR"]) +
          nbOfQuarkGens *(makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c2,
              2, L5, "FermionL"] +

             makeDeltaForSUNFermion[SUN, +1, +1, \[Mu]Val, c2, 2, L5,
              "FermionR"] +

             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c2, 2, L5,
              "FermionL"] +

             makeDeltaForSUNFermion[SUN, -1, -1, \[Mu]Val, c2, 2, L5,
              "FermionR"]
            ) +
          	makeDeltaForSUNFermion[SUN, -1, +1, \[Mu]Val, c0Prime, 2, L5,
           "FermionL"] +

          makeDeltaForSUNFermion[SUN, +1, -1, \[Mu]Val, c0Prime, 2, L5,
           "FermionR"];

        SU2RVectorCorr =
         makeDeltaForSUNVector[SUN, -1, +1, \[Mu]Val, 0, 3,
           L5, \[CapitalLambda]Cut] +

          makeDeltaForSUNVector[SUN, -1, -1, \[Mu]Val, 0, 2,
           L5, \[CapitalLambda]Cut]
        ;

        fullSU2RCorr = SU2RScalarCorr + SU2RFermionCorr + SU2RVectorCorr;
        Return[1/(8 \[Pi]^2) fullSU2RCorr];

        ]



gSU4 = Sqrt[4 \[Pi]  / \[Alpha]4inv];
gSU2L = Sqrt[4 \[Pi] /  \[Alpha]2Linv];
gSU2R = Sqrt[4 \[Pi] /  \[Alpha]2Rinv];

(*Determining \[CapitalLambda]Max for diferent couplings via the \
functional route.*)
(* Uncomment bellow for simple approximation *)
\
(*CtSU4Approx =\[Pi]/(2 gSU4) /.{\[Mu]Main\[Rule] mKK5};\
\[CapitalLambda]MaxApprox = CtSU4Approx * 16 \[Pi]  / L5;*)


CSU2RSol  =
FindRoot[1/gSU2R^2 -
   SU2RCorr[mKK5, \[CapitalLambda]Fake, L5, c0, c0Prime, c1, c2] -
   CSU2R /. {\[CapitalLambda]Fake ->
    16 \[Pi] / L5 * CSU2R}, {CSU2R, 100}];
\[CapitalLambda]MaxFunctional2R = CSU2R * 16 \[Pi]  / L5 /. CSU2RSol;
CSU2LSol  =
FindRoot[1/gSU2L^2 -
   SU2LCorr[mKK5, \[CapitalLambda]Fake, L5, c0, c0Prime, c1, c2] -
   CSU2L /. {\[CapitalLambda]Fake ->
    16 \[Pi] / L5 * CSU2L}, {CSU2L, 100}];
\[CapitalLambda]MaxFunctional2L = CSU2L * 16 \[Pi]  / L5 /. CSU2LSol;
CSU4Sol  =
FindRoot[1/gSU4^2 -
   SU4Corr[mKK5, \[CapitalLambda]Fake, L5, c0, c0Prime, c1, c2] -
   CSU4 /. {\[CapitalLambda]Fake -> 16 \[Pi] / L5 * CSU4}, {CSU4,
  100}];
\[CapitalLambda]MaxFunctional = CSU4 * 16 \[Pi]  / L5 /. CSU4Sol;



Print["Functional route \[CapitalLambda]Max : ", \
\[CapitalLambda]MaxFunctional];
(* ,
	"\nApproximation route for \[CapitalLambda]Max : ", \
\[CapitalLambda]MaxApprox]; *)


texStyle = {FontFamily -> "Latin Modern Roman", FontSize -> 13};

(**** Lambda Max Plot****)
(* \[CapitalLambda]Sol1 =
  TimeConstrained[
   Quiet[FindRoot[
     4 \[Pi] (CSU4  +
           SU4Corr[\[CapitalLambda]MaxSol, \
\[CapitalLambda]MaxFunctional, L5, c0, c0Prime, c1, c2] /. CSU4Sol) - nonPertBound ==
      0, {\[CapitalLambda]MaxSol, \[CapitalLambda]MaxFunctional}]],
   timeOutSolve];
Print["SU4 LambdaMax:", \[CapitalLambda]Sol1];


\[CapitalLambda]Sol2 =
  TimeConstrained[
   Quiet[FindRoot[
     4 \[Pi] (CSU2L  +
           SU2LCorr[\[CapitalLambda]MaxSol, \
\[CapitalLambda]MaxFunctional2L, L5, c0, c0Prime, c1, c2] /.
          CSU2LSol) - nonPertBound ==
      0, {\[CapitalLambda]MaxSol, \[CapitalLambda]MaxFunctional}]],
   timeOutSolve];
Print["SU2L LambdaMax:", \[CapitalLambda]Sol2];

\[CapitalLambda]Sol3 = TimeConstrained[
   Quiet[FindRoot[
     4 \[Pi] (CSU2R  +
           SU2RCorr[\[CapitalLambda]MaxSol, \
\[CapitalLambda]MaxFunctional2R, L5, c0, c0Prime, c1, c2] /.
          CSU2RSol) - nonPertBound ==
      0, {\[CapitalLambda]MaxSol, \[CapitalLambda]MaxFunctional}]],
   timeOutSolve];
Print["SU2R LambdaMax:", \[CapitalLambda]Sol3];

\[Lambda]sol1 = \[CapitalLambda]MaxSol /. \[CapitalLambda]Sol1;
\[Lambda]sol2 = \[CapitalLambda]MaxSol /. \[CapitalLambda]Sol2;
\[Lambda]sol3 = \[CapitalLambda]MaxSol /. \[CapitalLambda]Sol3;

\[Lambda]Sols = {\[Lambda]sol1, \[Lambda]sol2, \[Lambda]sol3};
\[Lambda]SolsFinal = {};
For[i = 1, i <= Length[\[Lambda]Sols], i++,
 If [ TrueQ[Head@\[Lambda]Sols[[i]] == Real],
  AppendTo[\[Lambda]SolsFinal, \[Lambda]Sols[[i]]];
  ]
 ]

If[TrueQ[Length[\[Lambda]SolsFinal] > 0],
 \[CapitalLambda]MaxPlot = Min[\[Lambda]SolsFinal];,
 \[CapitalLambda]MaxPlot =
   Min[\[CapitalLambda]MaxFunctional, \
\[CapitalLambda]MaxFunctional2L, \[CapitalLambda]MaxFunctional2R];
 ] *)


 \[CapitalLambda]MaxPlot =
   Min[\[CapitalLambda]MaxFunctional, \[CapitalLambda]MaxFunctional2L, \
 \[CapitalLambda]MaxFunctional2R];
 (************** Lambda Max ends here.***************)

If [ makePlots == True,
Print["--> Plotting 5D  RGES"];
  plot5D = Plot[
 	    {
  	 4 \[Pi] (CSU4  +
      SU4Corr[\[Mu]Plot, \[CapitalLambda]MaxFunctional, L5, c0,
       c0Prime, c1, c2] /. CSU4Sol),
  	    4 \[Pi] (CSU2L  +
      SU2LCorr[\[Mu]Plot, \[CapitalLambda]MaxFunctional2L, L5, c0,
       c0Prime, c1, c2] /. CSU2LSol),
  	    4 \[Pi] (CSU2R  +
      SU2RCorr[\[Mu]Plot, \[CapitalLambda]MaxFunctional2R, L5, c0,
       c0Prime, c1, c2] /. CSU2RSol)
  },

 	{\[Mu]Plot, mKK5, \[CapitalLambda]MaxPlot},

 	Frame -> True, FrameStyle -> BlackFrame,
 	FrameLabel -> MaTeX[{"\\mu (GeV)", " \\alpha^{-1}"}, FontSize -> 24],
 	BaseStyle -> texStyle,
 	ImageSize -> Large,
 	Epilog -> {Dashed,
   Thickness[0.001], { InfiniteLine[{{mKK5, 0.0}, {mKK5, 0.5}}]},


   Inset[MaTeX["m_{\\text{KK}_5}",
     Magnification -> 1.4], {mKK5 + 0.7*10^5, 10},
    Scaled[{0.5, 1}]]},
 	PlotRange -> {{mKK5 - 0.7*10^5, \[CapitalLambda]MaxPlot + 0.7*10^5}, Automatic},

 	PlotLegends ->
  MaTeX[{"\\frac{4\\pi}{g^2_{4C}}", "\\frac{4\\pi}{g^2_{2L}}",
    "\\frac{4\\pi}{g^2_{2R}}"}, FontSize -> 18],
 	PerformanceGoal -> "Speed",
 	PlotStyle -> {ColorData[97, "ColorList"][[15]],
   ColorData[97, "ColorList"][[2]], ColorData[97, "ColorList"][[14]]}
 ];

 Export["WeinbergPlots/GPS_5D.pdf", plot5D];
]

\[Alpha]2LinvatLambda =  N[4 \[Pi] (CSU2L  +    SU2LCorr[\[CapitalLambda]MaxPlot, \
            \[CapitalLambda]MaxFunctional2L, L5, c0, c0Prime, c1, c2] /. CSU2LSol)];
\[Alpha]2RinvatLambda =  N[4 \[Pi] (CSU2R  + SU2RCorr[\[CapitalLambda]MaxPlot, \
              \[CapitalLambda]MaxFunctional2R, L5, c0, c0Prime, c1, c2] /. CSU2RSol)];
\[Alpha]4CinvatLambda =  N[4 \[Pi] (CSU4  +  SU4Corr[\[CapitalLambda]MaxPlot, \[CapitalLambda]MaxFunctional,
               L5, c0, c0Prime, c1, c2] /. CSU4Sol)];

Print["-------------------------------------"];
Print["At \[CapitalLambda]Max (\[Alpha]4C)^-1 : ", \[Alpha]4CinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2L)^-1 : ", \[Alpha]2LinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2R)^-1 : ", \[Alpha]2RinvatLambda ];
Print["-------------------------------------"];


\[Alpha]2LWeinb = 1/ \[Alpha]2LinvatLambda;
\[Alpha]2RWeinb = 1/ \[Alpha]2RinvatLambda;
\[Alpha]4CWeinb = 1/ \[Alpha]4CinvatLambda;
CassDiff = 1/(12 \[Pi]) (3/5*2 + 2/5* 4);

sin2ThWEvolved = ((\[Alpha]2LWeinb \[Alpha]4CWeinb +
    2/3 \[Alpha]2LWeinb \[Alpha]2RWeinb + \[Alpha]2RWeinb \
\[Alpha]4CWeinb -
    5/3 \[Alpha]2LWeinb \[Alpha]2RWeinb \[Alpha]4CWeinb*
     CassDiff)/(\[Alpha]2RWeinb \[Alpha]4CWeinb))^(-1);
Print["@\[CapitalLambda]Max:",\[CapitalLambda]MaxPlot ,",  have evolved value of sin \[Theta]W :", sin2ThWEvolved];
Print["-------------------------------------"];




\[Alpha]2LinvatMkk5 =  N[4 \[Pi] (CSU2L  +    SU2LCorr[mKK5+1, \
            \[CapitalLambda]MaxFunctional2L, L5, c0, c0Prime, c1, c2] /. CSU2LSol)];
\[Alpha]2RinvatMkk5 =  N[4 \[Pi] (CSU2R  + SU2RCorr[mKK5+1, \
              \[CapitalLambda]MaxFunctional2R, L5, c0, c0Prime, c1, c2] /. CSU2RSol)];
\[Alpha]4CinvatMkk5 =  N[4 \[Pi] (CSU4  +  SU4Corr[mKK5+1, \[CapitalLambda]MaxFunctional,
               L5, c0, c0Prime, c1, c2] /. CSU4Sol)];

(* Print["-------------------------------------"];
Print["At \[CapitalLambda]Max (\[Alpha]4C)^-1 : ", \[Alpha]4CinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2L)^-1 : ", \[Alpha]2LinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2R)^-1 : ", \[Alpha]2RinvatLambda ];
Print["-------------------------------------"]; *)


\[Alpha]2LWeinb2 = 1/ \[Alpha]2LinvatMkk5;
\[Alpha]2RWeinb2 = 1/ \[Alpha]2RinvatMkk5;
\[Alpha]4CWeinb2 = 1/ \[Alpha]4CinvatMkk5;
CassDiff = 1/(12 \[Pi]) (3/5*2 + 2/5* 4);

sin2ThWEvolved2 = ((\[Alpha]2LWeinb2 \[Alpha]4CWeinb2 +
    2/3 \[Alpha]2LWeinb2 \[Alpha]2RWeinb2 + \[Alpha]2RWeinb2 \
\[Alpha]4CWeinb2 -
    5/3 \[Alpha]2LWeinb2 \[Alpha]2RWeinb2 \[Alpha]4CWeinb2 *
     CassDiff)/(\[Alpha]2RWeinb2 \[Alpha]4CWeinb2))^(-1);
Print["@MKK5, have evolved value of sin \[Theta]W :", sin2ThWEvolved2];
Print["-------------------------------------"];



(* EXPORT HERE!!**)
(* Print["At MKK5, (\[Alpha]1Y)^-1 :" ,\[Alpha]1YatMKK5inv];
Print["At MKK5, (\[Alpha]2L)^-1 :" ,\[Alpha]2LatMKK5inv];
Print["At MKK5, (\[Alpha]3C)^-1 :" ,\[Alpha]3CatMKK5inv];
Print["At \[CapitalLambda]Max (\[Alpha]4C)^-1 : ", \[Alpha]4CinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2L)^-1 : ", \[Alpha]2LinvatLambda ];
Print["At \[CapitalLambda]Max (\[Alpha]2R)^-1 : ", \[Alpha]2RinvatLambda ];
Print["@MKK5, have evolved value of sin \[Theta]W :", sin2ThWEvolved2];
 *)
Export[jsonNameOut, <|"sin2ThW-Lambda" -> sin2ThWEvolved, "LambdaMax" -> \[CapitalLambda]MaxPlot,
                      "sin2ThW-MKK5" -> sin2ThWEvolved2, "a1Yinv-MKK5" -> \[Alpha]1YatMKK5inv, "a2Linv-MKK5" -> \[Alpha]2LatMKK5inv, "a3Cinv-MKK5" -> \[Alpha]3CatMKK5inv, "a4Cinv-Lambda" -> \[Alpha]4CinvatLambda, "a2Linv-Lambda" -> \[Alpha]2LinvatLambda, "a2Rinv-Lambda" -> \[Alpha]2RinvatLambda |>];


(* Print[N[mKK5], "  ", \[CapitalLambda]MaxPlot] *)
(* Export[jsonNameOut, <|"sin2ThW" -> sin2ThWEvolved|> ]; *)
