(*Set to True to enable status printing*)
printStatus = True;
(*Set to True to find and export the Trilinear and Top Yukawa*)
addTrilin = True;
(*Enable to ignore lax constraints*)
overWriteLAXConstr = True;

Print["Adding trilinear: ", addTrilin];
Print["Printing Status: ", printStatus];
Print["Constraints Enabled: ", Not[overWriteLAXConstr] ];
(*---------------------------------------------------------------------------------------------------*)
(* Initialise the Machine 0 *)(* Set the output to match MachinePrecision*)

machineZero = 10^(-$MachinePrecision);
SetSystemOptions[ "MachineRealPrintPrecision" -> Round[$MachinePrecision]];


constrList = <|"Higgs" ->                   <|"Max" -> 750.0, "Min" -> 1.0|>,
         			 "mTop" ->                    <|"Max" -> 1050.0, "Min" -> 1.0|>,
         			"mTau" ->                     <|"Max" -> 70.0,  "Min" -> 0.0000001|>,
         			"mBottom" ->                  <|"Max" -> 90.0,  "Min" -> 0.0000001|>,
         			"ThetaHiggs" ->               <|"Max" -> 3.0,    "Min" -> 0.0|>,
         			"mZprime" ->                  <|"Max" -> 10^10,  "Min" -> 2420.0|>,
         			"mWpm" ->                     <|"Max" -> 600.0,  "Min" -> 1.0|>,
              "mTauTow" ->                  <|"Max" -> 10^10,  "Min" -> 560.0|>,
              "mPsiDark" ->                 <|"Max" -> 10^10,  "Min" -> 690.0|>,
              "mBottomTow" ->               <|"Max" -> 10^10,  "Min" -> 690.0|>
   |>;

(* Mass Approximations for all the Bessel Equations*)

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

(* Higgs Mass Auxiliaries*)

fHfunc[k_, zL_, sin2\[Theta]W_, \[Alpha]EM_] :=
  Sqrt[ (6* sin2\[Theta]W) /(4 * Pi * \[Alpha]EM)] * (k/
     Sqrt[ (1 - 1/zL)*(zL^3 - 1) ]);
find\[Theta]HMinRule[VeffFCT_, q_, \[Theta]H_] :=
  Quiet[FindMinimum[
     Quiet[
      NIntegrate [
       VeffFCT /. {\[Theta]H -> \[Theta]HLocal}, {q,
        0, \[Infinity]}]], {\[Theta]HLocal, 0, 1}
     ][[2]]];


 (**   Auxiliaries **)
 return0Masses[] := {"Higgs" -> 0.0, "mTop" -> 0.0, "mBottom" -> 0.0,
    "mTau" -> 0.0, "mNeutrino" -> 0.0, "mPsiDark" -> 0.0,
    "ThetaHiggs" -> 0.0, "mWpm" -> 0.0, "mZ0" -> 0.0,
    "mZprime" -> 0.0, "Triviality" -> 1};
 evalConstr[strParticle_, particleMass_] := If[
    constrList[strParticle]["Max"] > particleMass >
     constrList[strParticle]["Min"],
    Return[True];
    ,
    Return[False];
    ];

 esc = Association["reset" -> "\033[1;0m", "black" -> "\033[1;30m",
    "red" -> "\033[1;31m", "green" -> "\033[1;32m",
    "yellow" -> "\033[1;33m", "blue" -> "\033[1;34m",
    "magenta" -> "\033[1;35m"];

 printMessage[startStopBool_, particleStr_] := Module[{},
   If[startStopBool, colorANSI = esc["blue"];
    startStopStr = "Started";, colorANSI = esc["green"];
    startStopStr = "Finished"];

   Print["\n\
 ---------------------------------------------------------------------------------------------------\
 "];
   Print[colorANSI, startStopStr, esc["reset"], ": Finding ",
    esc["yellow"], particleStr, esc["reset"], " Solutions"];
   Print["---------------------------------------------------------------------------------------------------\
 \n"];
   ]

(*---------------Effective potential Subscript[V, eff](\[Theta]H) and Higgs Derivative ------------------------------------*)

VeffFCT = -((1/zL^4)*(0.20264236728467555*k^4*q^3*
     Log[1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
           BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
           BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
           BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
           BesselK[1/2 + c0Prime, q/zL]))])) +
  (0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*
    Log[(E^((4*q)/zL)*(1 + q)*(q - zL) + 2*E^((2*q*(1 + zL))/zL)*zL +
       E^(4*q)*(-1 + q)*(q + zL) - 2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*\[Theta]H])/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
        E^(2*q)*(-1 + q)*(q + zL)))])/zL^4 -
  (0.07599088773175333*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
          BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
           q/zL])*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL]))])/zL^4 -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
          BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
           q/zL])*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
         (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
            BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
           (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
            BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
           (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
             BesselK[1/2 + c1, q/zL]))/
          (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
           (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
              BesselK[1/2 + c1, q/zL])^2)))]) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
        (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
            BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
            BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
          ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
              BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2)))]) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*
    Log[1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/
          ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[
                1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
             (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*BesselK[
                1/2 + c0, q/zL]))) +
         1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
             (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
              (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
        Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
        ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
            BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
          (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
            BesselK[1/2 + c0, q/zL]))*Conjugate[
         (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
              BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) -
           2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])]) -
  (0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    Log[((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
         E^(2*q)*(-1 + q)*(q + zL)) + 2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
        E^(2*q)*(-1 + q)*(q + zL)))])/zL^4 -
  (0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    Log[1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
        (-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*
          (q + zL)))])/zL^4 - (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
         (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
             BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
              q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
            BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) -
         (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
         (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
              q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
             BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
          (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
             BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
          ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2,
                q/zL]*BesselK[1/2 + c2, q])^2 +
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[
                1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))])
HiggsDeriv = 0. - (1/zL^4)*(0.20264236728467555*k^4*q^3*
    (-((zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
          BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
          BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
          BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
          BesselK[1/2 + c0Prime, q/zL])*(1 + (zL*Cos[\[Theta]H/2]^2)/
          (q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
            BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
           (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
            BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL]))))) -
     (zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
      (q^4*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
         BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])^2*
       (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
         BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])^2*
       (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
             BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
             BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
             BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
             BesselK[1/2 + c0Prime, q/zL])))^2) + (zL*Sin[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
        BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
       (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
        BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])*
       (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
            BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
            BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
            BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
            BesselK[1/2 + c0Prime, q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])^2*
        (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
           BesselK[1/2 + c0, q/zL])^2*
        (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) +
                c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL])))^2)) + (zL*Cos[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
       (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL])*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL])))) - (zL*Sin[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
       (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL])*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])^2*
        (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
           BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
            (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
            (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
              BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*
                BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
                BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1,
                q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))^2*
        (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) +
                c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[
                  -(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0,
                  q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
                BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*(
                BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                 BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*
                   BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
                   BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[
                   1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^
                2))))^2)) + (zL*Cos[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
       (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL] +
        (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
          (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
           BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
            BesselK[1/2 + c1, q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
             BesselK[1/2 + c1, q/zL])^2))*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[
                -(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0,
                q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1,
                q]*BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*
                 BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
                 BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[
                 1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^
              2))))) - (zL*Sin[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
       (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL] +
        (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
          (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
           BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
            BesselK[1/2 + c1, q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
             BesselK[1/2 + c1, q/zL])^2))*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[
                -(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0,
                q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1,
                q]*BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*
                 BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
                 BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[
                 1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^
              2))))))) - (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
           BesselK[1/2 + c0, q/zL])^2*(BesselI[-(1/2) + c0, q/zL]*
           BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
            q/zL] + (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
             BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
             BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
            \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
               BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))^2*
        (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
             BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
            (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
             (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
                BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*(
                BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(
                BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
              ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                 BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*
                (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                  BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^2)) +
     (zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
        BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
       (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
        (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
           BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
           BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
          (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
            BesselK[1/2 + c2, q/zL]))/
         ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
          \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
             BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
              BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL]))/((BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                 -(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                  q/zL])^2))))) - (zL*Sin[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
        BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
       (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
        BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
        (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
           BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
           BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
          (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
            BesselK[1/2 + c2, q/zL]))/
         ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
          \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
             BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
              BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL]))/((BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                 -(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                  q/zL])^2))))))) - (1/zL^4)*(0.012665147955292222*k^4*q^3*
    (-(((1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*
                BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[
                 -(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
                 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
               BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                 2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
            1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
                BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*(
                2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
                (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                   BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
                 (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0,
                  q/zL])])*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]) -
         (8*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
          (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[
                1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
             (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*BesselK[
                1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[
                -(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*
                (k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^
        2/(1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*
                BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[
                 -(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
                 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
               BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                 2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
            1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
                BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*(
                2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
                (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                   BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
                 (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0,
                  q/zL])])*Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
          (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[
                1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
             (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*BesselK[
                1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[
                -(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*
                (k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^
        2) +
     ((1/q^2)*(zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) +
                c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
              (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
                BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*
                (k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
         Cos[\[Theta]H/2]^2) - (1/q^2)*(zL*((k*((-k)*q + I*M*zL))/
           ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
              (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
                BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*
                (k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
         Sin[\[Theta]H/2]^2) - (12*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
           (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0,
              q/zL])]) + (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
           (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))/
      (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) +
                c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
              (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
                BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*
                (k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
         Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*
           (mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0,
              q/zL])])))) - (1/zL^4)*(0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((16*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
       ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
           E^(2*q)*(-1 + q)*(q + zL)) + 2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)^2) +
     (4*E^((2*q*(1 + zL))/zL)*q^2*Cos[\[Theta]H]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
         E^(2*q)*(-1 + q)*(q + zL)) + 2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2) -
     (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
         E^(2*q)*(-1 + q)*(q + zL)) + 2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2))) -
  (1/zL^4)*(0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((16*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
       ((E^(2*q) - E^((2*q)/zL))^2*(-1 + sin2\[Theta]W)^2*
        ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
        (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
            (-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*
              (q + zL))))^2)) - (4*E^((2*q*(1 + zL))/zL)*q^2*Cos[\[Theta]H]^2)/
      ((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*(q - zL) +
        E^(2*q)*(-1 + q)*(q + zL))*(1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/
         ((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*
            (q - zL) + E^(2*q)*(-1 + q)*(q + zL))))) +
     (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
       (-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
          (-1 + sin2\[Theta]W)*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*
            (q + zL))))))) - (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((-((2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/
           (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(
                BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(
                BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
              (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
              (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
              ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                  BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
                (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
                  BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
               (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                  BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*
                (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                   BesselK[1/2 + c2, q/zL])^2)))) +
         (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]*(2*\[Mu]11Prime^2 +
            (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                 -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
             Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                 -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
             ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                 BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
               (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
                 BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
              (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                 BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*
               (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2))))^2/
       (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
            (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                 -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
             Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                 -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
             ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                 BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
               (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
                 BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
              (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                 BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*
               (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2))))^2) +
     (-((10*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2,
                q]*BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
               BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
            ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
              (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
                BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                  q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                 BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                  q/zL])^2)))) + (2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*
         Sin[\[Theta]H]^4)/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2,
                 q/zL]*BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) +
       (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2,
                 q/zL]*BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) -
       (2*zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2,
                 q/zL]*BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))/
      (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[
                 -(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2,
                 q/zL]*BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))))) +
  (1/zL^4)*(0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*
    ((8*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*\[Theta]H])/(E^((4*q)/zL)*(1 + q)*(q - zL) +
       2*E^((2*q*(1 + zL))/zL)*zL + E^(4*q)*(-1 + q)*(q + zL) -
       2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*\[Theta]H]) -
     (16*E^((4*q*(1 + zL))/zL)*q^4*Sin[2*\[Theta]H]^2)/
      (E^((4*q)/zL)*(1 + q)*(q - zL) + 2*E^((2*q*(1 + zL))/zL)*zL +
        E^(4*q)*(-1 + q)*(q + zL) - 2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*\[Theta]H])^2))
HiggsDeriv3 = -((1/zL^4)*(0.20264236728467555*k^4*q^3*
     (-((3*zL^2*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2])/(2*q^4*(BesselI[-(1/2) + c0Prime, q/zL]*
            BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime,
             q/zL])^2*(BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
           BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])^2*
         (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime,
                q] - BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
             (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime,
                q]*BesselK[1/2 + c0Prime, q/zL])))^2)) + (zL*Cos[\[Theta]H/2]*Sin[\[Theta]H/2])/
       (q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
         BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
        (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
          BesselK[1/2 + c0Prime, q/zL])*(1 + (zL*Cos[\[Theta]H/2]^2)/
          (q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
            BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
           (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
             BesselK[1/2 + c0Prime, q/zL])))) - (2*zL^3*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2]^3)/
       (q^6*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
          BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])^3*
        (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
           BesselK[1/2 + c0Prime, q/zL])^3*
        (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime,
               q] - BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
            (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
              BesselK[1/2 + c0Prime, q/zL])))^3) + (3*zL^2*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
       (2*q^4*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
          BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])^2*
        (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
           BesselK[1/2 + c0Prime, q/zL])^2*
        (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime,
               q] - BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
            (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
              BesselK[1/2 + c0Prime, q/zL])))^2)))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*((2*zL^3*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2]^3)/
      (q^6*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^3*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])^3*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])))^3) - (3*zL^2*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2])/
      (2*q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])))^2) + (3*zL^2*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
      (2*q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])))^2) - (zL*Cos[\[Theta]H/2]*Sin[\[Theta]H/2])/
      (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
        BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*((2*zL^3*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2]^3)/
      (q^6*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^3*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
         (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
             BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
            BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*(BesselI[1/2 + c1, q/zL]*
             BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
          (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
               BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
             BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))^3*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*
                BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1,
                   q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                 BesselK[1/2 + c1, q/zL])^2))))^3) - (3*zL^2*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2])/
      (2*q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
         (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
             BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
            BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*(BesselI[1/2 + c1, q/zL]*
             BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
          (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
               BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
             BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*
                BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1,
                   q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                 BesselK[1/2 + c1, q/zL])^2))))^2) + (3*zL^2*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
      (2*q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
         BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
         (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
             BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
            BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*(BesselI[1/2 + c1, q/zL]*
             BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
          (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
               BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
             BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*
                BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1,
                   q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
              (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                 BesselK[1/2 + c1, q/zL])^2))))^2) - (zL*Cos[\[Theta]H/2]*Sin[\[Theta]H/2])/
      (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
        BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
        (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
            BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
           BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*(BesselI[1/2 + c1, q/zL]*
            BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*
              BesselK[1/2 + c1, q])^2 + (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
            BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
              BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*BesselK[
                -(1/2) + c1, q/zL] + BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*BesselK[
                1/2 + c1, q/zL]))/(\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
                BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                BesselK[1/2 + c1, q/zL])^2))))))) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*((2*zL^3*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2]^3)/
      (q^6*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL])^3*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
             BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
             BesselK[1/2 + c2, q/zL]))/((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL])^2))^3*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
            (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
                BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
             ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                 BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^
        3) - (3*zL^2*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2])/
      (2*q^4*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL])^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
             BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
             BesselK[1/2 + c2, q/zL]))/((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL])^2))^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
            (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
                BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
             ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                 BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^
        2) + (3*zL^2*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
      (2*q^4*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL])^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
             BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
             BesselK[1/2 + c2, q/zL]))/((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2,
                q]*BesselK[1/2 + c2, q/zL])^2))^2*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
            (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*
                BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
             ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                 BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^
        2) - (zL*Cos[\[Theta]H/2]*Sin[\[Theta]H/2])/(q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
        BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
       (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL] + (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*
            BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
            BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
           BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
         ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
             BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
              BesselK[1/2 + c2, q/zL])^2))*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] + BesselI[-(1/2) + c0,
                q]*BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[
                -(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
            ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])^2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
               2))))))) - (1/zL^4)*(0.012665147955292222*k^4*q^3*
    ((2*((1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                 q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
              ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
                BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                 2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
            1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
                BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*
                 BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0,
                     q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
                 (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*Cos[\[Theta]H/2]*
           Sin[\[Theta]H/2]) - (8*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
          (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
             BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
              2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
            BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*BesselK[
                1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*
               BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
                 BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^3)/
      (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
         (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
             2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
           BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
              BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                 BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
              (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^3 +
     (-((1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Cos[\[Theta]H/2]*Sin[\[Theta]H/2])) -
       (12*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]^3*Sin[\[Theta]H/2])/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
           BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
            2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
                BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])]) + (20*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
           BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
            2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
                BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])]))/
      (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) +
                c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
              BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*
                k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
              (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                   q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
               BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
           BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
            2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
             BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*
                BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])])) -
     (3*((1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]) -
        (8*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]*Sin[\[Theta]H/2]^3)/
         (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
             2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
           BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
              BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                 BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
              (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))*
       ((1/q^2)*(zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Cos[\[Theta]H/2]^2) -
        (1/q^2)*(zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
              BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) - (12*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]^2*
          Sin[\[Theta]H/2]^2)/(q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
            BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
             2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*Conjugate[
           (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
              BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
             (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                  q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL])]) + (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
         (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
             2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
           BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
              BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                 BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
              (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])))/
      (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
             ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*
               BesselK[1/2 + c0, q] + BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
                2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
           1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
               (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                    q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
                BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
         (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*((mB^2*zL*BesselI[-(1/2) + c0, q/zL] -
             2*k*(k*q - I*M*zL)*BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
           BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] + 2*k*(k*q - I*M*zL)*
              BesselK[1/2 + c0, q/zL]))*Conjugate[(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                 BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
              (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^2)) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    ((128*E^((6*q*(1 + zL))/zL)*q^6*Cos[\[Theta]H]^3*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^3*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^3*
       (1 + (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^3) -
     (48*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]^3*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))^2*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 + (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) +
     (48*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^2*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 + (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) -
     (16*E^((2*q*(1 + zL))/zL)*q^2*Cos[\[Theta]H]*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 + (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
          ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))))) +
  (1/zL^4)*(0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*
    ((1024*E^((6*q*(1 + zL))/zL)*q^6*Cos[\[Theta]H]^3*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^3*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^3*
       (1 + (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^3) -
     (192*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]^3*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))^2*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 + (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) +
     (192*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^2*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 + (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) -
     (32*E^((2*q*(1 + zL))/zL)*q^2*Cos[\[Theta]H]*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 + (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*
          ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))))) -
  (1/zL^4)*(0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((128*E^((6*q*(1 + zL))/zL)*q^6*Cos[\[Theta]H]^3*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^3*
        (-1 + sin2\[Theta]W)^3*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^3*
        (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*
            ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^3)) -
     (48*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]^3*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))^2*
       (-1 + sin2\[Theta]W)^2*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) +
     (48*E^((4*q*(1 + zL))/zL)*q^4*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/((E^(2*q) - E^((2*q)/zL))^2*
       (-1 + sin2\[Theta]W)^2*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*
           ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))^2) +
     (16*E^((2*q*(1 + zL))/zL)*q^2*Cos[\[Theta]H]*Sin[\[Theta]H])/((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[\[Theta]H]^2)/((E^(2*q) - E^((2*q)/zL))*(-1 + sin2\[Theta]W)*
          ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))))))) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    ((2*(-((2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/
           (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(
                BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
              (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*
                   BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
                (BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                   BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                  BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                   BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))) +
         (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]*(2*\[Mu]11Prime^2 +
            (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) -
            (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) +
            (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^
                2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                  BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                 BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^3)/
      (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^3 -
     (3*(-((2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/
          (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
              (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^
                2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                  BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                 BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))) +
        (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]*(2*\[Mu]11Prime^2 +
           (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[
                1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) -
           (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
           (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[
                1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) +
           (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))*
       (-((10*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
          (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
              (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^
                2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                  BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                 BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))) +
        (2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^4)/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) +
        (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]^2*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) -
        (2*zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))))/
      (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[
                -(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))^2 +
     (-((24*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]^3*Sin[\[Theta]H])/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*BesselK[
                -(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[1/2 + c2,
                q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*
                 BesselK[-(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))) +
       (32*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]^3)/
        (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
               BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
            (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) -
       (8*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
        (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
               BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
            (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))/
      (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^2))/
        (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2, q/zL] +
               BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2, q])^2 +
            (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))))
(*---------------------------------------------------------------------------------------------------*)
(*---------------Higgs Mass auxiliaries and \[Mu] setting ------------------------------------*)

Clear[zL, k, \[Mu]11, \[Mu]1, \[Mu]11Prime, \[Mu]2Tilde, c0Prime, c0,
  c1, c2, \[Theta]HRule, \[Theta]Hmin, sin2\[Theta]W, \[Xi]Gauge, M,
  mB, \[Alpha]EM, fH, \[Lambda]ListTop];


timeOut = 60;
timeOutSolve = 60;

(* Uncomment below in the .m file *)
JsonNb=$ScriptCommandLine[[2]];
jsonName="dataIn"<>JsonNb<>".json";
jsonNameOut="massesOut"<>JsonNb<>".json";

(*Comment out in the .m, file*)
(* SetDirectory[NotebookDirectory[]]; *)
(* jsonName = "dataInThreadNb-0.json";
jsonNameOut = "massesOutThreadNb-0.json"; *)


(* Constant setting function.*)
sin2\[Theta]W = 0.2312;
\[Xi]Gauge = 0;
M = -10^7;
mB = 1.145*10^12;
\[Alpha]EM = 1/127.96;


dataRule=Import[jsonName];
k="k"/.dataRule;
zL="zL"/.dataRule;
c0="c0"/.dataRule;
c1="c1"/.dataRule;
c2="c2"/.dataRule;
c0Prime="c0Prime"/.dataRule;

\[Mu]1="Mu1"/.dataRule;
\[Mu]2Tilde="Mu2Tilde"/.dataRule;

\[Mu]11="Mu11"/.dataRule;
\[Mu]11Prime="Mu11Prime"/.dataRule;



(* k = 89130;zL = 35;c0 = 0.3325;c1 = 0.0;c2 = -0.7;c0Prime = 0.5224;\[Mu]1 = 11.18;\[Mu]2Tilde = 0.7091;
\[Mu]11 = 0.108;\[Mu]11Prime = 0.108; *)

(* k= 266559.06;zL =35;c0 =0.3289;c1 =0.0;c2 =-0.7;c0Prime =0.5960;\[Mu]11=0.1730;\[Mu]11Prime=0.1730;
\[Mu]2Tilde =1.138;\[Mu]1 =Sqrt[(1+2 c1)/(1 + 2 c0)]*zL^(c0-c1) * (172.44/4.18)*\[Mu]11; *)

fH = fHfunc[k, zL, sin2\[Theta]W, \[Alpha]EM];

mKK5 = \[Pi] k /(zL - 1);
replacementRules = {z -> 1, c0var -> c0,
   zLvar -> zL, \[Mu]1var -> \[Mu]1, \[Mu]11var -> \[Mu]11,
   c1var ->
    c1, \[Mu]2TildeVar -> \[Mu]2Tilde, \[Mu]11PrimeVar -> \
\[Mu]11Prime, c2var -> c2, c0PrimeVar -> c0Prime};

(*Tower Declarations*)
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


uQuarkSol = SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 /. replacementRules;
bQuarkSol =
  SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 + ((\[Mu]1var^2) SR0* CR0 *
       SL1 *CR1)/(\[Mu]11var^2 (CR1)^2 - (SL1)^2) /.
   replacementRules;
tauLeptonSol =
  SL0 * SR0  + (Sin[\[Theta]Hmin/
        2])^2 + (\[Mu]2TildeVar^2 SL0 CL0 SR2 CL2) / \
(\[Mu]11PrimeVar^2 (CL2)^2  - (SR2)^2) /. replacementRules;
psiDarkQuarkSol =
  SL0prime * SR0prime + (Cos[\[Theta]Hmin/2])^2 /. replacementRules;
WbosonSol = (2* Cprime[\[Lambda], z, zLvar]*
      Sbasis[\[Lambda], z,
       zLvar] + \[Lambda] (Sin[\[Theta]Hmin])^2) /. replacementRules;
WRSol = Cbasis[\[Lambda], z, zLvar] /. replacementRules;
\[Gamma]Tower = Cprime[\[Lambda], z, zLvar] /. replacementRules;


(*Find the Higgs Potential Minimum *)

If[printStatus == True, printMessage[True, "Higgs Minimum"];];
\[Theta]HRule =
  TimeConstrained [find\[Theta]HMinRule[VeffFCT, q, \[Theta]H][[1]] ,
   timeOut];
\[Theta]Hmin = Mod[\[Theta]HLocal /. \[Theta]HRule, \[Pi]];




If[\[Theta]Hmin == 0 ||
  Not[overWriteLAXConstr || evalConstr["ThetaHiggs", \[Theta]Hmin]],
 Print[esc["red"],
  "Trivial \[Theta]H. Aborting and Returning 0 Masses.", esc["reset"]];
 Export[jsonNameOut, return0Masses[]];
 Quit[];
 ,


 (*Higgs Mass Solution*)

 If[printStatus == True, printMessage[True, "Higgs Mass"];];

 mH = Sqrt[ 1/fH^2 *  Quiet[NIntegrate [
      HiggsDeriv /. {\[Theta]H -> \[Theta]HLocal /. \[Theta]HRule}, \
{q, 0, \[Infinity]}]]];




 If [Not[overWriteLAXConstr || evalConstr["Higgs", mH]] ,
  Print[esc["red"],
   "Higgs mass doesn't satisfy lax constraints. Aborting",
   esc["reset"]];
  Export[jsonNameOut, return0Masses[]];
  Quit[];
  ,

  If[printStatus == True, printMessage[False, "Higgs Mass"];];
  ]]

(*-------------------Tower Solvers --------------------------*)

solveLambdaEqFull[\[Lambda]Eqn_, split_: False, splitVal_: 0,
   topBound_: 2*mKK5/k, bottBound_: machineZero] := Block[
   {\[Lambda], \[Lambda]List, \[Lambda]ListSplitLow, \
\[Lambda]ListSplitHigh},

   If [split == False,
    \[Lambda]List =
      TimeConstrained[
       NSolve[\[Lambda]Eqn == 0 &&
         topBound >= \[Lambda] > bottBound, \[Lambda] \[Element]
         Reals], timeOutSolve];
    ,

    \[Lambda]ListSplitLow =
     TimeConstrained[
      NSolve[\[Lambda]Eqn == 0 &&
        splitVal >= \[Lambda] >= machineZero, \[Lambda] \[Element]
        Reals], timeOutSolve];
    \[Lambda]ListSplitHigh =
     TimeConstrained[
      NSolve[\[Lambda]Eqn == 0 &&
        topBound >= \[Lambda] >= splitVal, \[Lambda] \[Element]
        Reals], timeOutSolve];
    (*Print[\[Lambda]ListSplitLow];
    Print[\[Lambda]ListSplitHigh];*)
    \[Lambda]List =
     DeleteDuplicates[
      Join[\[Lambda]ListSplitLow, \[Lambda]ListSplitHigh]];
    ];
   Return[\[Lambda]List];
   ];

solveLambdaEqwFindRoot[\[Lambda]Eqn_, \[Lambda]Guess_, offSet_: 0] :=
  Block[
   {\[Lambda]Sol},
   \[Lambda]Sol =
    Quiet[ TimeConstrained[
      FindRoot[\[Lambda]Eqn, {\[Lambda], \[Lambda]Guess + offSet}]  ,
      timeOutSolve ]];
   Return[\[Lambda]Sol];
   ];
solveLambdaRoutine[\[Lambda]Eqn_, \[Lambda]Guess_, partName_,
  topBound_: 2*mKK5/k, split_: False, offSet_: 0] := Block[
  {\[Lambda]List, \[Lambda]ListatGuess, \[Lambda]ListatNextGuess},
  (*The full routine consists of trying to solve the equation in the \
full Nsolve range.
  If that fails the last resort is using the FindRoot alogrithm to \
find the first two towers and put them in. *)
  (*timeOutSolve=
  machineZero;*)
  \[Lambda]List =
   solveLambdaEqFull [\[Lambda]Eqn, split, \[Lambda]Guess, topBound];


  If[TrueQ[Head @ \[Lambda]List == Symbol] ,
   (*timeOutSolve=
   100;*)
   \[Lambda]ListatGuess =
    solveLambdaEqwFindRoot[\[Lambda]Eqn, \[Lambda]Guess];
   \[Lambda]ListatNextGuess =
    solveLambdaEqwFindRoot[\[Lambda]Eqn, \[Lambda]Guess + \[Pi]/zL];

   \[Lambda]List =
    DeleteDuplicates[
     Join[\[Lambda]ListatGuess, \[Lambda]ListatNextGuess]];

   ];


  If[TrueQ[Head @ \[Lambda]List == Symbol] ||
    Length[\[Lambda]List] == 0 ,
   Print[esc["red"], "Timeout reached in finding ", partName,
    " solution. Aborting.", esc["reset"]];
   Export[jsonNameOut, return0Masses[]];
   Quit[];
   ,

   Return[\[Lambda]List];
   ]
  ]


(*-------------------Running Solvers --------------------------*)
toSolveDict = <|

   "mWpm" -> <|"\[Lambda]Eqn" -> WbosonSol,
     "\[Lambda]Guess" -> mWboson[k, zL, \[Theta]Hmin]/k,
     "TopBound" -> 2.0*mKK5/k, "Split" -> False|>,
   "mTop" -> <|"\[Lambda]Eqn" -> uQuarkSol,
     "\[Lambda]Guess" -> muTypeQuark[c0, k, zL, \[Theta]Hmin ]/k ,
     "TopBound" -> 2.0*mKK5/k, "Split" -> False|>,
   "mPsiDark" -> <|"\[Lambda]Eqn" -> psiDarkQuarkSol,
     "\[Lambda]Guess" ->
      mPsiDarkType[c0Prime, k, zL, \[Theta]Hmin] /k ,
     "TopBound" -> 2.0*mKK5/k, "Split" -> False|>,
   "mTau" -> <|"\[Lambda]Eqn" -> tauLeptonSol,
     "\[Lambda]Guess" ->
      meTypeLepton[c2, c0, \[Mu]11Prime, \[Mu]2Tilde, k,
        zL, \[Theta]Hmin] / k , "TopBound" -> 2.0*mKK5/k,
     "Split" -> False|>,
   "mBottom" -> <|"\[Lambda]Eqn" -> bQuarkSol,
     "\[Lambda]Guess" ->
      mdTypeQuark[c0, c1, \[Mu]11, \[Mu]1, k, zL, \[Theta]Hmin]/k ,
     "TopBound" -> 2.0*mKK5/k, "Split" -> False|>
   |>;

solDict = <||>;
partListNames = Keys[toSolveDict];

For[partIdx = 1, partIdx <= Length[partListNames], partIdx++,
 (*Solver Bits*)
 partName = partListNames[[partIdx]];
 \[Lambda]Eqn = toSolveDict[[partName]][["\[Lambda]Eqn"]];
 \[Lambda]Guess =  toSolveDict[[partName]][["\[Lambda]Guess"]];
 TopBound =  toSolveDict[[partName]][["TopBound"]];
 split = toSolveDict[[partName]][["Split"]];

 (*Solve and Print*)

 If [printStatus == True, printMessage[True, partName];];

 \[Lambda]Sols =
  solveLambdaRoutine[\[Lambda]Eqn, \[Lambda]Guess, partName,
   TopBound, split];
 solDict[[partName]] = \[Lambda]Sols;

 If[printStatus == True, printMessage[False, partName];];

 ]

(*Adding the Z tower which is fixed by the weinberg angle and the W \
mass*)
mZList = {auxVar -> \[Lambda]/Sqrt[1 - sin2\[Theta]W] } /.
    solDict[["mWpm"]] /. {auxVar -> \[Lambda]};
solDict[["mZ0"]] = mZList;


(*-------------------Exporting and Checking --------------------------*)
exportDict = <|"mWpm" -> <|"Soln" -> "mWpm", "Idx" -> 1|>,
   		      "mTop" -> <|"Soln" -> "mTop", "Idx" -> 1|>,
   		      "mZ0" -> <|"Soln" -> "mZ0", "Idx" -> 1|>,
   			"mZprime" -> <|"Soln" -> "mZ0", "Idx" -> 2|>,
   			"mPsiDark" -> <|"Soln" -> "mPsiDark", "Idx" -> 1|>,
   			"mTau" -> <|"Soln" -> "mTau", "Idx" -> 1|>,
   			"mBottom" -> <|"Soln" -> "mBottom", "Idx" -> 1|>,
        (************ Tau and Bottom 2nd soln ***********)
        "mTauTow" -> <|"Soln" -> "mTau", "Idx" -> 2|>,
   			"mBottomTow" -> <|"Soln" -> "mBottom", "Idx" -> 2|>
   |>;
toExp = <||>;



Print["---------------------------------------------------------------------------------------------------\
\n"]
Print[esc["magenta"],"      Cheking and Exporting masses.     ", esc["reset"]]
Print["---------------------------------------------------------------------------------------------------\
\n"]

partListNames = Keys[exportDict];
For[partIdx = 1, partIdx <= Length[partListNames], partIdx++,
 partName = partListNames[[partIdx]];

 partSolnName = exportDict[[partName]][["Soln"]];
 partSolnIdx = exportDict[[partName]][["Idx"]];

 partMass = \[Lambda]*k /. solDict[[partSolnName]][[partSolnIdx]];
 toExp[[partName]] = partMass;

 (* Print["Checking ", partName, " mass of ", partMass] *)
 If[Not[overWriteLAXConstr || evalConstr[partName, partMass] ],
  Print[esc["red"], partName,
   " mass doesn't obey lax constraints. Aborting.", esc["reset"]];
  Export[jsonNameOut, return0Masses[]];
  Quit[];,
  Print[esc["green"], "✔ ", esc["reset"],  partName, " passed constraint ",  constrList[partName], " with a mass of ", partMass, "\n"];
  ]

 ]
 (*Exporting the Higgs, and neutrino masses*)
 If[
  Not[overWriteLAXConstr || evalConstr["Higgs", mH] ],
  Print[esc["red"],
   "Higgs mass doesn't obey lax constraints. Aborting.",
   esc["reset"]];
  Export[jsonNameOut, return0Masses[]]
  Quit[],
  toExp[["Higgs"]] = mH;
  ]

 If[Not[overWriteLAXConstr ||
    evalConstr["ThetaHiggs", \[Theta]Hmin] ],
  Print[esc["red"],
   " Higgs Minimum doesn't obey lax constraints. Aborting.",
   esc["reset"]];
  Export[jsonNameOut, return0Masses[]]
  Quit[],
  toExp[["ThetaHiggs"]] = \[Theta]Hmin;
  ]

(* Neutrino masses *)
toExp[["mNeutrino"]] =
  m\[Nu]Type[toExp[["mTop"]], M, mB, c0, zL] * 10^9;
toExp[["Triviality"]] = 0 ;

(*Adding the trilinear and yukawa couplings*)
(* Factor of 1/6 included in the Trilinear*)
If[addTrilin == True,
  \[Tau]Eff =
   0.16667/fH^3 NIntegrate[
     HiggsDeriv3 /. {\[Theta]H -> \[Theta]HLocal /. \[Theta]HRule}, \
{q, 0, \[Infinity]}];
  yTopSM = Sqrt[2]*172.44/246;
  yTopEff =
   yTopSM*Cos[\[Theta]H] /. {\[Theta]H -> \[Theta]HLocal /. \
\[Theta]HRule};
  toExp[["HiggsTrilin"]] = \[Tau]Eff;
  toExp[["TopYukawa"]] = yTopEff;];


Export[jsonNameOut, toExp];
Print[toExp];
