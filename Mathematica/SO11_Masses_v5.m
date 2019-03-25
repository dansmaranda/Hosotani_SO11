(* LAX CONSTRAINTS CONTROLL ACCEPTION / REJECTION!!!!*)

constrList = <|"Higgs" ->       <|"Max" -> 625.0, "Min" -> 1.0|>,
         			 "mTop" ->        <|"Max" -> 870.0, "Min" -> 1.0|>,
         			 "mTau" ->        <|"Max" -> 100.0, "Min" -> 0.0000001|>,
         			 "mBottom" ->     <|"Max" -> 100.0, "Min" -> 0.0000001|>,
         			 "ThetaHiggs"->   <|"Max" -> 1.0, "Min" -> 0.0|>
   |>;




(* Set the output to match MachinePrecision*)
SetSystemOptions[ "MachineRealPrintPrecision" -> Round[$MachinePrecision]];

(*--------------------------------------Fermion masses-------------------------------------------*)
(* Fermion Masses (just the neutrino )*)
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
m\[Nu]Type [mu_, M_, mB_, c0_, zL_] := -((
   2 (mu)^2 *M* ((zL)^(2 c0 + 1)))/((2 c0 + 1) *(mB^2)));

(*---------------------------------------------------------------------------------------------------*)
(*--------------------Fermion and Boson Basis functions S, L, SR, CR, SL, CL -----------------------------------------------*)

(*Bessel basis functions Subscript[FHat,\[Alpha]\[Beta]],Subscript[F, \[Alpha]\[Beta]] *)

Fab[\[Alpha]_, \[Beta]_, u_, v_] :=
  BesselJ[\[Alpha], u] BesselY[\[Beta], v] -   BesselY[\[Alpha], u] BesselJ[\[Beta], v];
FHat[\[Alpha]_, \[Beta]_, u_, v_] :=
  BesselI[\[Alpha], u] * BesselK[\[Beta], v] -   Exp[-I (\[Alpha] - \[Beta] ) * Pi] * BesselK[\[Alpha], u] *
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


(*---------------------------------------------------------------------------------------------------*)
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

(*---------------------------------------------------------------------------------------------------*)
(*---------------Higgs Mass auxiliaries and \[Mu] setting ------------------------------------*)

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
\[Mu]11Fct[\[Mu]1_, zL_, c0_,
   c1_] := (4.18/172.44)*\[Mu]1*(1/(zL^(c0 - c1)))*
   Sqrt[(1 + 2 c0)/(1 + 2 c1)];
\[Mu]11PrimeFct[\[Mu]2Tilde_, zL_, c0_, c2_] := Which[
   c0 < 0.5 && c2 < 0.5,
   (1.776/172.44)*\[Mu]2Tilde*(1/(zL^(c2 - c0)))*
    Sqrt[(1 - 2 c0)/(1 - 2 c2)],
   c0 < 0.5 && c2 > 0.5,
   (1.776/172.44) *\[Mu]2Tilde*(1/(zL^(0.5 - c0)))*
    Sqrt[(1 - 2 c0)/(2 c2 - 1)],
   c0 > 0.5 && c2 < 0.5,
   (1.776/172.44) *\[Mu]2Tilde*(1/(zL^(c2 - 0.5)))*
    Sqrt[(2 c0 - 1)/(1 - 2 c2)],
   c0 > 0.5 && c2 > 0.5,
   (1.776/172.44) *\[Mu]2Tilde*Sqrt[(2 c0 - 1)/(2 c2 - 1)]
   ];


(*---------------------------------------------------------------------------------------------------*)
(*---------------Loading the data / setting values. ------------------------------------*)

Clear[zL, k, \[Mu]11, \[Mu]1, \[Mu]11Prime, \[Mu]2Tilde, c0Prime, c0,
  c1, c2, \[Theta]HRule, \[Theta]Hmin, sin2\[Theta]W, \[Xi]Gauge, M,
  mB, \[Alpha]EM, fH, \[Lambda]ListTop];



timeOut = 25;
timeOutSolve = 15;



(* Uncomment below in the .m file *)

JsonNb=$ScriptCommandLine[[2]];
jsonName="dataIn"<>JsonNb<>".json";
jsonNameOut="massesOut"<>JsonNb<>".json";


dataRule = Import[jsonName];
overWrite = False;
(* Constant setting function.*)
sin2\[Theta]W = 0.2312;
\[Xi]Gauge = 0;
M = -10^7;
mB = 1.145*10^12;
\[Alpha]EM = 1/127.96;

(* k = "k" /. dataRule;
zL = "zL" /. dataRule;
c0 = "c0" /. dataRule;
c1 = "c1" /. dataRule;
c2 = "c2" /. dataRule;
c0Prime = "c0Prime" /. dataRule;

\[Mu]1 = "Mu1" /. dataRule;
\[Mu]2Tilde = "Mu2Tilde" /. dataRule;

\[Mu]11 = "Mu11" /. dataRule;
\[Mu]11Prime = "Mu11Prime" /. dataRule; *)


k = 89130;
zL = 35;
c0 = 0.3325;
c1 = 0.0;
c2 = -0.7;
c0Prime = 0.5224;

\[Mu]1 = 11.18;
\[Mu]2Tilde = 0.7091;

\[Mu]11 = 0.108;
\[Mu]11Prime = 0.108;

fH = fHfunc[k, zL, sin2\[Theta]W, \[Alpha]EM];

(*\[Mu]11=\[Mu]11Fct[\[Mu]1,zL,c0,c1];
\[Mu]11Prime=\[Mu]11PrimeFct[\[Mu]2Tilde,zL,c0,c2];
*)

mKK5 = \[Pi] k /(zL - 1);
replacementRules = {z -> 1, c0var -> c0,
   zLvar -> zL, \[Mu]1var -> \[Mu]1, \[Mu]11var -> \[Mu]11,
   c1var ->
    c1, \[Mu]2TildeVar -> \[Mu]2Tilde, \[Mu]11PrimeVar -> \
\[Mu]11Prime, c2var -> c2, c0PrimeVar -> c0Prime};

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

  Print["\n---------------------------------------------------------------------------------------------------\
"];
  Print[colorANSI, startStopStr, esc["reset"], ": Finding ",
   esc["yellow"], particleStr, esc["reset"], " Solutions"];
  Print["---------------------------------------------------------------------------------------------------\
\n"];
  ]



 printMessage[True, "Higgs Minimum"];
 \[Theta]HRule =
   TimeConstrained [find\[Theta]HMinRule[VeffFCT, q, \[Theta]H][[1]] ,
    timeOut];

 (*Plot[NIntegrate[VeffFCT, {q, 0, \[Infinity]}], {\[Theta]H, -0.25, \
 0.25}, PerformanceGoal\[Rule]"Speed"]*)

 If[ TrueQ[Head @ \[Theta]HRule == Symbol],
 Print[esc["red"],
 "Timeout Reached trying to find the Higgs potential minimum. \
Aborting.", esc["reset"]];
  Export[jsonNameOut, return0Masses[]];
  Quit[];
  ,

  printMessage[False, "Higgs Minimum"];
  \[Theta]Hmin = Mod[\[Theta]HLocal /. \[Theta]HRule, \[Pi]];


  If[\[Theta]Hmin == 0 ||
    Not[overWrite || evalConstr["ThetaHiggs", \[Theta]Hmin]],
   Print[esc["red"],
    "Trivial \[Theta]H. Aborting and Returning 0 Masses.",
    esc["reset"]];
   Export[jsonNameOut, return0Masses[]];
   Quit[];
   ,

   (*Higgs Mass Solution*)
   printMessage[True, "Higgs Mass"];
   mH = Sqrt[
     1/fH^2 *
      Quiet[NIntegrate [
        HiggsDeriv /. {\[Theta]H -> \[Theta]HLocal /. \[Theta]HRule}, \
 {q, 0, \[Infinity]}]]];


   If [Not[overWrite || evalConstr["Higgs", mH]] ,
    Print[esc["red"],
     "Higgs mass doesn't satisfy lax constraints. Aborting",
     esc["reset"]];
    Export[jsonNameOut, return0Masses[]];
    Quit[];
    ,

    printMessage[False, "Higgs Mass"];

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


    uQuarkSol =
     SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 /. replacementRules;
    bQuarkSol =
     SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 + ((\[Mu]1var^2) SR0* CR0 *
          SL1 *CR1)/(\[Mu]11var^2 (CR1)^2 - (SL1)^2) /.
      replacementRules;
    tauLeptonSol =
     SL0 * SR0  + (Sin[\[Theta]Hmin/
           2])^2 + (\[Mu]2TildeVar^2 SL0 CL0 SR2 CL2) / \
 (\[Mu]11PrimeVar^2 (CL2)^2  - (SR2)^2) /. replacementRules;
    psiDarkQuarkSol =
     SL0prime * SR0prime + (Cos[\[Theta]Hmin/2])^2 /.
      replacementRules;
    WbosonSol = (2* Cprime[\[Lambda], z, zLvar]*
         Sbasis[\[Lambda], z,
          zLvar] + \[Lambda] (Sin[\[Theta]Hmin])^2) /. replacementRules;
    WRSol = Cbasis[\[Lambda], z, zLvar] /. replacementRules;
    \[Gamma]Tower = Cprime[\[Lambda], z, zLvar] /. replacementRules;


    (*Tau Solution*)
    printMessage[True, "Tau"];
    tauGuess =
     meTypeLepton[c2, c0, \[Mu]11Prime, \[Mu]2Tilde, k,
       zL, \[Theta]Hmin] / k;
    \[Lambda]ListTauGuess =
     TimeConstrained[
      NSolve[tauLeptonSol == 0 &&
        tauGuess >= \[Lambda] >= 0.0, \[Lambda] \[Element] Reals],
      timeOutSolve];
    \[Lambda]ListTauFullRange  =
     TimeConstrained[
      NSolve[tauLeptonSol == 0 &&
        mKK5/k >= \[Lambda] >= tauGuess, \[Lambda] \[Element] Reals],
      timeOutSolve];
    \[Lambda]ListTau =
     DeleteDuplicates[
      Join[\[Lambda]ListTauGuess, \[Lambda]ListTauFullRange]];


    If[TrueQ[Head @ \[Lambda]ListTauGuess == Symbol] ||
      TrueQ[Head @ \[Lambda]ListTauFullRange == Symbol]  ||
      Length[\[Lambda]ListTau] == 0,
     Print[esc["red"],
      "Timeout reached in finding the \[Tau] solution. Aborting.",
      esc["reset"]];
     Export[jsonNameOut, return0Masses[]];
     Quit[];
     ,

     (*Tau Mass*)
     mTau = k*\[Lambda] /. \[Lambda]ListTau[[1]];

     If[Not[overWrite || evalConstr["mTau", mTau] ],
      Print[esc["red"],
       "\[Tau] mass doesn't obey lax constraints. Aborting.",
       esc["reset"]];
      Export[jsonNameOut, return0Masses[]];
      ,

      printMessage[False, "Tau"];

      (*Bottom Solution*)
      printMessage[True, "Bottom"];
      bottomGuess =
       mdTypeQuark[c0, c1, \[Mu]11, \[Mu]1, k, zL, \[Theta]Hmin]/k;
      \[Lambda]ListBottomGuess =
       TimeConstrained[
        NSolve[bQuarkSol == 0 &&
          bottomGuess >= \[Lambda] >= 0.0, \[Lambda] \[Element] Reals],
         timeOutSolve];
      \[Lambda]ListBottomFullRange  =
       TimeConstrained[
        NSolve[bQuarkSol == 0 &&
          mKK5/k >= \[Lambda] >= bottomGuess, \[Lambda] \[Element]
          Reals], timeOutSolve];
      \[Lambda]ListBottom =
       DeleteDuplicates[
        Join[\[Lambda]ListBottomGuess, \[Lambda]ListBottomFullRange]];


      If[TrueQ[Head @ \[Lambda]ListBottomGuess == Symbol] ||
        TrueQ[Head @ \[Lambda]ListBottomFullRange == Symbol] ||
        Length[\[Lambda]ListBottom] == 0,
       Print[esc["red"],
        "Timeout reached in finding the Bottom Quark solution. \
 Aborting.", esc["reset"]];
       Export[jsonNameOut, return0Masses[]];
       Quit[];
       ,

       (*Bottom Mass*)

       mBottom = \[Lambda] * k /. \[Lambda]ListBottom[[1]];

       If[Not[overWrite || evalConstr["mBottom", mBottom] ],
        Print[esc["red"],
         "Bottom quark mass doesn't obey lax constraints. Aborting.",
         esc["reset"]];
        Export[jsonNameOut, return0Masses[]];
        Quit[];
        ,
        printMessage[False, "Bottom"];



        (*Top Solution*)
        (*Note that the solution range for the  Top KK tower is from 0 to 2* Mkk5/k .
        This is for the future RGE codes.*)

        printMessage[True, "Top"];
        \[Lambda]ListTop =
         TimeConstrained[
          NSolve[uQuarkSol == 0 &&
            2*mKK5/k >= \[Lambda] > 0, \[Lambda] \[Element] Reals],
          timeOutSolve];

        If[
         TrueQ[Head @ \[Lambda]ListTop == Symbol] ||
          Length[\[Lambda]ListTop] == 0 ,
         Print[esc["red"],
          "Timeout reached in finding either the TopQuark solution. \
 Aborting.", esc["reset"]];
         Export[jsonNameOut, return0Masses[]];
         Quit[];
         ,


         mTopKKTower = k*\[Lambda] /. \[Lambda]ListTop;
         mTop = mTopKKTower[[1]];


         If [Not[overWrite || evalConstr["mTop", mTop]  ],
          Print[esc["red"],
           "Top quark mass solution doesn't obey lax constraints. \
 Aborting.", esc["reset"]];
          Export[jsonNameOut, return0Masses[]];
          Quit[];
          ,

          printMessage[False, "Top"];
          printMessage[True, "\[CapitalPsi]Dark & W+-"];

          (*Psi Dark and W boson solutions. Note the 2 mKK5/
          k range*)
          \[Lambda]ListPsiDark =
           TimeConstrained[
            NSolve[psiDarkQuarkSol == 0 &&
              2*mKK5/k >= \[Lambda] > 0, \[Lambda] \[Element] Reals],
            timeOutSolve];
          \[Lambda]ListW =
           TimeConstrained[
            NSolve[WbosonSol == 0 &&
              2*mKK5/k >= \[Lambda] > 0, \[Lambda] \[Element] Reals],
            timeOutSolve];


          If[
           TrueQ[Head @ \[Lambda]ListPsiDark == Symbol] ||
            TrueQ[Head @ \[Lambda]ListW == Symbol] ||
            Length[\[Lambda]ListPsiDark] == 0 ||
            Length[\[Lambda]ListW] == 0 ,

           Print[esc["red"],
            "Timeout reached in finding either the \[CapitalPsi] Dark \
 or the W+- solution. Aborting.", esc["reset"]];
           Export[jsonNameOut, return0Masses[]];
           Quit[];
           ,


           printMessage[False, "\[CapitalPsi]Dark & W+-"];

           mWKKTower = k*\[Lambda] /. \[Lambda]ListW;
           mPsiDark = k*\[Lambda] /. \[Lambda]ListPsiDark[[1]];
           mW = mWKKTower[[1]];
           mZ   =  mW/ Sqrt[1 - sin2\[Theta]W];
           mZprime = mWKKTower[[2]]/Sqrt[1 - sin2\[Theta]W];
           (*****  The 10^9 is here to give the neutrino in eV *)

            mTauNeutrino = m\[Nu]Type[mTop, M, mB, c0, zL] * 10^9;


           Print["---------------------------------------------------------------------------------------------------\
 "];

           Print["Higgs mass of :           ", Abs[mH], " (GeV)"];

           Print["Higgs minimum <\[Theta]H> is located at :     ",
            Abs[\[Theta]Hmin], " (rads)"];

           Print["---------------------------------------------------------------------------------------------------\
 "];

           Print["The 1st KK mode for the top quark: ", mTop,
            " (GeV)"];(*Print["The 2nd KK mode for the top quark: ",
           mTop2];*)Print["The 1st KK mode for the bottom quark: ",
            mBottom, " (GeV)"];
           Print["The 1st KK mode for the tau lepton: ", mTau,
            " (GeV)"];
           Print["The 1st KK mode for the dark fermion multiplet: ",
            mPsiDark, " (GeV)"];
           Print["The 1st KK mode for the W⁺⁻  bosons: ", mW,
            " (GeV)"];
           Print["The 1st KK mode for the Z⁰ bosons: ", mZ, " (GeV)"];
           Print["The 2nd KK mode for the Z⁰ bosons: ", mZprime,
            " (GeV)"];
           Print["The 1st KK mode for \[Nu]_\[Tau] neutrinos: ",
            mTauNeutrino, " (GeV)"];
           Print["---------------------------------------------------------------------------------------------------\
 "];(******************Export To JSON*********************************)
           Export[jsonNameOut, {"Higgs" -> Abs[mH], "mTop" -> mTop,
             "mBottom" -> mBottom, "mTau" -> mTau,
             "mNeutrino" -> mTauNeutrino, "mPsiDark" -> mPsiDark,
             "ThetaHiggs" -> Abs[\[Theta]Hmin], "mWpm" -> mW,
             "mZ0" -> mZ,(*"Mu11"\[Rule]\[Mu]11,
             "Mu11Prime"\[Rule]\[Mu]11Prime,*)"mZprime" -> mZprime,
             "Triviality" -> 0}];

           ]
          ]
         ]
        ]
       ]
      ]
     ]
    ]
   ]
  ]
