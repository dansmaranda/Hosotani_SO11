#!/usr/bin/env wolframscript



Subscript[sin2\[Theta], W] = 0.2312;

(* Get Json data *)
Print[$ScriptCommandLine[[2]]];
JsonNb = $ScriptCommandLine[[2]];
jsonName = "dataIn"  <> JsonNb <> ".json";
dataRule = Import [jsonName];



(********************************************************************)
k = "k" /. dataRule;
zL = "zL" /. dataRule;


\[Xi]Gauge = 0;

Subscript[c, 0] = "c0" /.dataRule;
Subscript[c, 1] = "c1" /.dataRule;
Subscript[c, 2] = "c2" /.dataRule;
c0Prime = "c0Prime" /.dataRule;

\[Mu]11 = "Mu11" /. dataRule;
\[Mu]1 = "Mu1" /. dataRule;

\[Mu]11Prime = "Mu11Prime" /.dataRule;
  \[Mu]2Tilde = "Mu2Tilde" /.dataRule  ;

M = -10^7;
mB = 1.145 * 10^12;

\[Alpha]EM = 1/127.96;

fH[k_, zL_, sin2\[Theta]W_, \[Alpha]EM_] :=
  Sqrt[ (6* sin2\[Theta]W) /(4 * Pi * \[Alpha]EM)] * (k/
     Sqrt[ (1 - 1/zL)*(zL^3 - 1) ]);
fH = fH[k, zL, Subscript[sin2\[Theta], W], \[Alpha]EM];

(********************************************************************)



muTypeQuark[c0_] := If[c0 < 0.5, (k/zL)*Sqrt[1 - 4*c0^2]*Sin[Subscript[\[Theta], H]/2],
   k*zL^(-2^(-1) - Subscript[c, 0])*Sqrt[4*c0^2 - 1]*Sin[Subscript[\[Theta], H]/2]]
mdTypeQuark[c0_, c1_, \[Mu]11_, \[Mu]1_] := If[c0 < 0.5, k*zL^(c0 - c1 - 1)*(\[Mu]11/\[Mu]1)*Sqrt[(1 - 2*c0)*(1 + 2*c1)]*
    Sin[Subscript[\[Theta], H]/2], k*zL^(-c1 - 1/2)*(\[Mu]11/\[Mu]1)*Sqrt[(2*c0 - 1)*(2*c1 + 1)*Sin[Subscript[\[Theta], H]/2]]]
meTypeLepton[c2_, c0_, \[Mu]11Prime_, \[Mu]2Tilde_] := If[c2 < 0.5, k*zL^(-1 + c2 - c0)*(\[Mu]11Prime/\[Mu]2Tilde)*
    Sqrt[(2*c0 + 1)*(1 - 2*c2)]*Sin[Subscript[\[Theta], H]/2], k*zL^(-2^(-1) - c0)*(\[Mu]11Prime/\[Mu]2Tilde)*
    Sqrt[(2*c0 + 1)*(2*c2 - 1)]*Sin[Subscript[\[Theta], H]/2]]
m\[Nu]Type[mu_, M_, mB_, c0_] := -((2*mu^2*M*zL^(2*c0 + 1))/((2*c0 + 1)*mB^2));
mPsiDarkType[c0Prime_] := If[c0Prime < 0.5, (k/zL)*Sqrt[1 - 4*c0Prime^2]*Cos[Subscript[\[Theta], H]/2],
    k*zL^(-2^(-1) - c0Prime)*Sqrt[4*c0Prime^2 - 1]*Cos[Subscript[\[Theta], H]/2]];
(********************************************************************)




VeffFct = -((1/zL^4)*(0.20264236728467555*k^4*q^3*Log[1 + (zL*Cos[Subscript[\[Theta], H]/2]^2)/
        (q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
           BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
          BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL]))])) +
  (0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*Log[(E^((4*q)/zL)*(1 + q)*(q - zL) + 2*E^((2*q*(1 + zL))/zL)*zL +
       E^(4*q)*(-1 + q)*(q + zL) - 2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*Subscript[\[Theta], H]])/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)))])/zL^4 -
  (0.07599088773175333*k^4*q^3*Log[1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
       (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
         BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
        (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
          BesselK[1/2 + Subscript[c, 0], q/zL]))])/zL^4 -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    Log[1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0],
           q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
        (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
          BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*
             BesselK[-(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*
             BesselK[1/2 + Subscript[c, 0], q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1],
              q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
           (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
             BesselK[1/2 + Subscript[c, 1], q/zL]))/(\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[
                -(1/2) + Subscript[c, 1], q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1],
                q])^2 + (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] -
             BesselI[1/2 + Subscript[c, 1], q]*BesselK[1/2 + Subscript[c, 1], q/zL])^2)))]) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
         BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
        (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
         BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
           (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
          ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2)))]) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*
    Log[1 + (1/q^2)*(2*((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*
             BesselK[-(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
              q/zL])*((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0],
                q/zL])*BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
             (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) +
                Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
             (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
             2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Sin[Subscript[\[Theta], H]/2]^2) -
      (4*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Sin[Subscript[\[Theta], H]/2]^4)/
       (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
         BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
        ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
          BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
          (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
        Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
          (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
           (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
              BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
           2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])]) -
  (0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    Log[((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)) +
       2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/((E^(2*q) - E^((2*q)/zL))*
       ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)))])/zL^4 -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2*(2*\[Mu]11Prime^2 +
         (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
             BesselK[1/2 + Subscript[c, 2], q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/
       (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
             BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2],
              q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] -
            BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) +
         (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
          (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
             BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
          ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
            (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] + BesselI[-(1/2) + Subscript[c, 2],
                q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 + (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[
                -(1/2) + Subscript[c, 2], q] + BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^
             2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2],
                q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2)))]) -
  (0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*Log[1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/
       ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
        (-1 + Subscript[sin2\[Theta], W]))])/zL^4 ;
(*******************************************************************************************************************)
HiggsDeriv = 0. - (1/zL^4)*(0.20264236728467555*k^4*q^3*
    (-((zL*Cos[Subscript[\[Theta], H]/2]^2)/(2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
         BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
        (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])*
        (1 + (zL*Cos[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
            BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
             BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL]))))) -
     (zL^2*Cos[Subscript[\[Theta], H]/2]^2*Sin[Subscript[\[Theta], H]/2]^2)/
      (q^4*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
          BesselK[-(1/2) + c0Prime, q/zL])^2*(BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
         BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])^2*
       (1 + (zL*Cos[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
            BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
             BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])))^2) +
     (zL*Sin[Subscript[\[Theta], H]/2]^2)/(2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
        BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*
       (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])*
       (1 + (zL*Cos[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime, q] -
           BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime, q/zL])*(BesselI[1/2 + c0Prime, q/zL]*
            BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*(-((zL^2*Cos[Subscript[\[Theta], H]/2]^2*Sin[Subscript[\[Theta], H]/2]^2)/
       (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
          BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])^2*
        (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
           BesselK[1/2 + Subscript[c, 0], q/zL])^2*(1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
           (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
             BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
              BesselK[1/2 + Subscript[c, 0], q/zL])))^2)) + (zL*Cos[Subscript[\[Theta], H]/2]^2)/
      (2*q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
       (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
         BesselK[1/2 + Subscript[c, 0], q/zL])*(1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
         (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
          (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[1/2 + Subscript[c, 0], q/zL])))) - (zL*Sin[Subscript[\[Theta], H]/2]^2)/
      (2*q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
       (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
         BesselK[1/2 + Subscript[c, 0], q/zL])*(1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
         (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
          (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[1/2 + Subscript[c, 0], q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*(-((zL^2*Cos[Subscript[\[Theta], H]/2]^2*Sin[Subscript[\[Theta], H]/2]^2)/
       (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
          BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])^2*
        (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
           BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) +
                Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q])*
            (BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
             BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
            (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
              BesselK[1/2 + Subscript[c, 1], q/zL]))/(\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*
                BesselK[-(1/2) + Subscript[c, 1], q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*
                BesselK[1/2 + Subscript[c, 1], q])^2 + (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[
                1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*BesselK[1/2 + Subscript[c, 1], q/zL])^2))^2*
        (1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) +
                Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
              BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[
                  -(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0],
                  q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
                BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*(BesselI[1/2 + Subscript[c, 1],
                  q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*BesselK[
                  1/2 + Subscript[c, 1], q/zL]))/(\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) +
                     Subscript[c, 1], q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])^
                 2 + (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] -
                 BesselI[1/2 + Subscript[c, 1], q]*BesselK[1/2 + Subscript[c, 1], q/zL])^2))))^2)) +
     (zL*Cos[Subscript[\[Theta], H]/2]^2)/(2*q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
       (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
         BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[-(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*
            BesselK[1/2 + Subscript[c, 0], q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1],
             q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
          (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
            BesselK[1/2 + Subscript[c, 1], q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
             BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])^2 +
          (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
             BesselK[1/2 + Subscript[c, 1], q/zL])^2))*(1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
         (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
          (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[
                -(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0],
                q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
              BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
             (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1],
                q]*BesselK[1/2 + Subscript[c, 1], q/zL]))/(\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*
                 BesselK[-(1/2) + Subscript[c, 1], q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*
                 BesselK[1/2 + Subscript[c, 1], q])^2 + (BesselI[1/2 + Subscript[c, 1], q/zL]*
                BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*BesselK[1/2 + Subscript[c, 1],
                 q/zL])^2))))) - (zL*Sin[Subscript[\[Theta], H]/2]^2)/
      (2*q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
       (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
         BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[-(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*
            BesselK[1/2 + Subscript[c, 0], q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1],
             q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
          (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
            BesselK[1/2 + Subscript[c, 1], q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
             BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])^2 +
          (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*
             BesselK[1/2 + Subscript[c, 1], q/zL])^2))*(1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/
         (q^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
          (BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
            BesselK[1/2 + Subscript[c, 0], q/zL] + (\[Mu]1^2*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[
                -(1/2) + Subscript[c, 0], q/zL] + BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0],
                q])*(BesselI[1/2 + Subscript[c, 1], q]*BesselK[-(1/2) + Subscript[c, 1], q/zL] +
              BesselI[-(1/2) + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q])*
             (BesselI[1/2 + Subscript[c, 1], q/zL]*BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1],
                q]*BesselK[1/2 + Subscript[c, 1], q/zL]))/(\[Mu]11^2*(BesselI[1/2 + Subscript[c, 1], q]*
                 BesselK[-(1/2) + Subscript[c, 1], q/zL] + BesselI[-(1/2) + Subscript[c, 1], q/zL]*
                 BesselK[1/2 + Subscript[c, 1], q])^2 + (BesselI[1/2 + Subscript[c, 1], q/zL]*
                BesselK[1/2 + Subscript[c, 1], q] - BesselI[1/2 + Subscript[c, 1], q]*BesselK[1/2 + Subscript[c, 1],
                 q/zL])^2))))))) - (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((zL^2*Cos[Subscript[\[Theta], H]/2]^2*Sin[Subscript[\[Theta], H]/2]^2)/
       (q^4*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] - BesselI[1/2 + Subscript[c, 0], q]*
           BesselK[1/2 + Subscript[c, 0], q/zL])^2*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*
           BesselK[-(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
            q/zL] + (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
             BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
            (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
           ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
            \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
               BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))^2*
        (1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
             BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
            (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
             BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
             (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
                BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*(
                BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*(
                BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
              ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
                 BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
                (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                  BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))))^2)) +
     (zL*Cos[Subscript[\[Theta], H]/2]^2)/(2*q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
        BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
       (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
        (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
          (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
           BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
          (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
           BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
         ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
          \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))*
       (1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
           BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
          (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
              BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
             (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
             (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
            ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
             \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))))) -
     (zL*Sin[Subscript[\[Theta], H]/2]^2)/(2*q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
        BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
       (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
        BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
        (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
          (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
           BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
          (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
           BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
         ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
            BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
          \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))*
       (1 + (zL*Sin[Subscript[\[Theta], H]/2]^2)/(q^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
           BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
          (BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
           BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] +
              BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])*
             (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
             (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL]))/
            ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2 +
             \[Mu]11Prime^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2))))))) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*
    (-(((1/q^2)*(2*((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0],
                 q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
              ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
                BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
                (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                   q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                  -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                  q/zL])*(2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
                (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                   BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
                2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Cos[Subscript[\[Theta], H]/2]*
           Sin[Subscript[\[Theta], H]/2]) - (8*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Cos[Subscript[\[Theta], H]/2]*
           Sin[Subscript[\[Theta], H]/2]^3)/(q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
             BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
             (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
           Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
              BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
             (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
              (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])]))^2/
       (1 + (1/q^2)*(2*((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                 -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                 q/zL])*((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[
                   1/2 + Subscript[c, 0], q/zL])*BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
                (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                   q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                  -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                  q/zL])*(2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
                (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                   BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
                2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Sin[Subscript[\[Theta], H]/2]^2) -
         (4*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Sin[Subscript[\[Theta], H]/2]^4)/
          (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
             BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
             (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
           Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
              BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
             (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
              (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])]))^2) +
     ((1/q^2)*(((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
             BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
              BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
              (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                 q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                q/zL])*(2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
              (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Cos[Subscript[\[Theta], H]/2]^2) -
       (1/q^2)*(((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
             BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
              BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
              (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                 q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                q/zL])*(2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
              (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Sin[Subscript[\[Theta], H]/2]^2) -
       (12*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Cos[Subscript[\[Theta], H]/2]^2*Sin[Subscript[\[Theta], H]/2]^2)/
        (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
          BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
           BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
           (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
         Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
           (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
            (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
               BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
            2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])]) +
       (4*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Sin[Subscript[\[Theta], H]/2]^4)/
        (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
          BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
           BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
           (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
         Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
           (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
            (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
               BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
            2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])]))/
      (1 + (1/q^2)*(2*((k*zL*((-k)*q + I*M*zL))/((BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) +
                Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
              BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
              (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0],
                 q/zL]))) + Conjugate[zL]/Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[
                -(1/2) + Subscript[c, 0], q] - BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0],
                q/zL])*(2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
              (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
              2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])*Sin[Subscript[\[Theta], H]/2]^2) -
       (4*k*zL*(k*q - I*M*zL)*Conjugate[zL]*Sin[Subscript[\[Theta], H]/2]^4)/
        (q^4*(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
          BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + Subscript[c, 0], q/zL] - 2*k*(k*q - I*M*zL)*BesselI[1/2 + Subscript[c, 0], q/zL])*
           BesselK[1/2 + Subscript[c, 0], q] + BesselI[1/2 + Subscript[c, 0], q]*
           (mB^2*zL*BesselK[-(1/2) + Subscript[c, 0], q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + Subscript[c, 0], q/zL]))*
         Conjugate[(BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[-(1/2) + Subscript[c, 0], q] -
            BesselI[-(1/2) + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL])*
           (2*BesselI[1/2 + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q] -
            (mB^2*zL*(BesselI[1/2 + Subscript[c, 0], q]*BesselK[-(1/2) + Subscript[c, 0], q/zL] +
               BesselI[-(1/2) + Subscript[c, 0], q/zL]*BesselK[1/2 + Subscript[c, 0], q]))/(k*(k*q - I*M*zL)) -
            2*BesselI[1/2 + Subscript[c, 0], q]*BesselK[1/2 + Subscript[c, 0], q/zL])])))) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((16*E^((4*q*(1 + zL))/zL)*q^4*Cos[Subscript[\[Theta], H]]^2*Sin[Subscript[\[Theta], H]]^2)/
       ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)) +
         2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)^2) + (4*E^((2*q*(1 + zL))/zL)*q^2*Cos[Subscript[\[Theta], H]]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)) +
       2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2) - (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL)) +
       2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2))) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((-((2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[Subscript[\[Theta], H]]*Sin[Subscript[\[Theta], H]]^3)/
           (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                 BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[
                  -(1/2) + Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2],
                  q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) +
             (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
              (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                 BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                   BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) +
                     Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2],
                    q/zL] + BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
               (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                  BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*
                (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2],
                    q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2)))) + (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[Subscript[\[Theta], H]]*
           Sin[Subscript[\[Theta], H]]*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[
                 -(1/2) + Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2],
                 q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) -
            (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
               BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
              (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                BesselK[1/2 + Subscript[c, 2], q/zL])) + (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[
                 -(1/2) + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2],
                 q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2 +
            \[Mu]11Prime^2*((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
                 BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
               (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] +
                 BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
              (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                 BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*
               (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                  BesselK[1/2 + Subscript[c, 2], q/zL])^2))))^2/
       (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2*(2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*
              (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2],
                 q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                BesselK[1/2 + Subscript[c, 2], q/zL])) - (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/
          (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[
                 -(1/2) + Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2],
                 q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) +
            (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
               BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
             (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*((BesselI[-(1/2) + Subscript[c, 2], q/zL]*
                  BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[
                   -(1/2) + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2],
                   q/zL] + BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
              (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                 BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*
               (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                  BesselK[1/2 + Subscript[c, 2], q/zL])^2))))^2) +
     (-((10*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[Subscript[\[Theta], H]]^2*Sin[Subscript[\[Theta], H]]^2)/
         (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[
                -(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2],
                q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] -
              BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) +
           (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
              BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2],
                q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
            ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
              (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] +
                BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
             (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] +
                BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*
              (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
                 BesselK[1/2 + Subscript[c, 2], q/zL])^2)))) + (2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*
         Sin[Subscript[\[Theta], H]]^4)/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*
            (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
             (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] + BesselI[
                 -(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] + BesselI[-(1/2) + Subscript[c, 2],
                 q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*
                BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2],
                 q/zL])^2))) + (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[Subscript[\[Theta], H]]^2*
         (2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
              BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) +
                Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] -
             BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) -
          (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
             (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] + BesselI[
                 -(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] + BesselI[-(1/2) + Subscript[c, 2],
                 q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*
                BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2],
                 q/zL])^2))) - (2*zL^2*(1 - \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2*
         (2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
              BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) +
                Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] -
             BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) -
          (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
             (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] + BesselI[
                 -(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] + BesselI[-(1/2) + Subscript[c, 2],
                 q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*
                BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2],
                 q/zL])^2))))/(1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2*
         (2*\[Mu]11Prime^2 + (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*
              BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) +
                Subscript[c, 2], q/zL])*(BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] -
             BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2], q/zL])) -
          (-1 + \[Mu]11Prime^2)*Sin[Subscript[\[Theta], H]]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])*
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] -
             BesselI[-(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
           (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*
              BesselK[1/2 + Subscript[c, 2], q/zL])^2 + \[Mu]11Prime^2*
           ((BesselI[-(1/2) + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] - BesselI[
                 -(1/2) + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL])^2*
             (BesselI[1/2 + Subscript[c, 2], q]*BesselK[-(1/2) + Subscript[c, 2], q/zL] + BesselI[
                 -(1/2) + Subscript[c, 2], q/zL]*BesselK[1/2 + Subscript[c, 2], q])^2 +
            (BesselI[1/2 + Subscript[c, 2], q/zL]*BesselK[-(1/2) + Subscript[c, 2], q] + BesselI[-(1/2) + Subscript[c, 2],
                 q]*BesselK[1/2 + Subscript[c, 2], q/zL])^2*(BesselI[1/2 + Subscript[c, 2], q/zL]*
                BesselK[1/2 + Subscript[c, 2], q] - BesselI[1/2 + Subscript[c, 2], q]*BesselK[1/2 + Subscript[c, 2],
                 q/zL])^2)))))) + (0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*
    ((8*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*Subscript[\[Theta], H]])/(E^((4*q)/zL)*(1 + q)*(q - zL) + 2*E^((2*q*(1 + zL))/zL)*zL +
       E^(4*q)*(-1 + q)*(q + zL) - 2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*Subscript[\[Theta], H]]) -
     (16*E^((4*q*(1 + zL))/zL)*q^4*Sin[2*Subscript[\[Theta], H]]^2)/
      (E^((4*q)/zL)*(1 + q)*(q - zL) + 2*E^((2*q*(1 + zL))/zL)*zL + E^(4*q)*(-1 + q)*(q + zL) -
        2*E^((2*q*(1 + zL))/zL)*q^2*Cos[2*Subscript[\[Theta], H]])^2))/zL^4 -
  (1/zL^4)*(0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((16*E^((4*q*(1 + zL))/zL)*q^4*Cos[Subscript[\[Theta], H]]^2*Sin[Subscript[\[Theta], H]]^2)/
       ((E^(2*q) - E^((2*q)/zL))^2*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))^2*
        (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/((E^(2*q) - E^((2*q)/zL))*
            ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*(-1 + Subscript[sin2\[Theta], W])))^2*
        (-1 + Subscript[sin2\[Theta], W])^2)) - (4*E^((2*q*(1 + zL))/zL)*q^2*Cos[Subscript[\[Theta], H]]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/((E^(2*q) - E^((2*q)/zL))*
          ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*(-1 + Subscript[sin2\[Theta], W])))*
       (-1 + Subscript[sin2\[Theta], W])) + (4*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/
      ((E^(2*q) - E^((2*q)/zL))*((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*
       (1 - (2*E^((2*q*(1 + zL))/zL)*q^2*Sin[Subscript[\[Theta], H]]^2)/((E^(2*q) - E^((2*q)/zL))*
          ((-E^((2*q)/zL))*(1 + q)*(q - zL) + E^(2*q)*(-1 + q)*(q + zL))*(-1 + Subscript[sin2\[Theta], W])))*
       (-1 + Subscript[sin2\[Theta], W]))));

\[Theta]HMinRule =
  Quiet[FindMinimum[
    Quiet [NIntegrate [
      VeffFct /. {Subscript[\[Theta], H] -> \[Theta]HLocal}, {q,
       0, \[Infinity]}]], {\[Theta]HLocal, 0, 2 }  ][[2]]];

Subscript[m, H] =
  Sqrt[1/fH^2 *
     NIntegrate [
      HiggsDeriv /. {Subscript[\[Theta],
          H] -> \[Theta]HLocal /. \[Theta]HMinRule}, {q,
       0, \[Infinity]}] ];


mHFinal = Abs[Subscript[m, H]];
jsonNameOut = "massesOut"  <> JsonNb <> ".json";




Subscript[m, t] =
muTypeQuark[Subscript[c,
0]] /. {Subscript[\[Theta],
H] -> \[Theta]HLocal /. \[Theta]HMinRule};

mTop = muTypeQuark[Subscript[c,
  0]] /. {Subscript[\[Theta],
    H] -> \[Theta]HLocal /. \[Theta]HMinRule} ;
mBottom = mdTypeQuark[Subscript[c, 0], Subscript[c,
  1], \[Mu]11, \[Mu]1] /. {Subscript[\[Theta],
    H] -> \[Theta]HLocal /. \[Theta]HMinRule} ;

mTau = meTypeLepton[Subscript[c, 2], Subscript[c,
      0], \[Mu]11Prime, \[Mu]2Tilde] /. {Subscript[\[Theta],
        H] -> \[Theta]HLocal /. \[Theta]HMinRule} ;

mNeutrino = m\[Nu]Type[Subscript[m, t], M, mB, Subscript[c, 0]] *
  10^9 /. {Subscript[\[Theta],
    H] -> \[Theta]HLocal /. \[Theta]HMinRule} ;
mPsiDark =  mPsiDarkType[
   c0Prime] /. {Subscript[\[Theta],
     H] -> \[Theta]HLocal /. \[Theta]HMinRule}

Print ["Top mass of : ", mTop , " GeV"]
Print ["Bottom mass of : ", mBottom , " GeV"]
Print ["Tau lepton mass of : ", mTau , " GeV"]
Print ["Tau Neutrino mass of : ", mNeutrino , " eV"]
Print ["Dark fermion mass of : ", mPsiDark  , " GeV"]



Export[jsonNameOut, {"Higgs" -> mHFinal, "mTop" -> mTop , "mBottom" -> mBottom, "mTau"-> mTau, "mNeutrino" -> mNeutrino , "mPsiDark" -> mPsiDark}];
Print["Higgs mass of :     ", mHFinal, " (GeV)"];
