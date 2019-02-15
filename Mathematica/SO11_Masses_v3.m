(****  Constants   **)
(*********** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! *****************)
(*  This version has 1/2 instead of 1 in Q0 for Wpm and Z0 effective potential contributions, as it should be? *)



(* Subscript[sin2\[Theta], W] = 0.2312;
\[Alpha]EM = 1/127.96; *)

timeOut = 20;

(**** Get Json data *****)

JsonNb = $ScriptCommandLine[[2]];
jsonName = "dataIn"  <> JsonNb <> ".json";
jsonNameOut = "massesOut"  <> JsonNb <> ".json";


(* jsonName = "dataIn.json"; *)
dataRule = Import [jsonName];

(* Print[ "GOT JSON" ]; *)


sin2\[Theta]W = 0.2312;
\[Xi]Gauge = 0;

M = -10^7;
mB = 1.145 * 10^12;
\[Alpha]EM = 1/127.96;


\[Mu]11Fct[\[Mu]1_, zL_, c0_,
   c1_] := (4.18 / 172.44) * \[Mu]1 * (1/(zL^(c0 - c1))) *
   Sqrt[(1 + 2 c0)/(1 + 2 c1)];
\[Mu]11PrimeFct[\[Mu]2Tilde_, zL_, c0_, c2_ ] := Which[
   c0 < 0.5 && c2 < 0.5,
   (1.776/172.44) * \[Mu]2Tilde * (1/(zL^(c2 - c0))) * Sqrt[(1 - 2 c0)/(1 - 2 c2)],
   c0 < 0.5 && c2 > 0.5,
   (1.776/172.44) *\[Mu]2Tilde * (1/(zL^(0.5 - c0))) * Sqrt[(1 - 2 c0)/(2 c2 - 1)],
   c0 > 0.5 && c2 < 0.5,
   (1.776/172.44) * \[Mu]2Tilde * (1/(zL^(c2 - 0.5))) * Sqrt[(2 c0 - 1)/(1 - 2 c2)],
   c0 > 0.5 && c2 > 0.5,
   (1.776/172.44) * \[Mu]2Tilde *  Sqrt[(2 c0 - 1)/(2 c2 - 1)]
   ];

(* k = 89130;
zL = 35;

c0 = 0.3325;
c1 = 0.0;
c2 = -0.7;
c0Prime = 0.5224;

\[Mu]11 = 0.108;
\[Mu]1 = 11.18;

\[Mu]11Prime = 0.108;
\[Mu]2Tilde = 0.7091; *)


(********************** Get stuff from JSON *************************)
k = "k" /. dataRule;
zL = "zL" /. dataRule;



c0 = "c0" /.dataRule;
c1 = "c1" /.dataRule;
c2 = "c2" /.dataRule;
c0Prime = "c0Prime" /.dataRule;

\[Mu]1 = "Mu1" /. dataRule;
\[Mu]2Tilde = "Mu2Tilde" /.dataRule  ;


\[Mu]11 = \[Mu]11Fct[\[Mu]1, zL, c0, c1];
\[Mu]11Prime = \[Mu]11PrimeFct[\[Mu]2Tilde, zL, c0, c2];

(********************** Potential declaration *************************)
VeffFCT = -((1/zL^4)*(0.20264236728467555*k^4*q^3*
     Log[1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
           BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
           BesselK[-(1/2) + c0Prime, q/zL])*
         (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
          BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL]))])) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
          BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
          BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL]))]) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
          BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
          BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL] +
         (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
            BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
           (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
            BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
           (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
            BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
          (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
           (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
             BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2)))]) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[1/2 + c0, q/zL]*
          BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
          BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*
          BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
          BesselK[-(1/2) + c0, q/zL] +
         (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
            BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
           (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
            BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
          ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
           \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
             2)))]) - (1/zL^4)*(0.012665147955292222*k^4*q^3*
    Log[1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/
          ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
               BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
            BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
              2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
         1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
               q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
             (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0,
                  q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0,
                  q]))/(k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) -
      (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
         BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
        ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
            BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
         BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
           2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
        Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
              BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
            (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL])])]) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    Log[1 + (zL*Sin[\[Theta]H]^2)/(2*q^2*(BesselI[0, q/zL]*BesselK[0, q] -
         BesselI[0, q]*BesselK[0, q/zL])*
        (BesselI[1, q/zL]*BesselK[1, q] - BesselI[1, q]*
          BesselK[1, q/zL]))]) + (1/zL^4)*(0.018997721932938333*k^4*q^3*
    \[Xi]Gauge^2*Log[1 + (zL*Sin[\[Theta]H]^2)/
       (q^2*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
          BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
         BesselI[1, q]*BesselK[1, q/zL]))]) -
  (1/zL^4)*(0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    Log[1 - (zL*Sin[\[Theta]H]^2)/(2*q^2*(-1 + sin2\[Theta]W)*
        (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
          BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
         BesselI[1, q]*BesselK[1, q/zL]))]) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    Log[1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
         (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
             BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
             BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
             BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
             BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
          Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
         (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
             BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
             BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
             BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
             BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
          (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
            BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
         \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
                q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2,
                q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2,
                q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2)))])

(********************** Higgs Derivative *************************)

HiggsDeriv = -((1/zL^4)*(0.20264236728467555*k^4*q^3*
     (-((zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
           BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
           BesselK[-(1/2) + c0Prime, q/zL])*
         (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
          BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])*
         (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
              BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime,
               q]*BesselK[-(1/2) + c0Prime, q/zL])*
            (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
             BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/
                zL]))))) - (zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[-(1/2) + c0Prime, q/zL]*BesselK[-(1/2) + c0Prime,
            q] - BesselI[-(1/2) + c0Prime, q]*BesselK[-(1/2) + c0Prime,
            q/zL])^2*(BesselI[1/2 + c0Prime, q/zL]*
           BesselK[1/2 + c0Prime, q] - BesselI[1/2 + c0Prime, q]*
           BesselK[1/2 + c0Prime, q/zL])^2*
        (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
              BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime,
               q]*BesselK[-(1/2) + c0Prime, q/zL])*
            (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
             BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])))^
         2) + (zL*Sin[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
          BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime, q]*
          BesselK[-(1/2) + c0Prime, q/zL])*
        (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
         BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime, q/zL])*
        (1 + (zL*Cos[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0Prime, q/zL]*
             BesselK[-(1/2) + c0Prime, q] - BesselI[-(1/2) + c0Prime,
              q]*BesselK[-(1/2) + c0Prime, q/zL])*
           (BesselI[1/2 + c0Prime, q/zL]*BesselK[1/2 + c0Prime, q] -
            BesselI[1/2 + c0Prime, q]*BesselK[1/2 + c0Prime,
              q/zL]))))))) - (1/zL^4)*(0.07599088773175333*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])^2*
        (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
          BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])^2*
        (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
              BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
              BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL])))^2)) +
     (zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0, q/zL]*
         BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
         BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL])*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])))) -
     (zL*Sin[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0, q/zL]*
         BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
         BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL])*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
          (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])))))) -
  (1/zL^4)*(0.07599088773175333*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])^2*
        (BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
          BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL] +
          (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
             BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
            (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
            (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
             BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
           (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
               BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
            (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
              BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))^2*
        (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
              BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
              BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
              BesselK[1/2 + c0, q/zL] + (\[Mu]1^2*(BesselI[1/2 + c0, q]*
                 BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                  q/zL]*BesselK[1/2 + c0, q])*(BesselI[1/2 + c1, q]*
                 BesselK[-(1/2) + c1, q/zL] + BesselI[-(1/2) + c1,
                  q/zL]*BesselK[1/2 + c1, q])*(BesselI[1/2 + c1, q/zL]*
                 BesselK[1/2 + c1, q] - BesselI[1/2 + c1, q]*
                 BesselK[1/2 + c1, q/zL]))/(\[Mu]11^2*
                (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
                  BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
               (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
                 BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))))^
         2)) + (zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0, q/zL]*
         BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
         BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL] +
        (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
          (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
           BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
           BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
            BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
            BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
            BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL] +
           (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
              BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
             (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
              BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
            (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
                BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
               BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))))) -
     (zL*Sin[\[Theta]H/2]^2)/(2*q^2*(BesselI[-(1/2) + c0, q/zL]*
         BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
         BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL] +
        (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
          (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
           BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
           BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
         (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
             BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
          (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
            BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))*
       (1 + (zL*Sin[\[Theta]H/2]^2)/(q^2*(BesselI[-(1/2) + c0, q/zL]*
            BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
            BesselK[-(1/2) + c0, q/zL])*(BesselI[1/2 + c0, q/zL]*
            BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
            BesselK[1/2 + c0, q/zL] +
           (\[Mu]1^2*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
              BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q])*
             (BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
              BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])*
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
              BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL]))/
            (\[Mu]11^2*(BesselI[1/2 + c1, q]*BesselK[-(1/2) + c1, q/zL] +
                BesselI[-(1/2) + c1, q/zL]*BesselK[1/2 + c1, q])^2 +
             (BesselI[1/2 + c1, q/zL]*BesselK[1/2 + c1, q] -
               BesselI[1/2 + c1, q]*BesselK[1/2 + c1, q/zL])^2))))))) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((zL^2*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
       (q^4*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
          BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])^2*
        (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
          (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0,
               q] + BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
             BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
           ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
            \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
              2))^2*(1 + (zL*Sin[\[Theta]H/2]^2)/
           (q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
             BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
            (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
             (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0,
                  q] + BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
               (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*(
                BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
                BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
              ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                 BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^
                2 + \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[
                    -(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                   BesselK[1/2 + c2, q/zL])^2))))^2)) +
     (zL*Cos[\[Theta]H/2]^2)/(2*q^2*(BesselI[1/2 + c0, q/zL]*
         BesselK[1/2 + c0, q] - BesselI[1/2 + c0, q]*
         BesselK[1/2 + c0, q/zL])*(BesselI[-(1/2) + c0, q/zL]*
         BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
         BesselK[-(1/2) + c0, q/zL] +
        (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0, q] +
           BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
           BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
          (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
           BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
         ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
            BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
          \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
               q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
            2))*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0,
                q] + BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
              BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
            ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
             \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                  q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
               2))))) - (zL*Sin[\[Theta]H/2]^2)/
      (2*q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
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
          \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
               q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
            2))*(1 + (zL*Sin[\[Theta]H/2]^2)/
         (q^2*(BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
           BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])*
          (BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
           BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL] +
           (\[Mu]2Tilde^2*(BesselI[1/2 + c0, q/zL]*BesselK[-(1/2) + c0,
                q] + BesselI[-(1/2) + c0, q]*BesselK[1/2 + c0, q/zL])*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
             (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2, q] +
              BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL]))/
            ((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2 +
             \[Mu]11Prime^2*(BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                  q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
               2))))))) - (1/zL^4)*(0.012665147955292222*k^4*q^3*
    (-(((1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0,
                 q/zL]*BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0,
                 q]*BesselK[-(1/2) + c0, q/zL])*
              ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                  BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
               BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0,
                   q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + c0,
                   q/zL]))) + 1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*
                 BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
                 BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*
                 BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                    BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                     q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) -
                2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
           Cos[\[Theta]H/2]*Sin[\[Theta]H/2]) - (8*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]*
           Sin[\[Theta]H/2]^3)/(q^4*(BesselI[-(1/2) + c0, q/zL]*
             BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
             BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
               BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
            BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
              2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
           Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0,
                q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                    c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[
                   1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[
                1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^2/
       (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/
             ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
               BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
              ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                  BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
               BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0,
                   q/zL] + 2*k*(k*q - I*M*zL)*BesselK[1/2 + c0,
                   q/zL]))) + 1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*
                 BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
                 BesselK[-(1/2) + c0, q/zL])*(2*BesselI[1/2 + c0, q/zL]*
                 BesselK[1/2 + c0, q] - (mB^2*zL*(BesselI[1/2 + c0, q]*
                    BesselK[-(1/2) + c0, q/zL] + BesselI[-(1/2) + c0,
                     q/zL]*BesselK[1/2 + c0, q]))/(k*(k*q - I*M*zL)) -
                2*BesselI[1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*
           Sin[\[Theta]H/2]^2) - (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
          (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
            BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
               BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
            BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
              2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
           Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0,
                q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                    c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[
                   1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[
                1/2 + c0, q]*BesselK[1/2 + c0, q/zL])]))^2) +
     ((1/q^2)*(zL*((k*((-k)*q + I*M*zL))/((BesselI[-(1/2) + c0, q/zL]*
              BesselK[-(1/2) + c0, q] - BesselI[-(1/2) + c0, q]*
              BesselK[-(1/2) + c0, q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
             BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
               2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0,
                q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                    c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[
                   1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[
                1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*Cos[\[Theta]H/2]^2) -
       (1/q^2)*(zL*((k*((-k)*q + I*M*zL))/
           ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
             BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
               2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0,
                q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                    c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[
                   1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[
                1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) -
       (12*k*zL^2*(k*q - I*M*zL)*Cos[\[Theta]H/2]^2*Sin[\[Theta]H/2]^2)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
            2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
         Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
              q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])]) +
       (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
            2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
         Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
              q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])]))/
      (1 + (1/q^2)*(2*zL*((k*((-k)*q + I*M*zL))/
           ((BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
             BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
            ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
                BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
             BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
               2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))) +
          1/Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
                q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0,
                q/zL])*(2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0,
                q] - (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) +
                    c0, q/zL] + BesselI[-(1/2) + c0, q/zL]*BesselK[
                   1/2 + c0, q]))/(k*(k*q - I*M*zL)) - 2*BesselI[
                1/2 + c0, q]*BesselK[1/2 + c0, q/zL])])*Sin[\[Theta]H/2]^2) -
       (4*k*zL^2*(k*q - I*M*zL)*Sin[\[Theta]H/2]^4)/
        (q^4*(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0, q] -
          BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
         ((mB^2*zL*BesselI[-(1/2) + c0, q/zL] - 2*k*(k*q - I*M*zL)*
             BesselI[1/2 + c0, q/zL])*BesselK[1/2 + c0, q] +
          BesselI[1/2 + c0, q]*(mB^2*zL*BesselK[-(1/2) + c0, q/zL] +
            2*k*(k*q - I*M*zL)*BesselK[1/2 + c0, q/zL]))*
         Conjugate[(BesselI[-(1/2) + c0, q/zL]*BesselK[-(1/2) + c0,
              q] - BesselI[-(1/2) + c0, q]*BesselK[-(1/2) + c0, q/zL])*
           (2*BesselI[1/2 + c0, q/zL]*BesselK[1/2 + c0, q] -
            (mB^2*zL*(BesselI[1/2 + c0, q]*BesselK[-(1/2) + c0, q/zL] +
               BesselI[-(1/2) + c0, q/zL]*BesselK[1/2 + c0, q]))/
             (k*(k*q - I*M*zL)) - 2*BesselI[1/2 + c0, q]*
             BesselK[1/2 + c0, q/zL])])))) -
  (1/zL^4)*(0.012665147955292222*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((zL^2*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
       (q^4*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
           BesselK[0, q/zL])^2*(BesselI[1, q/zL]*BesselK[1, q] -
          BesselI[1, q]*BesselK[1, q/zL])^2*
        (1 + (zL*Sin[\[Theta]H]^2)/(2*q^2*(BesselI[0, q/zL]*BesselK[0, q] -
             BesselI[0, q]*BesselK[0, q/zL])*(BesselI[1, q/zL]*
              BesselK[1, q] - BesselI[1, q]*BesselK[1, q/zL])))^2)) +
     (zL*Cos[\[Theta]H]^2)/(q^2*(BesselI[0, q/zL]*BesselK[0, q] -
        BesselI[0, q]*BesselK[0, q/zL])*
       (BesselI[1, q/zL]*BesselK[1, q] - BesselI[1, q]*
         BesselK[1, q/zL])*(1 + (zL*Sin[\[Theta]H]^2)/
         (2*q^2*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))) -
     (zL*Sin[\[Theta]H]^2)/(q^2*(BesselI[0, q/zL]*BesselK[0, q] -
        BesselI[0, q]*BesselK[0, q/zL])*
       (BesselI[1, q/zL]*BesselK[1, q] - BesselI[1, q]*
         BesselK[1, q/zL])*(1 + (zL*Sin[\[Theta]H]^2)/
         (2*q^2*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))))) +
  (1/zL^4)*(0.018997721932938333*k^4*q^3*\[Xi]Gauge^2*
    (-((4*zL^2*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/
       (q^4*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
           BesselK[0, q/zL])^2*(BesselI[1, q/zL]*BesselK[1, q] -
          BesselI[1, q]*BesselK[1, q/zL])^2*
        (1 + (zL*Sin[\[Theta]H]^2)/(q^2*(BesselI[0, q/zL]*BesselK[0, q] -
             BesselI[0, q]*BesselK[0, q/zL])*(BesselI[1, q/zL]*
              BesselK[1, q] - BesselI[1, q]*BesselK[1, q/zL])))^2)) +
     (2*zL*Cos[\[Theta]H]^2)/(q^2*(BesselI[0, q/zL]*BesselK[0, q] -
        BesselI[0, q]*BesselK[0, q/zL])*
       (BesselI[1, q/zL]*BesselK[1, q] - BesselI[1, q]*
         BesselK[1, q/zL])*(1 + (zL*Sin[\[Theta]H]^2)/
         (q^2*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))) -
     (2*zL*Sin[\[Theta]H]^2)/(q^2*(BesselI[0, q/zL]*BesselK[0, q] -
        BesselI[0, q]*BesselK[0, q/zL])*
       (BesselI[1, q/zL]*BesselK[1, q] - BesselI[1, q]*
         BesselK[1, q/zL])*(1 + (zL*Sin[\[Theta]H]^2)/
         (q^2*(BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))))) -
  (1/zL^4)*(0.006332573977646111*k^4*q^3*(-3. + \[Xi]Gauge^2)*
    (-((zL^2*Cos[\[Theta]H]^2*Sin[\[Theta]H]^2)/(q^4*(-1 + sin2\[Theta]W)^2*
        (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
           BesselK[0, q/zL])^2*(BesselI[1, q/zL]*BesselK[1, q] -
          BesselI[1, q]*BesselK[1, q/zL])^2*
        (1 - (zL*Sin[\[Theta]H]^2)/(2*q^2*(-1 + sin2\[Theta]W)*
            (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
              BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
             BesselI[1, q]*BesselK[1, q/zL])))^2)) -
     (zL*Cos[\[Theta]H]^2)/(q^2*(-1 + sin2\[Theta]W)*
       (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
         BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
        BesselI[1, q]*BesselK[1, q/zL])*
       (1 - (zL*Sin[\[Theta]H]^2)/(2*q^2*(-1 + sin2\[Theta]W)*
          (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))) +
     (zL*Sin[\[Theta]H]^2)/(q^2*(-1 + sin2\[Theta]W)*
       (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
         BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
        BesselI[1, q]*BesselK[1, q/zL])*
       (1 - (zL*Sin[\[Theta]H]^2)/(2*q^2*(-1 + sin2\[Theta]W)*
          (BesselI[0, q/zL]*BesselK[0, q] - BesselI[0, q]*
            BesselK[0, q/zL])*(BesselI[1, q/zL]*BesselK[1, q] -
           BesselI[1, q]*BesselK[1, q/zL])))))) -
  (1/zL^4)*(0.025330295910584444*k^4*q^3*
    (-((-((2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]*
            Sin[\[Theta]H]^3)/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
             (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                 BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                 BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
              (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
                BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
              (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
                BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
             \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                     c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                     c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[
                    -(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                   BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                   BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                   BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                   BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                   BesselK[1/2 + c2, q/zL])^2)))) +
         (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]*Sin[\[Theta]H]*(2*\[Mu]11Prime^2 +
            (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
             Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
            \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                    c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                    c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[
                   -(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                  BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2))))^2/
       (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
            (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
             Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
            (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
                BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
                BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
                BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
             (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
               BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
             (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
            \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                    c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                    c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[
                   -(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                  BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                  BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                  BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                  BesselK[1/2 + c2, q/zL])^2))))^2) +
     (-((10*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Cos[\[Theta]H]^2*
          Sin[\[Theta]H]^2)/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
           (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
               BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
               BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
               BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*BesselK[
                1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
              BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
              BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
           \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) +
                   c2, q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) +
                   c2, q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[
                  -(1/2) + c2, q/zL] + BesselI[-(1/2) + c2, q/zL]*
                 BesselK[1/2 + c2, q])^2 + (BesselI[1/2 + c2, q/zL]*
                 BesselK[-(1/2) + c2, q] + BesselI[-(1/2) + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2*(BesselI[1/2 + c2, q/zL]*
                 BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
                 BesselK[1/2 + c2, q/zL])^2)))) +
       (2*zL^2*(1 - \[Mu]11Prime^2)*(-1 + \[Mu]11Prime^2)*Sin[\[Theta]H]^4)/
        (q^4*((zL^2*\[Mu]11Prime^4)/q^4 + (1/q^2)*(2*zL*\[Mu]11Prime^4*
            (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])*
            (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])) +
          (1 + \[Mu]11Prime^4)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
          \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2,
                 q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2,
                 q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
              2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) +
       (2*zL^2*(1 - \[Mu]11Prime^2)*Cos[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
          \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2,
                 q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2,
                 q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
              2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))) -
       (2*zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
          \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2,
                 q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2,
                 q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
              2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))/
      (1 + (zL^2*(1 - \[Mu]11Prime^2)*Sin[\[Theta]H]^2*(2*\[Mu]11Prime^2 +
          (1/zL)*(2*q^2*(1 + \[Mu]11Prime^2)*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) - (-1 + \[Mu]11Prime^2)*
           Sin[\[Theta]H]^2))/(q^4*((zL^2*\[Mu]11Prime^4)/q^4 +
          (1/q^2)*(2*zL*\[Mu]11Prime^4*(BesselI[-(1/2) + c2, q/zL]*
              BesselK[-(1/2) + c2, q] - BesselI[-(1/2) + c2, q]*
              BesselK[-(1/2) + c2, q/zL])*(BesselI[1/2 + c2, q/zL]*
              BesselK[1/2 + c2, q] - BesselI[1/2 + c2, q]*
              BesselK[1/2 + c2, q/zL])) + (1 + \[Mu]11Prime^4)*
           (BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2, q] -
             BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2, q/zL])^2*
           (BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
             BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2 +
          \[Mu]11Prime^2*((BesselI[-(1/2) + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] - BesselI[-(1/2) + c2, q]*BesselK[-(1/2) + c2,
                 q/zL])^2*(BesselI[1/2 + c2, q]*BesselK[-(1/2) + c2,
                 q/zL] + BesselI[-(1/2) + c2, q/zL]*BesselK[1/2 + c2,
                 q])^2 + (BesselI[1/2 + c2, q/zL]*BesselK[-(1/2) + c2,
                 q] + BesselI[-(1/2) + c2, q]*BesselK[1/2 + c2, q/zL])^
              2*(BesselI[1/2 + c2, q/zL]*BesselK[1/2 + c2, q] -
               BesselI[1/2 + c2, q]*BesselK[1/2 + c2, q/zL])^2))))))

m\[Nu]Type [mu_, M_, mB_, c0_, zL_] := -((
  2 (mu)^2 *M* ((zL)^(2 c0 + 1)))/((2 c0 + 1) *(mB^2)));
(**********************************************)
Fab[\[Alpha]_, \[Beta]_, u_, v_] :=
  BesselJ[\[Alpha], u] BesselY[\[Beta], v] -
   BesselY[\[Alpha], u] BesselJ[\[Beta], v];
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
    2 \[Lambda] z zL Fab[1, 0, \[Lambda] z, \[Lambda] zL];
 Sbasis[\[Lambda]_, z_, zL_] := -\[Pi]/
    2 \[Lambda] z Fab[1, 1, \[Lambda] z , \[Lambda] zL];
 Cprime[\[Lambda]_, z_, zL_] := \[Pi]/
    2 (\[Lambda]^2) z zL Fab[0, 0, \[Lambda] z, \[Lambda] zL];
 Sprime[\[Lambda]_, z_, zL_] := -\[Pi]/
    2 (\[Lambda]^2) z Fab[0, 1, \[Lambda] z, \[Lambda] zL];
(*********************Basis functions*******************************)
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

(**************************************************************)

(******************  Higgs auxiliary functions  *********************************)
fHfunc[k_, zL_, sin2\[Theta]W_, \[Alpha]EM_] :=
  Sqrt[ (6* sin2\[Theta]W) /(4 * Pi * \[Alpha]EM)] * (k/
     Sqrt[ (1 - 1/zL)*(zL^3 - 1) ]);
find\[Theta]HMinRule[VeffFCT_, q_, \[Theta]H_] :=
  Quiet[FindMinimum[
     Quiet[
      NIntegrate [
       VeffFCT /. {\[Theta]H -> \[Theta]HLocal}, {q,
        0, \[Infinity]}]    ]   , {\[Theta]HLocal, 0, 2 }
     ][[2]]];
(******************  Matter mass  *********************************)

fH = fHfunc[k, zL, sin2\[Theta]W, \[Alpha]EM];

replacementRules = {z -> 1, c0var -> c0,
   zLvar -> zL, \[Mu]1var -> \[Mu]1, \[Mu]11var -> \[Mu]11,
   c1var ->
    c1, \[Mu]2TildeVar -> \[Mu]2Tilde, \[Mu]11PrimeVar -> \
\[Mu]11Prime, c2var -> c2, c0PrimeVar -> c0Prime};

\[Theta]HRule =
  TimeConstrained [find\[Theta]HMinRule[VeffFCT, q, \[Theta]H][[1]] ,
   timeOut];

If[ TrueQ[Head @ \[Theta]HRule == Symbol],
Print["Timeout Reached. Aborting."];
mH = 0.0; \[Theta]Hmin = 0.0; mTop =  0.0; mBottom = 0.0;
mTau = 0.0;
mTauNeutrino = 0.0; mPsiDark =  0.0; mW = 0.0; mZ = 0.0;mZprime = 0.0;
,


\[Theta]Hmin = \[Theta]HLocal /. \[Theta]HRule;

If[\[Theta]Hmin == 0,
 Print["Trivial potential. Aborting"];
 mH = 0.0; \[Theta]Hmin = 0.0; mTop =  0.0; mBottom = 0.0; mTau = 0.0;
 mTauNeutrino = 0.0; mPsiDark =  0.0; mW = 0.0; mZ = 0.0; mZprime = 0.0;
 ,

 mH = Sqrt[
   1/fH^2 *
    Quiet[NIntegrate [
      HiggsDeriv /. {\[Theta]H -> \[Theta]HLocal /. \[Theta]HRule}, \
{q, 0, \[Infinity]}]]];

 uQuarkSol = SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 /. replacementRules;
 bQuarkSol =
  SL0 * SR0  + (Sin[\[Theta]Hmin/2])^2 + ((\[Mu]1var^2) SR0* CR0 *
       SL1 *CR1)/(\[Mu]11var^2 (CR1)^2 - (SL1)^2) /. replacementRules;
 tauLeptonSol =
  SL0 * SR0  + (Sin[\[Theta]Hmin/
        2])^2 + (\[Mu]2TildeVar^2 SL0 CL0 SR2 CL2) / \
(\[Mu]11PrimeVar^2 (CL2)^2  - (SR2)^2) /. replacementRules;
 psiDarkQuarkSol =
  SL0prime * SR0prime + (Cos[\[Theta]Hmin/2])^2 /. replacementRules;

 WbosonSol = (2* Cprime[\[Lambda], z, zLvar]*
      Sbasis[\[Lambda], z,
       zLvar] + \[Lambda] (Sin[\[Theta]Hmin])^2) /.
   replacementRules;


 (*Number of \[Pi]/
 zL units to look for a solution (maybe set to 3 to be sure?)*)
 n = 3;

 \[Lambda]ListTop =
  Quiet[Reduce[
    uQuarkSol == 0 && 0.0 < \[Lambda] <  n \[Pi] /zL, \[Lambda]]];
 \[Lambda]ListWpm =
  Quiet[Reduce[
    WbosonSol == 0 && 0 < \[Lambda] <  n \[Pi] /zL, \[Lambda]]];
 \[Lambda]ListBottom =
  Quiet[FindRoot[bQuarkSol , {\[Lambda], 10^(-30)}]];
 \[Lambda]ListTau =
  Quiet[FindRoot[tauLeptonSol , {\[Lambda], 10^(-30)}]];

 \[Lambda]ListPsiDark =
  Quiet[Reduce[
    psiDarkQuarkSol == 0 &&
     0.0 < \[Lambda] < n \[Pi] /zL, \[Lambda]]];



 mTop =  Extract[Extract[\[Lambda]ListTop , {1}], {2}]*k;
 mBottom =  \[Lambda] * k /. \[Lambda]ListBottom;
 mTau =  \[Lambda] * k /. \[Lambda]ListTau;
 mPsiDark =  Extract[Extract[\[Lambda]ListPsiDark , {1}], {2}]*k;
 mW      =  Extract[Extract[\[Lambda]ListWpm , {1}], {2}]*k;
 mZ   =  mW/ Sqrt[1 - sin2\[Theta]W];
 mZprime =  Extract[Extract[\[Lambda]ListWpm , {2}], {2}]*
   k /Sqrt[1 - sin2\[Theta]W];
 (*****  The 10^9 is here to give us the neutrino in eV *)

 mTauNeutrino = m\[Nu]Type[mTop, M, mB, c0, zL] * 10^9;

 ]
]


(*Uncomment to debug*)

(* Print["Higgs mass of :           ", Abs[mH], " (GeV)"];
Print["Higgs minimum <\[Theta]H> is \
located at :     ", Abs[\[Theta]Hmin], " (GeV)"];
Print ["Top mass of :            ", mTop, " (GeV)"];
Print ["Bottom mass of :         ", mBottom, " GeV"];
Print ["Tau lepton mass of :     ", mTau, " GeV"];
Print ["Tau Neutrino mass of :   ",     mTauNeutrino, " eV"];
Print ["Dark fermion mass of :   ", mPsiDark, " GeV"];
Print["W mass of  :   ", mW,  " GeV"];
Print["Z mass of  :   ",   mZ,  " GeV"];
Print["The 2nd KK mode for the \!\(\*SuperscriptBox[\(Z\), \(0\)]\) \
bosons: ", mZprime]; *)





(******************  Export To JSON  *********************************)
Export[jsonNameOut, {"Higgs" -> Abs[mH],
                      "mTop" -> mTop ,
                      "mBottom" -> mBottom,
                      "mTau"-> mTau,
                      "mNeutrino" -> mTauNeutrino ,
                      "mPsiDark" -> mPsiDark,
                      "ThetaHiggs" -> Abs[\[Theta]Hmin],
                      "mWpm" -> mW,
                      "mZ0" -> mZ,
                      "Mu11" -> \[Mu]11,
                      "Mu11Prime" -> \[Mu]11Prime,
                      "mZprime" -> mZprime
                      } ];
