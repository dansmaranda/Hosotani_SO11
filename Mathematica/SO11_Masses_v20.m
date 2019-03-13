fullList =
 DeleteDuplicates[
  Join[\[Lambda]ListTauGuess, \[Lambda]ListTauFullRange]]
Length[fullList]
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
printMessage[True, "\[CapitalPsi]Dark/ W+-"]
printMessage[False, "\[CapitalPsi]Dark/ W+-"]
