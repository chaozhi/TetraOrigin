(* Mathematica Package *)

(* :Context: TetraOrigin`dosePreprocess`*)
(* :Author: Chaozhi Zheng <chaozhi@gmail.com>*)
(* :Mathematica Version: 9.0.1.0 *)
(* :Description: A package for checking input dossage data.*)

BeginPackage["TetraOrigin`dosePreprocess`"]
(* Exported symbols added here with SymbolName::usage *)  

doseValidation::usage = "doseValidation  "

haploValidation::usage = "haploValidation  "

transformSNPDose::usage = "transformSNPDose  "


Begin["`Private`"] (* Begin Private Context *) 

checkgeneticmap[SNPDose_] :=
    Module[ {ls, ls2,chr,res = True},
        ls = Transpose[SNPDose[[2 ;; 3, 2 ;;]]];
        ls = SplitBy[ls, First];
        chr = ls[[All, 1, 1]];
        If[ Length[chr] != Length[Union[chr]],
            Print["Markers on the same chromosome must be in neighbor columns"];
            res = False;
        ];
        ls2 = ls[[All, All, 2]];
        ls2 = OrderedQ[#] & /@ ls2;
        ls2 = Pick[ls[[All, 1, 1]], ls2, False];
        If[ ls2 =!= {},
            Print["The genetic distances must be strictly increasing on chromosomes ", ls2];
            res=False
        ];
        res
    ]
    
checkDose[SNPDose_,ploidy_] := Module[{ls, dose, pos},
  ls = Union[#] & /@ SNPDose[[4 ;;, 2 ;;]];
  ls = Map[ToString, ls, {2}];
  dose = Join[ToString[#] & /@ Range[0, ploidy], {"NA"}];
  pos = Position[Complement[#, dose] == {} & /@ ls, False];
  If[pos === {}, True,
   Print["The dosages must be ", dose, 
    "; the dosages for individuals ", SNPDose[[4 ;;, 1]][[Flatten[pos]]], 
    " are wrong!"];
   False
   ]
  ]    
  
doseValidation[SNPDose_, ploidy_] := 
 If[(checkgeneticmap[SNPDose] && checkDose[SNPDose, ploidy]),True,
  Abort[]
 ]  
 
haploValidation[haplo_] := Module[{chr, ls, res = True},
  chr = Split[haplo[[2, 2 ;;]]][[All, 1]];
  If[Length[chr] != Length[Union[chr]],
   Print["Markers on the same chromosome must be in neighbor columns"];
   res = False;
   ];
  ls = Complement[ToString[#] & /@ Union[Flatten[haplo[[3 ;;, 2 ;;]]]], {"1", "2"}];
  If[ls =!= {},
   Print["Alleles at markers must be labeled as 1 or 2!"];
   res = False;
   ];
  If[res,True,Abort[]]
  ] 



transformDose[dose_] :=
    Module[ {dosage},
        dosage = Rest[Transpose[Rest[dose]]];
        dosage = SplitBy[dosage, First];
        dosage = Transpose[Rest[Transpose[#]] & /@ dosage];
        dosage
    ]
    
jittersnp[snploc_] :=
    Module[ {factor = 5, gsnp, len, pos, j, delt},
        gsnp = SplitBy[snploc];
        len = Length[#] & /@ gsnp;
        pos = Flatten[Position[len - 1, _?Positive]];
        Do[
         Which[
          1 < j < Length[gsnp],
          delt = Differences[gsnp[[{j - 1, j, j + 1}, -1]]]/(len[[j]] factor);
          gsnp[[j]] +=Join[delt[[1]] Range[-Floor[len[[j]]/2], 0], 
            delt[[2]] Range[1, Ceiling[len[[j]]/2] - 1]],
          j == 1,
          delt = (gsnp[[j + 1, -1]] - gsnp[[j, -1]])/(len[[j]] factor);
          gsnp[[j]] += Range[0, len[[j]] - 1] delt,
          j == Length[gsnp],
          delt = (gsnp[[j, -1]] - gsnp[[j - 1, -1]])/(len[[j]] factor);
          gsnp[[j]] += Range[-len[[j]] + 1, 0] delt
          ], {j, pos}];
        gsnp = Flatten[gsnp];
        If[ Length[Union[gsnp]] < Length[gsnp],
            Print["Wrong in jittering overlapped SNPs!"];
            Abort[],
            gsnp
        ]
    ]
  
jittermap[snpMap_] :=
    Module[ {map = snpMap, ls, pos},
        ls = SplitBy[snpMap[[2 ;;, 2 ;;]], First];
        pos = Flatten[Position[Length[Union[#]] < Length[#] & /@ ls[[All, All, 2]],True]];
        If[ pos =!= {},
            Print["Warning: jittering the overlapped SNPs on chromosomes ", ls[[pos, 1, 1]],"!"];
            ls[[pos, All, 2]] = jittersnp[#] & /@ ls[[pos, All, 2]];
            map[[2 ;;, 2 ;;]] = Flatten[ls, 1];
        ];
        map
    ]    
    
transformSNPDose[SNPDose_List,chrsubset_, snpsubset_] :=
    Module[ {parentID,sibID,founderdose, sibdose,snpseq,chrsubset2, snpsubset2},
        parentID = SNPDose[[4 ;; 5, 1]];
        sibID = SNPDose[[6 ;;, 1]];
        snpseq = Transpose[SNPDose[[;;3]]];
        snpseq = jittermap[snpseq]; (*jitter the genetic map of SNPDose[[3]]*)
        founderdose = Join[SNPDose[[;;2]],SNPDose[[4 ;; 5]]];
        sibdose = Join[SNPDose[[;;2]],SNPDose[[6 ;;]]];
        {founderdose, sibdose} = transformDose[#] & /@ {founderdose, sibdose};
        snpseq = SplitBy[snpseq[[2 ;;]], #[[2]] &];
        chrsubset2 = If[ chrsubset === "All",
                         Range[Length[snpseq]],
                         chrsubset
                     ];
        snpsubset2 = If[ snpsubset === "All",
                         All,
                         snpsubset
                     ];
        {founderdose, sibdose} = #[[All, chrsubset2, snpsubset2]] & /@ {founderdose, sibdose};
        snpseq = snpseq[[chrsubset2, snpsubset2]];
        {parentID,sibID,snpseq, founderdose, sibdose,chrsubset2}
    ]    
          
  
End[] (* End Private Context *)

EndPackage[]