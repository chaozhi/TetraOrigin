(* Mathematica Package *)
(* Created by the Wolfram Workbench Nov 30, 2014 *)

(* :Context: TetraOrigin`*)
(* :Author: Chaozhi Zheng <chaozhi@gmail.com>*)
(* :Mathematica Version: 9.0.1.0 *)
(* :Description: A package defines user interface for haplotype reconstruciton in teraploids.*)

BeginPackage["TetraOrigin`",{"TetraOrigin`TetraModel`","TetraOrigin`dosePreprocess`","TetraOrigin`TetraAlgorithm`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["TetraOrigin`*"];
ClearAll @@ Names["TetraOrigin`*"];

maxStuck::usage = "maxStuck is an option to specifiy the max number of consecutive iterations that are rejected during each phasing run"

maxIteration::usage = "maxIteration is an option to specify the max number of iterations during each phasing run."

minRepeatRun::usage = "minRepeatRun is an option to specifiy the min number of phasing runs at the same local maximimum, so that phasing algorithm will stop if it statifies the minRepeatRun before reaching the maxPhasingRun."

maxPhasingRun::usage = "maxPhasingRun is an option to specify th max number of phasing runs, so that phasing algorithm will stop if it reaches the maxPhasingRun, even not statifying the minRepeatRun."

bivalentPhasing::usage = "bivalentPhasing is an option for phasing algorithm. If bivalentPhasing = True, phasing algorithm accounts for only bivalent chromosome pairing but not qudrivalent pairing. If bivalentPhasing = False, phasing algorithm accounts for both bivalent chromosome pairing and qudrivalent pairing."

bivalentDecoding::usage = "bivalentPhasing is an option for posterior decoding conditional on given parental haplotypes. If bivalentDecoding = True, posterior decoding accounts for only bivalent chromosome pairing but not qudrivalent pairing. If bivalentDecoding = False, posterior decoding accounts for both bivalent chromosome pairing and qudrivalent pairing."

inferTetraPhase::usage = "inferTetraPhase[SNPDose, chrsubset, snpsubset, eps, epsF,ploidy] returns {founderhaplo, loglhistory} where founderhaplo is the estimated parental haplotypes and loglhistory is the records of log likelihood for each proposed parental haplotypes. The ploidy = 4 for tetraploid species. epsF and eps are the dosage error probabilities for parents and siblings, respectively. SNPDose is the input marker data include the genetic map and the dosages for two parents and their full sibs. chrsubset is a list of indices for linkage groups, e.g. chrsubset={1,3} the first and the third linkage groups will be analyzed, and chrsubset = \"All\" all the likage groups will be analyzed. snpsubset is a list of indices for SNP markers of each linkage group, e.g. snpsubset ={2, 5, 10} the second, the fifth, and the tenth markers of each linkage groups will be analyzed, and snpsubset = \"All\" all the markers will be analyzed. "

relabelHaplo::usage = "relabelHaplo[esthaplo, refhaplo] returns the reordered parental haplotypes of esthaplo, so that the total number of allele mismatches between esthaplo and refhaplo is minimized. esthaplo and refhaplo have the same format as the estimation returned by inferTetraPHase."

loglTetraOrigin::usage = "loglTetraOrigin[SNPDose, chrsubset, snpsubset, eps, founderhaplo, ploidy] returns the log likelihood given the parental haplotype founderhaplo. Refer to inferTetraPhase for estimating parental haplotype and descriptions of other paremters."

inferTetraOrigin::usage = "inferTetraOrigin[SNPDose, chrsubset, snpsubset, eps, founderhaplo, ploidy, outputid] calculates the posterior probabilities for each sib at each SNP marker given the the parental haplotype founderhaplo, and the results are saved in the file \"TetraOrigin_Output_outputid_LinkageGroupA.txt\" for the linkage group A, and so on for the rest linkage groups. Refer to inferTetraPhase for estimating parental haplotype and descriptions of other paremters. \ninferTetraOrigin[SNPDose, chrsubset, snpsubset, eps, epsF, ploidy, outputid] is a combination of {founderhaplo, loglhistory}= inferTetraPhase[SNPDose, chrsubset, snpsubset, eps, epsF,ploidy] and inferTetraOrigin[SNPDose, chrsubset, snpsubset, eps, epsF,ploidy, founderhaplo, outputid]."

saveAsSummaryITO::usage = "saveAsSummaryITO[tetraResultFile, summaryFile] summarizes the results in tetraResultFile produced by inferTetraOrigin for each linkage group, and save the summary in summaryFile. \nsaveAsSummaryITO[tetraResultFile, summaryFile, refhaploFile] performs an addition reordering of the estimated parental haplotypes in tetraResultFile, with respect to the reference haplotypes in refhaploFile; see relabelHaplo[esthaplo, refhaplo]."

getSummaryITO::usage = "getSummaryITO[summaryFile] imports the summary produced by saveAsSummaryITO for furthur analysis in mathematica, and returns {description, snpmap, esthaplo, refhaplo, logllist, siblogl, genotypes, estgenoprob, esthaploprob}, where description is a list of explainations for the rest."

toGridProb::usage = "toGridProb[prob,snploc,grid] transfers the summarized genoprob or haploprob values at the locations of snploc to the probability values at the specified grid locations."

(*mergediploprob::usage = "mergediploprob  "

getchrrule::usage = "getchrrule  "*)

Begin["`Private`"]
(* Implementation of the package *)

Options[inferTetraPhase] = {
    maxStuck -> 10,
    maxIteration -> 100, 
    minRepeatRun -> 3,
    maxPhasingRun -> 20,
    bivalentPhasing ->True
}

Options[loglTetraOrigin] = {bivalentDecoding -> False}      

Options[inferTetraOrigin] = Join[Options[inferTetraPhase], Options[loglTetraOrigin]]

(*dataprobset is a very very large matrix, it is designed as a global constant, to avoid as function paremeter*)        
inferTetraOrigin[inputSNPDose_?(ListQ[#] ||StringQ[#]&), chrsubset_, snpsubset_, inputeps_?NonNegative, epsF_?NonNegative, ploidy_Integer, outputid_String, opts : OptionsPattern[]] :=
    Module[ {SNPDose = inputSNPDose,eps = inputeps,loglhistory,founderhaplo},
        If[ eps==0,
            eps=10^(-10.)
        ];
        If[ StringQ[inputSNPDose],
            If[ !FileExistsQ[inputSNPDose],
                Print["File ", inputSNPDose," does not exist!"];
                Return[$Failed]
            ];
            SNPDose = Import[inputSNPDose,"CSV"];
        ];
        doseValidation[SNPDose, ploidy];
        {founderhaplo,loglhistory} = inferTetraPhase[SNPDose, chrsubset, snpsubset, eps, epsF, ploidy, FilterRules[{opts}, Options[inferTetraPhase]]];
        inferTetraOrigin[SNPDose, chrsubset, snpsubset, eps, ploidy, founderhaplo, outputid,opts];
    ]

inferTetraOrigin[inputSNPDose_?(ListQ[#] ||StringQ[#]&), chrsubset_, snpsubset_, inputeps_?NonNegative, inputfounderhaplo_?(ListQ[#] ||StringQ[#]&), ploidy_Integer,outputid_String, opts : OptionsPattern[]] :=
    Module[ {SNPDose = inputSNPDose,founderhaplo = inputfounderhaplo, eps = inputeps,onlybivalent,parentID,sibID,chrnames,snpseq, founderdose, sibdose,chrsubset2,outputfile, startprob, tranprob, condstates,snpmap,ii,starttime},
        (*data validation*)
        If[ eps==0,
            eps=10^(-10.)
        ];
        If[ StringQ[inputSNPDose],
            If[ !FileExistsQ[inputSNPDose],
                Print["File ", inputSNPDose," does not exist!"];
                Return[$Failed]
            ];
            SNPDose = Import[inputSNPDose,"CSV"];
        ];
        If[ StringQ[inputfounderhaplo],
            If[ !FileExistsQ[inputfounderhaplo],
                Print["File ", inputfounderhaplo," does not exist!"];
                Return[$Failed]
            ];
            founderhaplo = Import[inputfounderhaplo,"CSV"];
        ];
        doseValidation[SNPDose, ploidy];
        haploValidation[founderhaplo];
        {parentID,sibID,snpseq, founderdose, sibdose,chrsubset2} = transformSNPDose[SNPDose,chrsubset, snpsubset];
        If[ founderhaplo[[1, 2 ;;]] != Flatten[snpseq[[All, All, 1]]],
            Print[Style["Warning: the marker IDs of founderhaplo are not consisitent with those in SNPDose!",Red]]
        ];
        founderhaplo = SplitBy[Transpose[founderhaplo[[2 ;;, 2 ;;]]], First][[All, All, 2 ;;]];
        If[ Length[founderhaplo]!=Length[snpseq],
            Print["The number of linkage group in founderhaplo has to ", Length[snpseq]];
            Abort[],
            Do[If[ Dimensions[founderhaplo[[ii]]]!={Length[snpseq[[ii]]],2 ploidy},
                   Print["The dimensions of founderhaplo are wrong!"];
                   Abort[]
               ],{ii,Length[founderhaplo]}]
        ];
        (*inference*)
        onlybivalent = OptionValue[bivalentDecoding];
        condstates = zygoteCondStates[ploidy, onlybivalent];
        chrnames = ToString[#]&/@snpseq[[All, 1, 2]];
        Do[
         starttime = SessionTime[];
         outputfile = "TetraOrigin_Output_" <> outputid <> "_LinkageGroup" <> chrnames[[ii]] <> ".txt";
         Print["Start Date =", DateString[], ". Outputfile=", outputfile];
         PrintTemporary["Pre-computing data likelihood..."];
         {startprob, tranprob} = tetraPriorProcess[snpseq, ii, ploidy, onlybivalent];
         (*dataprobset is calculated as a global variable*)
         tetraDoseLikelihood2[founderhaplo, sibdose, ii, eps, ploidy,onlybivalent];
         (*SNPDose[[;; 3, 1]] the column names*)
         snpmap = Join[{SNPDose[[;; 3, 1]]},snpseq[[ii]]];
         Put[{parentID,sibID,snpmap,Transpose[founderhaplo[[ii]]],eps,ploidy,onlybivalent}, outputfile];
         PrintTemporary["Posterior decoding independently for each offspring conditional on the MAP estimations..."];
         tetraPosteriorDecoding[startprob, tranprob, condstates,outputfile];
         Print["Finish Date =", DateString[], ". Time used in Posteriordecoding = ", Round[SessionTime[] - starttime, 0.1], " Seconds."];
         Print["--------------------------------------------------------------------------------------"], {ii, Length[snpseq]}];
        ClearAll[dataprobset];
    ]

inferTetraPhase[inputSNPDose_?(ListQ[#] ||StringQ[#]&), chrsubset_, snpsubset_, inputeps_?NonNegative, epsF_?NonNegative, ploidy_Integer,opts : OptionsPattern[]] :=
    Module[ {SNPDose = inputSNPDose,eps = inputeps,onlybivalent,minrun,maxrun,maxstuck,maxiteration,parentID,sibID,snpseq, founderdose, sibdose,chrsubset2,startprob, tranprob, 
      fhaploset,fhaploweight, fhaplo, founderhaplo,ii,starttime,phaselogl,phasingbook,snpID,chrID,haploID,loglhistory},        
        (*data validation*)
        If[ eps==0,
            eps=10^(-10.)
        ];
        If[ StringQ[inputSNPDose],
            If[ !FileExistsQ[inputSNPDose],
                Print["File ", inputSNPDose," does not exist!"];
                Return[$Failed]
            ];
            SNPDose = Import[inputSNPDose,"CSV"];
        ];
        doseValidation[SNPDose, ploidy];
        (*inference*)
        {parentID,sibID,snpseq, founderdose, sibdose,chrsubset2} = transformSNPDose[SNPDose,chrsubset, snpsubset];
        {onlybivalent,minrun,maxrun,maxstuck,maxiteration} = OptionValue@{bivalentPhasing,minRepeatRun,maxPhasingRun,maxStuck,maxIteration};
        {founderhaplo,loglhistory} = Transpose[Table[
         starttime = SessionTime[];
         Print["Start Date =", DateString[],". Chromosome  = ",snpseq[[ii,1,2]]];
         PrintTemporary["Pre-computing data likelihood..."];
         {startprob, tranprob} = tetraPriorProcess[snpseq, ii, ploidy, onlybivalent];
         (*dataprobset is calculated as a global variable*)
         {fhaploset,fhaploweight} =  tetraDoseLikelihood[founderdose, sibdose, ii, eps, epsF,ploidy,onlybivalent];
         PrintTemporary["Estimating parent phases by maximizing a posteriori..."];
         phasingbook = tetraPhasing[fhaploweight,startprob, tranprob, ploidy, onlybivalent,minrun,maxrun,maxstuck,maxiteration];
         phaselogl = phasingbook[[All,3]];
         loglhistory = phasingbook[[All,4]];
         fhaplo = First[phasingbook[[First[Ordering[phaselogl,-1]]]]];
         founderhaplo = MapThread[#1[[#2]] &, {fhaploset, fhaplo}];
         Print["Finish Phasing Date =", DateString[],". Time used in Phasing = ", Round[SessionTime[] - starttime, 0.1], " Seconds. log posterior of phasing runs = ",phaselogl];
         {founderhaplo,loglhistory}, {ii, Length[snpseq]}]];
        snpID = Flatten[snpseq[[All, All, 1]]];
        chrID = Flatten[snpseq[[All, All, 2]]];
        haploID = Flatten[Outer[StringJoin, parentID, "_" <> ToString[#] & /@ Range[ploidy]]];
        founderhaplo = Join[Transpose[{Join[{"SNP", "Chromosome"}, haploID]}], Join[{snpID, chrID}, Transpose[Join[Sequence @@ founderhaplo]]], 2];
        ClearAll[dataprobset];
        {founderhaplo,loglhistory}
    ]

loglTetraOrigin[inputSNPDose_?(ListQ[#] ||StringQ[#]&), chrsubset_, snpsubset_, inputeps_?NonNegative, inputfounderhaplo_?(ListQ[#] ||StringQ[#]&),ploidy_Integer,opts : OptionsPattern[]] :=
    Module[ {SNPDose = inputSNPDose,founderhaplo = inputfounderhaplo,eps = inputeps,onlybivalent,parentID,sibID,snpseq, founderdose, sibdose,chrsubset2,startprob, tranprob, condstates,ii,logllist,res},        
        (*data validation*)
        If[ eps==0,
            eps=10^(-10.)
        ];
        If[ StringQ[inputSNPDose],
            If[ !FileExistsQ[inputSNPDose],
                Print["File ", inputSNPDose," does not exist!"];
                Return[$Failed]
            ];
            SNPDose = Import[inputSNPDose,"CSV"];
        ];
        If[ StringQ[inputfounderhaplo],
            If[ !FileExistsQ[inputfounderhaplo],
                Print["File ", inputfounderhaplo," does not exist!"];
                Return[$Failed]
            ];
            founderhaplo = Import[inputfounderhaplo,"CSV"];
        ];
        doseValidation[SNPDose, ploidy];
        haploValidation[founderhaplo];
        {parentID,sibID,snpseq, founderdose, sibdose,chrsubset2} = transformSNPDose[SNPDose,chrsubset, snpsubset];
        If[ founderhaplo[[1, 2 ;;]] != Flatten[snpseq[[All, All, 1]]],
            Print[Style["Warning: the marker IDs of founderhaplo are not consisitent with those in SNPDose!",Red]]
        ];
        founderhaplo = SplitBy[Transpose[founderhaplo[[2 ;;, 2 ;;]]], First][[All, All, 2 ;;]];
        (*calculation*)
        onlybivalent = OptionValue[bivalentDecoding];
        If[ Length[founderhaplo]!=Length[snpseq],
            Print["The number of linkage group in founderhaplo has to ", Length[snpseq]];
            Abort[],
            Do[If[ Dimensions[founderhaplo[[ii]]]!={Length[snpseq[[ii]]],2 ploidy},
                   Print["The dimensions of founderhaplo[[",ii,"]] has to ", {Length[snpseq[[ii]]],ploidy}];
                   Abort[]
               ],{ii,Length[founderhaplo]}]
        ];
        condstates = zygoteCondStates[ploidy, onlybivalent];
        res = Table[              
            {startprob, tranprob} = tetraPriorProcess[snpseq, ii, ploidy, onlybivalent];
             (*dataprobset is calculated as a global variable*)
            tetraDoseLikelihood2[founderhaplo, sibdose, ii, eps, ploidy,onlybivalent];
            logllist = tetraMargLiklihood[startprob, tranprob,condstates];
            Total[Last[calsiblogl[logllist, ploidy, onlybivalent]]], {ii, Length[snpseq]}];
        ClearAll[dataprobset];
        Total[res]
    ]

mergediploprob[postdecode_, ploidy_, onlybivalent_] :=
    Module[ {sibtype,sibpos, siblogl,sibtypeprob,logllist, vvorigprob, states,condstates,condstates2,valents,nstate, nind, nvv,nsnp, 
      diploprob, vvprob,weight, vv,ind,i,j,type,ref,foc,order},
        logllist = postdecode[[All, All, 1]];
        {sibtypeprob,sibtype,sibpos, siblogl} = calsiblogl[logllist, ploidy, onlybivalent];
        states = zygoteStates[ploidy, onlybivalent];
        condstates = zygoteCondStates[ploidy, onlybivalent];
        valents = condstates[[All, 1, 1]];
        vvorigprob = postdecode[[All, All, 2]];
        nstate = Max[condstates[[All, All, 2]]];
        {nind, nvv, nsnp} = Take[Dimensions[vvorigprob],3];
        diploprob = ConstantArray[0, {nind, nsnp, nstate}];
        vvprob = Exp[# - Max[#]] & /@ logllist;
        (*the diplotype states are ordered according to the unordeded genotypes, 
        so that the data likelihoood is the same for the unordered the valent combinations*)
        condstates2 = condstates;
        Do[
          type = Map[Length, condstates2[[i, 1, 1]], {2}][[All, 1]];
          ref = states[[condstates2[[i, 1, 2]]]];
          If[ type[[1]] == 2,
              ref[[All, ;; 2]] = Sort[#] & /@ ref[[All, ;; 2]]
          ];
          If[ type[[2]] == 2,
              ref[[All, 3 ;;]] = Sort[#] & /@ ref[[All, 3 ;;]]
          ];
          Do[
           foc = states[[condstates2[[i, j, 2]]]];
           If[ type[[1]] == 2,
               foc[[All, ;; 2]] = Sort[#] & /@ foc[[All, ;; 2]]
           ];
           If[ type[[2]] == 2,
               foc[[All, 3 ;;]] = Sort[#] & /@ foc[[All, 3 ;;]]
           ];
           order = Flatten[Position[foc, #, {1}, 1, Heads -> False] & /@ ref];
           condstates2[[i, j, 2]] = condstates2[[i, j, 2]][[order]], {j, 2, Length[condstates2[[i]]]}], {i, Length[condstates2]}];
        Do[
         vv = sibpos[[ind]];
         weight = Normalize[vvprob[[ind, vv]], Total];
         Do[             
              diploprob[[ind, All, condstates2[[vv[[i]], j, 2]]]] += (weight[[i]]/Length[condstates2[[vv[[i]]]]]) vvorigprob[[ind, vv[[i]]]], {i,
            Length[vv]}, {j, Length[condstates2[[vv[[i]]]]]}], {ind, nind}];
        {logllist,sibtypeprob,sibtype,siblogl, valents,diploprob}
    ]

toGenoprob[diploprob_, ploidy_, onlybivalent_] :=
    Module[ {states, geno2diplo, prob},
        states = Sort[#] & /@ zygoteStates[ploidy, onlybivalent];
        geno2diplo = Flatten[Position[states, #]] & /@ zygoteGroupStates[ploidy, onlybivalent];
        prob = Transpose[diploprob, {2, 3, 1}];
        prob = Transpose[Total[prob[[#]]] & /@ geno2diplo, {3, 1, 2}];
        prob
    ]

toHaploprob[genoprob_, ploidy_, onlybivalent_] :=
    Module[ {groupstates, mtx, depth, haploprob, i, j},
        groupstates = zygoteGroupStates[ploidy, onlybivalent];
        mtx = ConstantArray[0, {Length[groupstates], 2 ploidy}];
        Do[mtx[[i, groupstates[[i, j]]]] += 1, {i, Length[groupstates]}, {j,Dimensions[groupstates][[2]]}];
        mtx = Transpose[mtx]/(ploidy);
        depth = Depth[genoprob];
        haploprob = Map[mtx.# &, genoprob, {depth - 2}];
        haploprob
    ]
    
getchrrule[esthaplo_, truehaplo_] :=
    Module[ {ploidy, est, true, set, dis, rule,p,ii,res, isreverse,est00},
        ploidy = Length[esthaplo]/2;
        est = Partition[esthaplo, ploidy];
        true = Partition[truehaplo, ploidy];
        set = Permutations[Range[ploidy]];
        (*return rules: a list of i->j means the haplotype i in the newhaplo is the haplotype j in the esthaplo*)
        res = Table[
            est00 = If[ isreverse,
                        Reverse[est],
                        est
                    ];
            dis = Table[Total[Abs[Flatten[true[[p]] - est00[[p, ii]]]]], {ii, set}];
            rule = First[Pick[set, dis, Min[dis]]] + If[ isreverse,
                                                         ploidy (2 - p),
                                                         ploidy (p - 1)
                                                     ];
            rule = Thread[Range[ploidy] + ploidy (p - 1) -> rule];
            {Min[dis], rule}, {isreverse, {False, True}}, {p, Length[est]}];
        dis = Total[res[[All, All, 1]], {2}];
        res = Flatten[#] & /@ res[[All, All, 2]];
        res[[Position[dis, Min[dis]][[1, 1]]]]
    ]

relabelHaplo[esthaplo_, refhaplo_] :=
    Module[ {est, ref, rules, est2, res},
        haploValidation[esthaplo];
        haploValidation[refhaplo];
        If[ esthaplo[[1, 2 ;;]] != refhaplo[[1, 2 ;;]],
            Print[Style["Warning: the marker IDs of esthaplo are not consisitent with those of refhaplo!", Red]]
        ];
        If[ esthaplo[[2, 2 ;;]] != refhaplo[[2, 2 ;;]],
            Print[Style["Warning: the chromosome IDs of esthaplo are not consisitent with those of refhaplo!", Red]]
        ];
        est = Transpose[#] & /@ SplitBy[Transpose[esthaplo[[2 ;;, 2 ;;]]], First][[All, All, 2 ;;]];
        ref = Transpose[#] & /@ SplitBy[Transpose[refhaplo[[2 ;;, 2 ;;]]], First][[All, All, 2 ;;]];
        rules = MapThread[getchrrule, {est, ref}];
        est2 = MapThread[#1[[#2]] &, {est, rules[[All, All, 2]]}];
        res = esthaplo;
        res[[3 ;;, 2 ;;]] = Transpose[Flatten[Transpose[#] & /@ est2, 1]];
        res
    ]
  
saveAsSummaryITO::wrongRefHaplo :=
    "The matrix dimension of the reference haplotypes of two parents has to  =  `1` !"
    
saveAsSummaryITO::RefHaploID :=
    "The SNP and chromosome IDs of the input refernece haplotypes are inconsistent with those in the estimated haplotypes!"    
    
saveAsSummaryITO[tetraResultFile_String?FileExistsQ, summaryFile_String,referencehaplo_:None] :=
    Module[ {res,postdecode,parentID, sibID, snpmap, founderhaplo, eps, ploidy, onlybivalent,logllist, sibtypeprob,sibtype, siblogl, valents,
        diploprob,genoprob,haploprob,refhaplo2,refhaplo,chrrule,haplodis,pos,states,groupstates,groupstates2,valents2,siblogl2,temp,
        snpID,chrID,haploID,vvID,founderhaplo2,logllist2,genotypes,rowID,genoprob2,haplotypes,haploprob2, summary, key = "inferTetraOrigin-Summary"},
        res = ReadList[tetraResultFile];
        postdecode = res[[2 ;;]];
        {parentID, sibID, snpmap, founderhaplo, eps, ploidy, onlybivalent} = res[[1, ;; 7]];
        {logllist, sibtypeprob,sibtype,siblogl,valents,diploprob} = mergediploprob[postdecode, ploidy, onlybivalent];            
        (*relalel 2 ploidy haplotypes of the two parents*)
        If[ referencehaplo === None,
            refhaplo2 = {{"None"}},
            refhaplo2 = If[ Head[referencehaplo]===String,
                            Import[referencehaplo, "CSV"],
                            referencehaplo
                        ];
            haploValidation[refhaplo2];
            If[ refhaplo2[[1, 2 ;;]] != snpmap[[2 ;;, 1]],
                Print[Style["Warning: the marker IDs of refhaplo are not consisitent with those in snpmap!",Red]]
            ];
            If[ ! ((Dimensions[founderhaplo] + {2, 1} == Dimensions[refhaplo2])),
                Message[saveAsSummaryITO::wrongRefHaplo, Dimensions[founderhaplo]+ {2, 1}];
                Abort[];
                Beep[];
            ];
            If[ refhaplo2[[;; 2, 2 ;;]] =!= Transpose[snpmap[[2 ;;, ;; 2]]],
                Message[saveAsSummaryITO::RefHaploID]
            ];
            refhaplo = refhaplo2[[3 ;;, 2 ;;]];
            chrrule = getchrrule[founderhaplo, refhaplo];
            founderhaplo = founderhaplo[[chrrule[[All, 2]]]];
            haplodis = Total[Abs[Flatten[founderhaplo - refhaplo]]];
            Print["The number of mismatches between estimated parental haplotypes and reference haplotypes: "<>ToString[haplodis]];
            (*change the order of two parents, so that the chrosomes in P1 is always labeled as 1234, and 5678 for P2*)
            If[ chrrule[[1, 2]] > ploidy,
                chrrule[[;; ploidy, 2]] -= ploidy;
                chrrule[[ploidy + 1 ;;, 2]] += ploidy;
            ];
            pos = Flatten[Position[valents, #, {1}, 1, Heads -> False] & /@ Map[Sort, valents /. chrrule, {2, 3}]];
            logllist = logllist[[All, pos]];
            states = zygoteStates[ploidy, onlybivalent];
            pos = Flatten[Position[states, #, {1}, 1, Heads -> False] & /@ (states /.chrrule)];
            diploprob = diploprob[[All, All, pos]];
        ];
        genoprob = toGenoprob[diploprob, ploidy, onlybivalent];
        haploprob = toHaploprob[genoprob, ploidy, onlybivalent];
        {logllist,genoprob,haploprob} = N[Round[#,10^(-5)]]&/@{logllist,genoprob,haploprob};          
        (*to add row and col names to {groupstates, founderhaplo, bvprob, genoprob, haploprob}*)
        groupstates = zygoteGroupStates[ploidy, onlybivalent];
        groupstates2 = StringJoin @@ # & /@ Map[ToString, groupstates, {2}];
        groupstates2 = Transpose[{"genotype" <> ToString[#] & /@ Range[Length[groupstates]],groupstates2}];
        groupstates2 = Join[{{"Genotype", "Code"}}, groupstates2];
        valents2 = Map[Flatten[#, 1] &, valents, {1}];
        valents2 = Map[StringJoin @@ # &, Map[ToString, valents2, {3}], {2}];
        valents2 = StringJoin @@ Riffle[#, "-"] & /@ valents2;
        valents2 = Transpose[{"Valent" <> ToString[#] & /@ Range[Length[valents]], valents2}];
        valents2 = Join[{{"Valent", "Code"}}, valents2];
        siblogl2 = Transpose[{sibID, sibtype, Round[siblogl,10^(-5.)]}];
        siblogl2 = Join[siblogl2,Round[sibtypeprob,10^(-5.)],2];
        temp = Take["Pr("<>#<>")"&/@{"22","24","42","44"},Dimensions[sibtypeprob][[2]]];
        siblogl2 = Prepend[siblogl2, Join[{"SIB", "sibtype","siblogl"},temp]];
        snpID = snpmap[[2 ;;, 1]];
        chrID = snpmap[[2 ;;, 2]];
        haploID = Flatten[Outer[StringJoin, parentID, "_" <> ToString[#] & /@ Range[ploidy]]];
        vvID = "Valent" <> ToString[#] & /@ Range[Dimensions[logllist][[2]]];
        founderhaplo2 = Join[Transpose[{Join[{"SNP", "Chromosome"}, haploID]}], Join[{snpID, chrID}, founderhaplo], 2];
        logllist2 = Join[Transpose[valents2], Join[Transpose[{sibID}], logllist, 2]];
        (*to output genoprob2 and haploprob2*)
        genotypes = "_genotype" <> ToString[#] & /@ Range[Length[groupstates]];
        rowID = Join[{"SNP", "Chromosome"}, Flatten[Outer[StringJoin, sibID, genotypes]]];
        genoprob2 = Flatten[Transpose[#] & /@ genoprob, 1];
        genoprob2 = Join[Transpose[{rowID}], Join[{snpID, chrID}, genoprob2], 2];
        haplotypes = "_haplotype" <> ToString[#] & /@ Range[2 ploidy];
        rowID = Join[{"SNP", "Chromosome"}, Flatten[Outer[StringJoin, sibID, haplotypes]]];
        haploprob2 = Join[Transpose[{rowID}],Join[{snpID, chrID}, Flatten[Transpose[#] & /@ haploprob, 1]], 2];
        (*export*)
        (*, where sibtype = 1,2,3,4 refer to no multivalent, multivalent in the Parent2 gamete, multivalent in the Parent1 gamete, multivalents in both gametes, respectively.*)
        summary = Join[
                     {{key,"Genetic map of biallelic markers"}},snpmap,                     
                     {{key,"MAP of parental haplotypes"}},founderhaplo2,
                     {{key,"Reference haplotypes"}},refhaplo2,                     
                     {{key,"ln marginal likelihood of each valent of each sib"}},logllist2,
                     {{key,"ln marginal likelihood given the LG type of each sib"}},siblogl2,
                     {{key,"Genotypes in order"}},groupstates2, 
                     {{key,"Conditonal genotype probability"}},genoprob2,
                     {{key,"Conditonal haplotype probability"}},haploprob2];
        Export[summaryFile, ExportString[summary, "CSV"], "Table"]
    ]                  


(*http://stackoverflow.com/questions/7525782/import-big-files-arrays-with-mathematica*)
readTable[file_String?FileExistsQ, format_String, chunkSize_: 1000] :=
    Module[ {stream, dataChunk, result, linkedList, add},
        SetAttributes[linkedList, HoldAllComplete];
        add[ll_, value_] :=
            linkedList[ll, value];
        stream = StringToStream[Import[file, "String"]];
        Internal`WithLocalSettings[Null,(*main code*)
            result = linkedList[];
            While[dataChunk =!= {}, 
             dataChunk = ImportString[StringJoin[Riffle[ReadList[stream, "String", chunkSize], "\n"]], format];
             result = add[result, dataChunk];
            ];
            result = Flatten[result, Infinity, linkedList],(*clean-up*)
         Close[stream]
         ];
        Join @@ result
    ]
   
          
getSummaryITO[summaryFile_String?FileExistsQ] :=
    Module[ {res, res2, description, snpmap, esthaplo, refhaplo, logllist, sibtypes, 
        genotypes, estgenoprob, esthaploprob,key = "inferTetraOrigin-Summary"},
        (*res = Import[summaryFile, "CSV"];*)
        res = readTable[summaryFile, "CSV"];
        res2 = Partition[Split[res, #1[[1]] != key && #2[[1]] != key &], 2];
        description = Flatten[#] & /@ res2[[All, 1]];
        res2 = res2[[All, 2]];
        {snpmap, esthaplo, refhaplo, logllist, sibtypes, genotypes} = res2[[{1, 2, 3, 4, 5, 6}]];
        estgenoprob = Partition[res2[[7, 3 ;;, 2 ;;]], Length[genotypes]-1];
        estgenoprob = Transpose[#] & /@ estgenoprob;
        esthaploprob = Partition[res2[[8, 3 ;;, 2 ;;]], Length[esthaplo]-2];
        esthaploprob = Transpose[#] & /@ esthaploprob;
        {description, snpmap, esthaplo, refhaplo, logllist, sibtypes, genotypes, estgenoprob, esthaploprob}
    ]    

indexByInterpolation[data_?(OrderedQ[#] && VectorQ[#, NumericQ] &)] :=
    Module[ {ls, f},
        ls = Transpose[{data, Range[Length[data]] - 1}];
        f = Interpolation[ls, InterpolationOrder -> 0];
        Function[x,
          Round[f[x]] + Switch[Depth[x],
               1, Total[Boole[Thread[Rest[data] == x]]],
               2, Total[Boole[Outer[Equal, x, Rest[data]]], {2}]
              ]
         ]
    ]
    
toGridProb[prob_, snploc_, grid_] :=
    Module[ {ls},
        ls = Append[snploc, 2 snploc[[-1]] - snploc[[1]]];
        ls = indexByInterpolation[ls][grid];
        prob[[All, ls]]
    ]/;Length[prob[[1]]]==Length[snploc]
  
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["TetraOrigin`*"];

EndPackage[]