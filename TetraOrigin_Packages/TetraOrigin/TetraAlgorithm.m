(* Mathematica Package *)

BeginPackage["TetraOrigin`TetraAlgorithm`",{"ContinuousTimeHmm`","TetraOrigin`TetraModel`"}]

(* Exported symbols added here with SymbolName::usage *)  

Unprotect @@ Names["TetraOrigin`TetraAlgorithm`*"];
ClearAll @@ Names["TetraOrigin`TetraAlgorithm`*"];

tetraMargLiklihood::usage = "tetraMargLiklihood  "

tetraPosteriorDecoding::usage = "tetraPosteriorDecoding  "

tetraPhasing::usage = "tetraPhasing  "

calsiblogl::usage = "calsiblogl  "

tetraMaxPosterior::usage = "tetraMaxPosterior  "

tetraForward::usage = "tetraForward  "

calbvpairlogl::usage = "calbvpairlogl  "

Begin["`Private`"] (* Begin Private Context *) 

tetraMargLiklihood[fhaplo_, startprob_, tranprob_, condstates_] :=
    Module[ {nsib = Length[dataprobset[[1, 1]]],logllist,dataprob,bvss,i,prob,ind},
        bvss = condstates[[All, 1, 2]];
        Monitor[
        logllist = Table[
           dataprob = dataprobset[[All, All, ind]];
           Table[
            prob = MapThread[#1[[#2]] &, {dataprob[[All, All, bvss[[i]]]], fhaplo}];
            CtLogLiklihood[startprob[[i]], tranprob[[i]], prob], {i,Length[bvss]}], {ind, nsib}], ProgressIndicator[ind, {1, nsib}]];
        logllist
    ]   
    
tetraMargLiklihood[startprob_, tranprob_, condstates_] :=
    Module[ {nsib = Length[dataprobset[[1]]],logllist,bvss,i,prob,ind},
        bvss = condstates[[All, 1, 2]];
        Monitor[
        logllist = Table[           
           Table[
            prob = dataprobset[[All, ind, bvss[[i]]]];
            CtLogLiklihood[startprob[[i]], tranprob[[i]], prob], {i,Length[bvss]}], {ind, nsib}], ProgressIndicator[ind, {1, nsib}]];
        logllist
    ]    
                      
tetraPosteriorDecoding[startprob_, tranprob_, condstates_, outputfile_] :=
    Module[ {outputstream,nsib = Length[dataprobset[[1]]],ls,prob,ind,i,bvss},
        Quiet[Close[outputfile]];
        outputstream = OpenAppend[outputfile];
        bvss = condstates[[All, 1, 2]];
        Monitor[
        Do[
            ls = Table[
                prob = dataprobset[[All, ind, bvss[[i]]]];
                CtPosteriorDecoding[startprob[[i]], tranprob[[i]], prob], {i,Length[bvss]}];
            ls[[All,2]] = Round[ls[[All,2]], 10^(-6.)];
            PutAppend[ls,outputstream], {ind, nsib}], ProgressIndicator[ind, {1, nsib}]];
        Close[outputstream];
    ]                  
               
(*fwlogpost[[H=i,Sib=o,OriginState=j]], priorhaploweight[[H=i]]*)
calfwlogpost[fwlogpost_, priorhaploweight_] :=
    Module[ {logscale, sibscale, logpgeno, logpcondorig, loggeno},
        sibscale = Log[Map[Total[Exp[#]] &, fwlogpost, {2}]];
        loggeno = Total[sibscale, {2}] + Log[priorhaploweight];
        logscale = Log[Total[Exp[loggeno]]];
        logpgeno = loggeno - logscale;
        logpcondorig = fwlogpost - sibscale;
        {logscale, sibscale, logpgeno, logpcondorig}
    ]

tetraForward[fhaploweight_,startprob_, tranprob_, dataprob_,bvpair_] :=
    Module[ {nseq = Length[dataprob], ngeno = Length[#] & /@ dataprob, nsib = Length[dataprob[[1, 1]]],
      fwlogpgeno, fwlogpost, fwlogscale, fwsibscale, ind,t,g},
        fwlogpgeno = fwlogpost = fwlogscale = fwsibscale = Table[0, {nseq}];
        (*fwlogpost[[1]] = Log[Map[# startprob &, dataprob[[1]], {2}]];*)
        fwlogpost[[1]] = Log[Table[dataprob[[1, g, ind]]  startprob[[bvpair[[ind]]]], {g,Length[dataprob[[1]]]}, {ind, nsib}]];
        {fwlogscale[[1]], fwsibscale[[1]], fwlogpgeno[[1]], fwlogpost[[1]]} = calfwlogpost[fwlogpost[[1]], fhaploweight[[1]]];
        (*PrintTemporary["Forward calculating..."];*)
        Monitor[Do[
          fwlogpost[[t]] = Log[Total[Table[Exp[fwlogpgeno[[t - 1, g]]] Table[Exp[fwlogpost[[t - 1, g, ind]]].tranprob[[bvpair[[ind]],t - 1]], {ind, nsib}], {g, ngeno[[t - 1]]}]]];
          fwlogpost[[t]] = Table[fwlogpost[[t]], {ngeno[[t]]}] + Log[dataprob[[t]]];
          {fwlogscale[[t]], fwsibscale[[t]], fwlogpgeno[[t]], fwlogpost[[t]]} = calfwlogpost[fwlogpost[[t]], fhaploweight[[t]]], {t, 2, nseq}], 
         ProgressIndicator[t, {2, nseq}]];
        {fwlogpgeno, fwlogpost, fwlogscale, fwsibscale}
    ]


tetraPathSampling[fhaploweight_,startprob_, tranprob_, dataprob_,bvpair_] :=
    Module[ {nseq = Length[dataprob], nsib = Length[dataprob[[1, 1]]], ngeno = 
    Length[#] & /@ dataprob, states = Range[Length[#]] & /@ startprob, fwpgeno, fwprob, 
      fwlogscale, fwsibscale, geno, orig, ls, weight, t,ind},
        {fwpgeno, fwprob, fwlogscale, fwsibscale} = tetraForward[fhaploweight,startprob, tranprob, dataprob,bvpair];
        fwpgeno = Exp[fwpgeno];
        fwprob = Exp[fwprob];
        geno = orig = Table[0, {nseq}];
        geno[[-1]] = RandomChoice[fwpgeno[[-1]] -> Range[ngeno[[-1]]]];
        orig[[-1]] = Table[RandomChoice[fwprob[[-1, geno[[-1]], ind]] -> states[[bvpair[[ind]]]]], {ind, nsib}];
        (*PrintTemporary["Backward sampling..."];*)
        Monitor[Do[
          ls = Table[tranprob[[bvpair[[ind]], t, All, orig[[t + 1, ind]]]], {ind, nsib}];
          ls = ls # & /@ fwprob[[t]];
          weight = Map[Total, ls, {2}];
          weight = fwpgeno[[t]] (Times @@ # & /@ weight);
          weight /= Max[weight];
          geno[[t]] = RandomChoice[weight -> Range[ngeno[[t]]]];
          orig[[t]] = Table[RandomChoice[ls[[geno[[t]], ind]] -> states[[bvpair[[ind]]]]], {ind, nsib}], {t,
            nseq - 1, 1, -1}], ProgressIndicator[t, {1, nseq - 1}]];
        {geno, orig}
    ]
    
calsibloglold[logllist_, ploidy_, onlybivalent_] :=
    Module[ {type, pos, vvprob, ls, sibtype, siblogl,sibpos},
        type = zygoteType[ploidy, onlybivalent];
        If[ Dimensions[logllist][[2]] != Length[type],
            Print["calsiblogl: onlybivalent does not match the dimension of logllist!",{type,Dimensions[logllist]}];
            Abort[]
        ];
        pos = Flatten[Position[type, #, {1}, Heads -> False]] & /@Range[Max[type]];
        vvprob = Exp[# - Max[#]] & /@ logllist;
        ls = Transpose[Total[Transpose[vvprob][[#]]] & /@ pos];
        sibtype = Flatten[Ordering[#, -1] & /@ ls];
        sibpos = pos[[sibtype]];
        siblogl = MapThread[#1[[#2]] &, {logllist, sibpos}];
        (*Mean = equally probable of the vialents given the type*)
        siblogl = Log[Mean[#] & /@ Exp[siblogl]];
        {sibtype/.{1->22,2->24,3->42,4->44}, sibpos,siblogl}
    ]   
    
calsiblogl[logllist_, ploidy_, onlybivalent_] :=
    Module[ {type, pos, prior,logposterior,sibtypeprob,sibtype, siblogl,sibpos},
        type = zygoteType[ploidy, onlybivalent];
        If[ Dimensions[logllist][[2]] != Length[type],
            Print["calsiblogl: onlybivalent does not match the dimension of logllist!",{type,Dimensions[logllist]}];
            Abort[]
        ];
        pos = Flatten[Position[type, #, {1}, Heads -> False]] & /@Range[Max[type]];
        (*equally probable for each of 9 or 16 the chromosome pairing*)
        (*prior(V_0|H)=1/9 for bvModel, or 1/16 for fullModel*)
        prior=1/Dimensions[logllist][[2]];
        logposterior=logllist+Log[prior];
        siblogl=Log[Total[Exp[logposterior],{2}]];
        logposterior-=siblogl;
        sibtypeprob = Transpose[Total[Transpose[Exp[logposterior]][[#]]] & /@ pos];
        sibtype = Flatten[Ordering[#, -1] & /@ sibtypeprob];
        sibpos = pos[[sibtype]];        
        {sibtypeprob,sibtype/.{1->22,2->24,3->42,4->44}, sibpos,siblogl}
    ]       

calbvpairlogl[fhaplo_, startprob_, tranprob_,ploidy_, onlybivalent_] :=
    Module[ {condstates,logllist,sibtype, sibpos,siblogl,sibtypeprob, bvpair,logl},
        condstates = zygoteCondStates[ploidy, onlybivalent];
        logllist = tetraMargLiklihood[fhaplo, startprob, tranprob,condstates];
        {sibtypeprob,sibtype, sibpos,siblogl} = calsiblogl[logllist, ploidy, onlybivalent];
        bvpair = MapThread[RandomChoice[#1[[#2]]->#2]&,{Exp[(#-Max[#])&/@logllist],sibpos}];
        logl = Round[Total[siblogl],10^(-2.)];
        {bvpair,logl}
    ]
    
randomfhaplo[fhaploweight_,bvpair_,startprob_, tranprob_,ploidy_, onlybivalent_] :=
    Module[ {condstates,bvss, dataprob, newbvphase, newfhaplo, neworig, ind,i},
        condstates = zygoteCondStates[ploidy, onlybivalent];
        newbvphase = RandomInteger[{1, Length[condstates[[#]]]}] & /@ bvpair;
        bvss = Extract[condstates[[All, All, 2]],Transpose[{bvpair, newbvphase}]];
        dataprob = Table[dataprobset[[All, All, ind, bvss[[ind]]]], {ind, Length[bvss]}];
        (*dataprob = Transpose[#] & /@ Transpose[dataprob];*)
        dataprob =Transpose[dataprob];
        Do[dataprob[[i]]=Transpose[dataprob[[i]]],{i,Length[dataprob]}];
        {newfhaplo, neworig} = tetraPathSampling[fhaploweight,startprob, tranprob, dataprob,bvpair];
        newfhaplo
    ]
    
tetraMaxPosterior[fhaploweight_,startprob_, tranprob_, ploidy_, onlybivalent_,maxstuck_,maxiteration_] :=
    Module[ {ngeno,fhaplo, bvpair,logl, newbvpair, newfhaplo, newlogl, loglhistory,accept,breakcond,count = 0, it},        
        (*ngeno: number of combinations of possible haplotypes of two parents, given the dosages of two parents;*)
        ngeno = Length[#] & /@ fhaploweight;
        fhaplo = MapThread[RandomChoice[#1 -> Range[#2]] &, {fhaploweight, ngeno}];
        {bvpair,logl} = calbvpairlogl[fhaplo, startprob, tranprob,ploidy, onlybivalent];
        PrintTemporary["iteartion = 0. ln[posterior|parent phases] \[Proportional] " <> ToString[logl]];
        loglhistory = {logl};
        Do[
         newfhaplo = randomfhaplo[fhaploweight,bvpair,startprob, tranprob,ploidy, onlybivalent];
         If[newfhaplo===fhaplo,
         	{newfhaplo,newbvpair,newlogl} = {fhaplo,bvpair,logl},
         	{newbvpair,newlogl} = calbvpairlogl[newfhaplo, startprob, tranprob,ploidy, onlybivalent];
         ];
         AppendTo[loglhistory,newlogl];
         Which[
             newlogl>logl,
             count = 0;
             accept = True;             
             {fhaplo,bvpair,logl} = {newfhaplo,newbvpair,newlogl},
             newlogl==logl,
             breakcond = True;
             accept = True,
             newlogl<logl,
             count++;             
             accept = False;
         ];
         PrintTemporary["iteartion = " <> ToString[it] <>". accept = "<>ToString[accept]<> ". ln[posterior|parent phases] \[Proportional] " <> ToString[logl]];
         If[ breakcond || count==maxstuck,
             Break[]
         ], {it, maxiteration}];
        {fhaplo, bvpair, logl,loglhistory}
    ]     
    
(*fhaploweight: tetraForward --> tetraPathSampling -->randomfhaplo --> tetraMaxPosterior --> tetraPhasing*)    
tetraPhasing[fhaploweight_,startprob_, tranprob_, ploidy_,onlybivalent_,minrun_,maxrun_,maxstuck_,maxiteration_] :=
    Module[ {phaseres,phaselogl,count,run,loglhistory,len},
        phaseres = phaselogl = Table[0,{maxrun}];
        Do[
            PrintTemporary["Phasing algorithm run "<>ToString[run]];
            phaseres[[run]] = tetraMaxPosterior[fhaploweight,startprob, tranprob, ploidy, onlybivalent,maxstuck,maxiteration];
            phaselogl[[run]] = Round[phaseres[[run,3]],10^(-2.)];
            count = Count[phaselogl[[;;run]],Max[phaselogl[[;;run]]]];
            If[ count == minrun,
                phaseres = Take[phaseres,run];
                phaselogl = Take[phaselogl,run];
                Break[]
            ],{run,maxrun}];
        PrintTemporary["Log posterior of phasing runs = "<>ToString[phaselogl]];
        loglhistory = phaseres[[All,-1]];
        len = Length[#]&/@loglhistory;
        If[ Length[len]===maxrun,        	
            Print[Style["Warning: the returned phasing may not be global optimized! Suggest to increase the option maxPhasingRun value.", Red]];
            If[ Mean[Flatten[len]]===maxiteration+1,
                Print[Style["Warning2: suggest to increase the option maxIteration value and/or re-examine the genetic map.", Red]];
            ]
        ];
        phaseres
    ]   
           
End[] (* End Private Context *)

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["TetraOrigin`TetraAlgorithm`*"];

EndPackage[]