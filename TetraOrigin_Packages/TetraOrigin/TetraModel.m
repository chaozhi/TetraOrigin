(* Mathematica Package *)

(* :Context: TetraOrigin`TetraModel`*)
(* :Author: Chaozhi Zheng <chaozhi@gmail.com>*)
(* :Mathematica Version: 9.0.1.0 *)
(* :Description: A package defines teraploid model for haplotype reconstruction.*)

BeginPackage["TetraOrigin`TetraModel`"]

(* Exported symbols added here with SymbolName::usage *)  

Unprotect @@ DeleteCases[Names["TetraOrigin`TetraModel`*"], "dataprobset"];
ClearAll @@ DeleteCases[Names["TetraOrigin`TetraModel`*"], "dataprobset"];

zygoteStates::usage = "zygoteStates  "

zygoteCondStates::usage = "zygoteCondStates  "

zygoteGroupStates::usage = "zygoteGroupStates  "

zygoteCondGroupStates::usage = "zygoteCondGroupStates  "

zygoteType::usage = "zygoteType  "

tetraPriorProcess::usage = "tetraPriorProcess  "

tetraDoseLikelihood::usage = "tetraDoseLikelihood  "

tetraDoseLikelihood2::usage = "tetraDoseLikelihood2  "

dataprobset::usage = "dataprobset "

calfhaploset::usage = "calfhaploset  "

calDataProb::usage = "calDataProb  "


Begin["`Private`"] (* Begin Private Context *) 

(******************************************PRIOR*******************************************)

zygoteStates[ploidy_, onlybivalent_] :=
    Module[ {s},
        s = If[ onlybivalent,
                Permutations[Range[ploidy], {2}],
                Tuples[{Range[ploidy], Range[ploidy]}]
            ];
        Flatten[#] & /@ Tuples[{s, s + ploidy}]
    ]

zygoteGroupStates[ploidy_, onlybivalent_] :=
    Module[ {s, i, j},
        s = Flatten[Table[{i, j}, {i, ploidy}, {j, i + Boole[onlybivalent] , ploidy}], 1];
        Flatten[#] & /@ Tuples[{s, s + ploidy}]
    ]

(*Only the largest mutivalent is included, so it works only for \
ploidy=4*)
gametevalents[ploidy_, onlybivalent_] :=
    Module[ {set, bivalents},
        set = Range[ploidy];
        bivalents = Subsets[Subsets[set, {2}], {ploidy/2}];
        bivalents = Pick[bivalents, (Union[Flatten[#]] == Range[ploidy] & /@bivalents)];
        If[ onlybivalent,
            bivalents,
            Join[bivalents, {{set}}]
        ]
    ]

zygotevalents[ploidy_, onlybivalent_] :=
    Module[ {s},
        s = gametevalents[ploidy, onlybivalent];
        Tuples[{s, s + ploidy }]
    ]

valent2gamete[x : {{_Integer, _Integer}, {_Integer, _Integer} ..}] :=
    Tuples[x]
    
valent2gamete[x : {{_Integer ..}}] :=
    Tuples[First[x], 2]

zygoteCondStates[ploidy_, onlybivalent_] :=
    Module[ {states, valents, allvalents, ls},
        states = zygoteStates[ploidy, onlybivalent];
        valents = zygotevalents[ploidy, onlybivalent];
        allvalents = Tuples[#] & /@ Map[Permutations, valents, {2}];
        ls = Map[valent2gamete, allvalents, {3}];
        ls = Map[Flatten, Map[Tuples, ls, {2}], {3}];
        ls = Replace[ls, Thread[states -> Range[Length[states]]], {3}];
        MapThread[Transpose[{#1, #2}] &, {allvalents, ls}]
    ]    

zygoteCondGroupStates[ploidy_, onlybivalent_] :=
    Module[ {groupstates, valents, ls},
        groupstates = zygoteGroupStates[ploidy, onlybivalent];
        valents = zygotevalents[ploidy, onlybivalent];
        ls = Map[valent2gamete, valents, {2}];
        ls = Map[Flatten, Tuples[#] & /@ ls, {2}];
        ls = Union[#] & /@ Map[Sort, ls, {2}];
        ls = Replace[ls, Thread[groupstates -> Range[Length[groupstates]]], {2}];
        Transpose[{valents, ls}]
    ]

zygoteType[ploidy_, onlybivalent_] :=
    Module[ {condstates, valents, nvalents, rule},
        condstates = zygoteCondStates[ploidy, onlybivalent];
        valents = condstates[[All, 1, 1]];
        nvalents = Replace[Map[Length, valents, {3}], {ploidy} -> {ploidy, ploidy}, {2}];
        nvalents = nvalents[[All, All, 1]];
        rule = {{2, 2} -> 1, {2, ploidy} -> 2, {ploidy, 2} ->3, {ploidy, ploidy} -> 4};
        nvalents /. rule
    ]
  
(*r=recombination fraction between two loci*)
valentTranMtx[r_, nvalent_] :=
    r/(nvalent - 1) + (1 - r (1 + 1/(nvalent - 1))) IdentityMatrix[nvalent]

gameteTranMtx[r_, nvalents_] :=
    KroneckerProduct @@ (valentTranMtx[r, #] & /@ nvalents)

zygoteTranMtx[r_, nvalentsP1_, nvalentsP2_] :=
    KroneckerProduct[gameteTranMtx[r, nvalentsP1], gameteTranMtx[r, nvalentsP2]]

(*return e.g.for 75 markers of a given linkagegroup and onlybivalent =False
Length[#] & /@ startprob:
    {16, 16, 16, 64, 16, 16, 16, 64, 16, 16, 16, 64, 64, 64, 64, 256}
Dimensions[#] & /@ tranprob: 
      {{74, 16,16}, {74, 16, 16}, {74, 16, 16}, {74, 64, 64},
       {74, 16,16},{74, 16, 16}, {74, 16, 16}, {74, 64, 64},{74, 16, 16}, {74,16, 16},{74, 16, 16}, 
       {74, 64, 64}, {74, 64, 64}, {74, 64,64}, {74, 64, 64}, {74, 256, 256}}
*)
tetraPriorProcess[snpseq_, linkagegroup_, ploidy_, onlybivalent_] :=
    Module[ {rfseq, condgroupstates, valents, nvalents, fun, pos, dim, nv,tranprob, startprob},
        (*genetic distance in Morgan*)
        rfseq = Differences[snpseq[[linkagegroup, All, -1]]]/100;
        (*Haldane map function, here rfreq defined one chromsome pairing, 
        so the map function is the same as the diploid case*)
        rfseq = 1/2 (1 - Exp[-2 rfseq]);
        condgroupstates = zygoteCondGroupStates[ploidy, onlybivalent];
        valents = condgroupstates[[All, 1]];
        nvalents = Replace[Map[Length, valents, {3}], {ploidy} -> {ploidy, ploidy}, {2}];
        tranprob = ConstantArray[0, Length[nvalents]];
        Do[
         fun = Function[{r}, Evaluate[zygoteTranMtx[r, Sequence @@ nv]]];
         pos = Flatten[Position[nvalents, nv, {1}, Heads -> False]];
         tranprob[[pos]] = Table[Map[fun, rfseq], {Length[pos]}], {nv, Union[nvalents]}];
        dim = Length[#] & /@ tranprob[[All, 1]];
        startprob = ConstantArray[1./#, #] & /@ dim;
        {startprob, tranprob}
    ] 
    
(********************************DATA Likelihood*******************************************)
(*mapping dosages to phased genotype*)
ruleToGeno[ploidy_] :=
    Module[ {temp, geno, i, j,rules},
        rules = Table[temp = Subsets[Range[ploidy], {i}];
                      geno = ConstantArray[1, {Length[temp], ploidy}];
                      Do[geno[[j, temp[[j]]]] = 2, {j, Length[temp]}];
                      i -> geno, {i, 0, ploidy}];
        Join[rules, {"NA" -> Flatten[rules[[All, 2]], 1]}]
    ]


(*for a particular linkage group*)
calfhaploset[fdose_, ploidy_] :=
    Module[ {genorule, fhaploset, fhaploweight},
        genorule = ruleToGeno[ploidy];
        fhaploset = Replace[fdose, genorule, {2}];
        fhaploset = 
         Map[Flatten[Outer[List, #[[1]], #[[2]], 1], 1] &, fhaploset];
        fhaploset = Map[Flatten, fhaploset, {2}];
        fhaploweight = Table[1/Length[#], {Length[#]}] & /@ fhaploset;
        {fhaploset, fhaploweight}
    ]

calfhaploset[fdose_, ploidy_, epsF_] :=
    Module[ {set, fhaploset, dose, freq, weight, fhaploweight,i,k},
        set = "NA" /. ruleToGeno[ploidy];
        set = Flatten[Outer[List, set, set, 1], 1];
        fhaploset = Table[set, {Length[fdose]}];
        dose = Map[Total, fhaploset - 1, {3}];
        (*#haplotypes correspond to dosages 0,1,2,3,4*)
        freq = Table[Binomial[ploidy, i], {i, 0, ploidy}];
        fhaploweight = Table[
          weight = Transpose[dose[[i]]];
          Do[If[ fdose[[i, k]] === "NA",
                 weight[[k]] = 1/((ploidy + 1) freq[[dose[[i, All, k]] + 1]]),
                 weight[[k]] = (weight[[k]] /. {fdose[[i, k]] -> 1 - epsF, _?(# != fdose[[i, k]] &) :> epsF/ploidy});
                 weight[[k]] /= freq[[dose[[i, All, k]] + 1]]
             ], {k, Length[weight]}];
          Times @@ weight, {i, Length[dose]}];
        {Map[Flatten, fhaploset, {2}], fhaploweight}
    ]

(*for a particular linkage group*)
calDataProb[founderdata_,offdose_,  eps_, ploidy_, onlybivalent_] :=
    Module[ {states, ss, dgeno, ddose,dataprob, offdose2,i, j, ls},
     (*dgeno:Genotype derived from states of parental origins and founder haplotype;
      ddose:the dosages based on the derived genotype (dgeno);
      dataprob: the probability of obsvered dosage at a marker for all siblings;*)
        states = zygoteStates[ploidy, onlybivalent];
        ss = Transpose[states];
        dgeno = Switch[Depth[founderdata],
            4,
            (*founderdata=fhaploset, list of all possible haplotypes at each markers*)            
            Table[Transpose[founderdata[[i, j, #]] & /@ ss], {i, Length[founderdata]}, {j, Length[founderdata[[i]]]}],
            3,
            (*founderdata=founderhaplo, list of haplotypes at each markers*)
            Table[Transpose[founderdata[[i, #]] & /@ ss], {i, Length[founderdata]}],
            _,
            Print["Wrong input of founderdata in calDataProb!"];
            Abort[]
        ];
        ddose = Map[Total, dgeno - 1, {Depth[dgeno]-2}];
        offdose2 = offdose/.{"NA" ->Infinity};
        dataprob = Monitor[
            Table[ls = (ddose[[i]] - #) & /@ offdose2[[i]];
                  ls = Replace[ls, {0 -> 1 - eps, -Infinity -> 1, _ -> eps/ploidy}, {Length[Dimensions[ls]]}];
                  ls, {i, Length[ddose]}], ProgressIndicator[i, {1, Length[ddose]}]];
        If[ Depth[founderdata]==4,
            Do[dataprob[[i]]=Transpose[dataprob[[i]]],{i,Length[dataprob]}]
        ];
        dataprob
    ]
        
(*tetraDoseLikelihood returns fhaploset,fhaploweight, and the global variable dataprobset.
e.g. 75 markers and 50 sibs
Dimensions[#] & /@ fhaploweight: 
    {64,96,96,...}
Dimensions[#] & /@ fhaploset: 
    {{64, 8}, {96, 8}, {96, 8}, {6, 8}, {36, 8}, {64, 8}, {256, 8}, {64,8}, {256, 8},...}
Dimensions[#] & /@ dataprobset: 
    {{64, 50, 256}, {96, 50, 256}, {96, 50, 256}, {6, 50, 256}, {36, 50, 256}, 
     {64, 50, 256}, {256, 50, 256}, {64, 50, 256}, {256, 50,256}, ...}
*)

tetraDoseLikelihood[founderdose_, sibdose_, linkagegroup_, eps_, epsF_,ploidy_, onlybivalent_] :=
    Module[ {fdose, offdose, fhaploset,fhaploweight},
        ClearAll[dataprobset];
        offdose = Transpose[sibdose[[All, linkagegroup]]];
        (*Dimensions[founderdose] == ReplacePart[Dimensions[sibdose], 1 -> 2]*)
        fdose = Transpose[founderdose[[All, linkagegroup]]];
        {fhaploset, fhaploweight} = If[epsF==0,calfhaploset[fdose, ploidy],calfhaploset[fdose, ploidy,epsF]];
        dataprobset = calDataProb[fhaploset, offdose, eps, ploidy, onlybivalent];
        {fhaploset,fhaploweight}
    ]
    
tetraDoseLikelihood2[founderhaplo_, sibdose_, linkagegroup_, eps_, ploidy_, onlybivalent_] :=
    Module[ {offdose},
        ClearAll[dataprobset];
        offdose = Transpose[sibdose[[All, linkagegroup]]];
        dataprobset = calDataProb[founderhaplo[[linkagegroup]], offdose, eps, ploidy, onlybivalent];
    ]

End[] (* End Private Context *)

SetAttributes[#, {Protected,ReadProtected}]&/@ DeleteCases[Names["TetraOrigin`TetraModel`*"], "dataprobset"];

EndPackage[]