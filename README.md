# lheToTree
Quick documentation on writeLHEToTree.cc and reweightTree.cc :


0) Usage: both files need to be compiled using CMSSW_7_4_11 and adding the root libraries glibs i.e.

g++ writeLHEToTree.cc  -o writeLHEToTree `root-config --cflags --glibs`
to compile reweightTree the libraries --libs --evelibs have to be added

the headerfiles AllHadronicConstants.h  Eff_lvjj.h  WWConstants.h  emEff.h  hadEff.h  had_ZZ_jjjj_Eff.h are needed

0.1) check pdgId of your graviton candidate: if it is different from 39 change the variable pdgIdGraviton, then compile again

1) after producing the .lhe file, use writeLHEToTree with three parameters:

    1.1) outputdirectory/outputfilename
    1.2) decay mode: WW_lvjj or WZ_lvjj or WW_jjjj or ZZ_jjjj
    1.3) inputdirectory
    

1.1) this input paramter is the name of the directory and output root file.

1.2) the decay mode paramter changes which tree is going to be written. Depending on the decay mode different variables get written in the root file: 
If WW_lvjj :
    Wlep/Whad pt eta phi e
    lepton pt eta phi e
    neutrino pt eta phi e
    tau pt eta phi e
    leptonID = 11 or 13
    isTauEvent = True or False
    neve = event-number
    mWW = invariant mass of Wlep+Whad 
    
If WZ_lvjj :
   writes same variables as in WW_lvjj only this time Whad is actually the hadronic Z
   
If WW_jjjj :
   Wplus/Wminus pt eta phi e
   neve = event-number
   mWW = invariant mass of Wplus+Wminus
   
If WZ_jjjj :
   Wplus/Wminus pt eta phi e
   neve = event-number
   mWW = invariant mass of Wplus+Wminus
   this time the Z  gets written in variable named Wplus

   
If ZZ_jjjj :
   Z1/Z2 pt eta phi e
   neve = event-number
   mWW = invariant mass of Z1+Z2
   
   
1.3)  this has to be the name of the directory in which the input .lhe file lies. If the name of the lhe file is different from cmsgrid_final.lhe this has to be changed in the code

2)Use reweightTree with three arguments:

    2.1) treedirectory/treename.root
    2.2) decay mode: use same decay mode with which the tree was produced
    2.3) polarisation : if in the model the W/Z are transversally polarised, this string must contain "trans" 
    
    
reweightTree adds to new branches to the tree: 
    passedAcceptance
    eventWeight
    
to find the acceptance/efficiency just loop over the tree and weight each event with eventweight (for the efficiencie numerator)
each event with passedAcceptance = true has passed the acceptance cuts of the analysis
