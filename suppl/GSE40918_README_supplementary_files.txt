Inferelator_Immgen.txt: 
Inferelator, with default parameters, applied to Immgen dataset from March 2011. Column 1 is the gene name, and each remaining column is for a specific TF with the values being the confidence score in that TF's regulatory interaction for each gene.


pCRMs_5TFs_th17.txt: 
pCRMs computed over ChiP-Seq datasets for BATF, IRF4, STAT3, c-MAF, RORC. Columns are: 
l - pCRM number, sorted based on chr and position
chr - chromosome of pCRM
expt - ChIP peaks comprising this pCRM, ordered by bp position
expt.alphanum.sorted - same as expt column, but sorted alphabetically on TF name
start,end - the start and end of the pCRM; s - median summit of peaks within pCRM
pval - p-value of peaks, sorted as in expt column
pval.mean - mean p-value
span.tfs - span of each peak in pCRM
span.l - span of entire pCRM
peak.ids - peak IDs from MACS file
trgt.prox - genes associated with pCRM (pCRM within +/- 5kb of TSS)
trgt.dist - genes associated with pCRM (pCRM within gene or +/- 10kb of gene)
dtss.prox - distances from pCRM to TSS of genes in column trgt.prox
dtss.dist - distances from pCRM to TSS of genes in column trgt.dist.


DEseq files (named Th17.*.txt): 
Differential expression between wild type and knockouts. Output is directly from the DEseq program. The following datasets were used in each file: 
Th17.batf.wt.vs.Th17.batf.ko_Aug_2_2012.txt - SL3301_SL4222::SL3303_SL4224
Th17.fosl2.wt.vs.Th17.fosl2.ko_Aug_2_2012.txt - SL12528_SL11894::SL12530_SL11896
Th17.hif1a.wt.vs.Th17.hif1a.ko_Aug_2_2012.txt - SL11890_SL12364::SL11892_SL12366
Th17.ikzf3.wt.vs.Th17.ikzf3.ko_Aug_2_2012.txt - SL6494::SL6496
Th17.irf4.wt.vs.Th17.irf4.ko_Aug_2_2012.txt - SL2775_SL3308::SL2777_SL3310
Th17.maf.wt.vs.Th17.maf.ko_Aug_2_2012.txt - SL2763_SL4218::SL2765_SL4220
Th17.rorc.wt.vs.Th17.rorc.ko_Aug_2_2012.txt - SL1844_SL2680::SL1846_SL2682
Th17.siNonTar.24hRPMI.vs.Th17.siEtv6.24hRPMI_Aug_2_2012.txt - SL11677_SL11682::SL11678_SL11683
Th17.siNonTar.24hScreen.vs.Th17.siAes.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13401_SL13402
Th17.siNonTar.24hScreen.vs.Th17.siAtf6.24hScreen_Aug_2_2012.txt - SL13370_SL13374::SL13381_SL13382
Th17.siNonTar.24hScreen.vs.Th17.siBcl11b.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13409_SL13410
Th17.siNonTar.24hScreen.vs.Th17.siCcr6.24hScreen_Aug_2_2012.txt - SL13369_SL13370_SL13371::SL13419_SL13420_SL13421
Th17.siNonTar.24hScreen.vs.Th17.siCrem.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13397_SL13398
Th17.siNonTar.24hScreen.vs.Th17.siDdit3.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13399_SL13400
Th17.siNonTar.24hScreen.vs.Th17.siEgr2.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13407_SL13408
Th17.siNonTar.24hScreen.vs.Th17.siId2.24hScreen_Aug_2_2012.txt - SL13369_SL13374::SL13375_SL13376
Th17.siNonTar.24hScreen.vs.Th17.siInhba.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13403_SL13404
Th17.siNonTar.24hScreen.vs.Th17.siJmjd6.24hScreen_Aug_2_2012.txt - SL13372_SL13373::SL13413_SL13414
Th17.siNonTar.24hScreen.vs.Th17.siKdm6b.24hScreen_Aug_2_2012.txt - SL11677_SL11682::SL11679_SL11684
Th17.siNonTar.24hScreen.vs.Th17.siLef1.24hScreen_Aug_2_2012.txt - SL13369_SL13373::SL13377_SL13378
Th17.siNonTar.24hScreen.vs.Th17.siNfatc2.24hScreen_Aug_2_2012.txt - SL13370_SL13373::SL13393_SL13394
Th17.siNonTar.24hScreen.vs.Th17.siPrkrir.24hScreen_Aug_2_2012.txt - SL13370_SL13374::SL13387_SL13388
Th17.siNonTar.24hScreen.vs.Th17.siRorc.24hScreen_Aug_2_2012.txt - SL13369_SL13370_SL13371_SL13374::SL13415_SL13416_SL13417_SL13418
Th17.siNonTar.24hScreen.vs.Th17.siSatb1.24hScreen_Aug_2_2012.txt - SL13369_SL13374::SL13389_SL13390
Th17.siNonTar.24hScreen.vs.Th17.siSirt2.24hScreen_Aug_2_2012.txt - SL13369_SL13373::SL13385_SL13386
Th17.siNonTar.24hScreen.vs.Th17.siSki.24hScreen_Aug_2_2012.txt - SL13369_SL13374::SL13379_SL13380
Th17.siNonTar.24hScreen.vs.Th17.siTrib3.24hScreen_Aug_2_2012.txt - SL13369_SL13374::SL13383_SL13384
Th17.siNonTar.24hScreen.vs.Th17.siZfp36l2.24hScreen_Aug_2_2012.txt - SL13371_SL13374::SL13405_SL13406
Th17.stat3.wt.vs.Th17.stat3.ko_Aug_2_2012.txt - SL1848_SL3541::SL1850_SL2675.
