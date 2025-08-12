//Parameters for the coalescence simulation program : fastsimcoal.exe
6 samples to simulate :
//Population effective sizes (number of genes)
POP0
POP1
POP2
POP3
POP4
POP5
//Samples sizes and samples age
2
2
6
6
6
6
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
1
//Migration matrix 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
17 historical event
11000000 1 0 1 RSANC01 0 0
TDIVNILLT 2 1 1 RSANC12 0 0
TDIVLTLT 4 2 1 LTANC0 0 0
TDIVLT1LT2 3 2 1 LTANC1 0 0
TDIVLT3LT4 5 4 1 LTANC2 0 0
TADM16 1 2 MIG16 1 0 0
TADM16 2 1 MIG61 1 0 0
TADM12 2 1 MIG12 1 0 0
TADM12 1 2 MIG21 1 0 0
TADM14 4 1 MIG14 1 0 0
TADM14 1 4 MIG41 1 0 0
TADM24 2 4 MIG24 1 0 0
TADM24 4 2 MIG42 1 0 0
TADM25 2 5 MIG25 1 0 0
TADM25 5 2 MIG52 1 0 0
TADM34 3 4 MIG34 1 0 0
TADM34 4 3 MIG43 1 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.5e-9 OUTEXP
