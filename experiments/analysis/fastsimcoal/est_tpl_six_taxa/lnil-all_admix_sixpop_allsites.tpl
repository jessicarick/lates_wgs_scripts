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
10 historical event
11000000 1 0 1 RSANC01 0 0
TDIVNILLT 2 1 1 RSANC12 0 0
TDIVLTLT 4 2 1 LTANC0 0 0
TDIVLT1LT2 3 2 1 LTANC1 0 0
TDIVLT3LT4 5 4 1 LTANC2 0 0
TADMNILLT 2 1 MIGNILLT 1 0 0
TADM12 2 1 MIG12 1 0 0
TADM13 3 1 MIG13 1 0 0
TADM14 4 1 MIG14 1 0 0
TADM15 5 1 MIG15 1 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.5e-9 OUTEXP
