#!/bin/bash
multi_rna\
	--mapfile ./mapfile.txt\
	--genomeDir /home/groups/singlecell/mabdalfttah/CeleScope/hs_ensembl_99\
	--thread 8\
	--mod shell\
        --expected_cell_num 5000 
