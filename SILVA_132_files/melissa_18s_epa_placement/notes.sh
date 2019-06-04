/home/morien/bin/sina-1.2.11/sina -i NODE_REP_FORDOWNSTREAM.fasta --intype fasta -o NODE_REP_FORDOWNSTREAM_18s.aligned.fasta --outtype fasta --ptdb /parfreylab/morien/taxonomy_DBs/silva_128_18s_SINA_db/SSURef_NR99_128_SILVA_07_09_16_opt.arb --overhang remove --insertion forbid

#combined alignment with 97 rep set
cat ~/Desktop/lab_member_files/taxonomy_databases/SILVA128/18s/rep_set_aligned/97_SILVA_128_aligned_rep_set.plus_entamoeba_blastocystis.fasta NODE_REP_FORDOWNSTREAM_18s.aligned.cleaned.fasta > full_18s_seq_set.aligned.fasta

#filtered for gaps and highly entropic bases
filter_alignment.py -i full_18s_seq_set.aligned.fasta -o ./ -s -e 0.05 -g 0.95

#raxml EPA placement
/home/morien/bin/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f v -p 5260 -m GTRCAT -t /parfreylab/shared/databases/SILVA_128/18s/trees/SILVA_128_18s_97_rep_set.tre -s /parfreylab/morien/epa_placement_18s/full_18s_seq_set.aligned_pfiltered.fasta -n EPA_placement_18s -w /parfreylab/morien/epa_placement_18s/raxml_epa/ -T 6

#fix tree for figtree viewing
sed -e 's/\[I[[:digit:]]*\]//g' raxml_epa/RAxML_labelledTree.EPA_placement_18s > raxml_epa/RAxML_labelledTree.EPA_placement_18s.figtree.tre

#OPTIONAL: remove QUERY label from placed sequences (ony do if necessary, these labels are useful!)
sed -e 's/QUERY___//g' input.tre > iput.no_query_labels.tre