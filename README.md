# TriC
Scripts for the analysis of Tri-C experiments

The TriC_MO.pl script performs the analysis of a Tri-C experiment. It depends on the bam file output of the CCseqBasic pipeline "COMBINED_reported_capture_reads_CS5.bam" (F6 folder), ran using an oligo file without proximity exclusion (see instructions here: http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/), and converted to sam format. This files contains all mapped reads containing 1 capture fragment post PCR duplicate, blat and ploidy filtering. The scripts also requires the oligo file and a file with the coordinates of the restriction fragments (of correspondong enzyme) in the genome. 

The script selects the reporter reads in cis, performs a proximity exclusion, maps the reads to the restriction fragments and counts multi-way interactions. These are outputted in a text file with suffix "_TriC_interactions.txt". 

Output:
1. Report with read and interaction counts
2. Text file containing reported multi-way interactions per viewpoint in format: RF1 \t RF2 \t count \n
3. Wig file containing all reported interactions in cis per viewpoint (Capture-C like track; optional)

Example of a minimal run command:
nohup perl TriC_MO.pl -sam COMBINED_reported_capture_reads_CS5.sam -o /my_path/tri-c_oligo_file_noprex.txt -r /my_path/mm9_nlaIII_coordinates.txt -name X &

For more instructions, see http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/

The "_TriC_interactions.txt" file can be used as input for the TriC_matrix_MO.py script to generate a basic interaction matrix.
