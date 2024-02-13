#Script to automate identification of contaminant reads for any virus based on read number and clustering
#Emma Thomson 12th February 2024 

#1. File format is as a fasta file separated by virus of interest including only sequences with a certain percentage of coverage, based on SAM_STATS file
#It is best to run a mapping pipeline for this so that all the files are ready in the correct format
#The tree will be available following the mapping pipeline
#2. Second step is to run cluster-picker using the fasta and tree as input
#3. Using the cluster-picker output, select the best parameters

#Dependencies include: SAM_STATS and SAM2COVERAGE scripts by Sreenu Vattipally, ClusterPicker by Sam Lycett, ClusterPicker,Java, and IQTree software

clear

echo "
CLUSTER-PICKER ANALYSIS SET-UP
This is a tidying script that removes obvious cross-contaminants or mis-aligned sequences from a fasta file made from your sam files. It doesn't replace good laboratory practice and precautions! The setup requires you to use a mapping pipeline with the virus of interest prior to running this script as you will need several input files from this process, including:
-a fasta file of consensus sequences
-a samstats file (produced by mapping pipeline)
-an allconfile (produced by mapping pipeline)
-a genotype check file where you have checked that the mapped reference matches the phylogenetic reference. One will be prepared for you if you don't have one  

There are sevieral dependencies for this script. You need to have clusterpicker, java and iqtree installed
"

#Create some necessary variables
read -p "What coverage threshold do you want to use? [50]" coverage_threshold
coverage=${coverage:-50}
read -p "What is the name of the virus you are interested in (add segment also if appropriate)?" virus
ls *sam
read -p "What is the name of the accession number of the virus? (Make sure that your sam files contain this accession number or the script will fail)." accession

#Define cluster picker settings
echo "Cluster picker settings help you to define the size of clusters you want to use."
#Defaults are inserted for three of the cluster picker variables. It would be possible to adjust these if needed. See cluster picker literature.
cp1=0
cp2=0
read -p "What genetic distance threshold for the clusters do you want to use? Enter a cluster size that you want to test e.g. 0.05 or more/less (but not 0!) [0.05]: " cp3
cp3=${cp3:-0.05}
cp4=0

#Define what % of reads you want as the read number filter (% of highest number of reads in cluster)
read -p "What setting do you want to use to filter the reads within clusters? If you want to remove all reads in a cluster that are less than the highest number of reads in that cluter, select 100. 10 (meaning 10% of the highest number of mapped reads in that cluster) is often a good setting here for metagenomic and TE runs.[10]: " filter_reads
filter_reads=${filter_reads:-10} 

gc1=gc1


read -p "Do you want to set a read threashold for your samples e.g. at least 500 mapped reads per sample? It's a good idea to base this on the mean number of reads in any water samples [0]: " read_threshold
read_threshold=${read_threshold:-0}

#Run sam stats
for st in *"$accession"*.sam; do echo $st>>"$accession"-samstats; SAM_STATS $st>>"$accession"-samstats; done
for st in *"$accession"*.sam; do SAM2CONSENSUS -i $st> $st.fa; done

#Create the sam_con file that you will need
rm -rf $accession.samstats
for samstats in "$accession"-samstats; do grep -B 1 "Coverage" $samstats > $accession.samstats; done
rm -rf sam_cov 
grep "sam" $accession.samstats > sam_temp
grep "Cov" $accession.samstats > cov_temp
paste sam_temp cov_temp > sam_cov

# Define the output file names
rm -rf "$virus"-"$accession"-"$coverage_threshold"-all.fa

# Process the sam_cov file and concatenate the files into a new file
sed -i 's/Coverage://g' sam_cov
sed -i 's/\%//g' sam_cov
cat sam_cov | awk -F'\t' -v threshold="$coverage_threshold" '{
    if ($4 + 0 > threshold) {
        print $1 ".fa";
    }
}' | while read -r filename; do
    if [ -f "$filename" ]; then
        cat "$filename" >> "$virus"-"$accession"-"$coverage_threshold".fa
    else
        echo "File $filename not found." >&2
    fi
done

allconfile="$virus"-"$accession"-"$coverage_threshold".fa
echo "Concatenation complete. Result is in $allconfile"
sed -i 's/Cannot open output\. Sending the output to STDOUT//g' $allconfile
sed -i 's/There_are_no_mapped_reads//g' $allconfile
sed -i '/^$/d' $allconfile
sed -i 's/ /_/g' $allconfile
rm -rf sam_temp cov_temp sam_cov

#read -p "What is the name of the run you are analysing. NB this script version will only recognise 1 run [1]: " run
run=${run:-1}


echo "You are working on "$virus", have selected a read cut-off of "$filter_reads"%, coverage of "$coverage"% and cluster size of "$cp3". You also selected a read threshold of "$read_threshold" reads."

echo "Contamination report log file will contain your results at the end of the run in the following file: cross-contamination-report-""$$"

echo "The analysis should take around 10 minutes - a bit longer if you have lots of samples. Don't disconnect and wait for further prompts in about 2 minutes..."

sed -i 's/ /_/g' $allconfile

#Logging starting number of sequences in cross-contamination report
echo "CROSS-CONTAMINATION REPORT" > cross-contamination-report-$$
acf=$(grep ">" $allconfile| wc -l)
cs=$(grep ">" $allconfile | wc -l)
echo "Allconfile is " "$allconfile"

#Put the sequences on one line
rm -rf $allconfile.one
cat $allconfile | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'>>$allconfile.one

#Remove empty fasta sequences from the allconfile
rm -rf newconfile
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' $allconfile.one  > newconfile

mafft --auto --thread -8 newconfile > $allconfile.fa_aligned
fasta=$allconfile.fa_aligned

#Fasttree is run to make the phyogeny
echo "Running fasttree..."; fasttree $fasta > $fasta\.treefile 

##THE 3 LINES BELOW CAN BE UNBLOCKED TO USE A PRE-EXISTING TREE
#ls *treefile *tree
#read -p "Do you have an underlying tree already? Insert the name of the tree if you have one. Otherwise press returun [$fasta\.treefile] " tree
tree=${tree:-$fasta\.treefile}

#Second step is to run cluster picker
rm -rf *clusterPicks*
echo "Running cluster-picker v.1.2.5..."
java -jar /software/cluster-picker-v1.2.5/release/ClusterPicker_1.2.5.jar $fasta $tree $cp1 $cp2 $cp3 $cp4
mv $fasta\_clusterPicks_list.txt $fasta\_clusterPicks_$cp3\_list.txt
mv $fasta\_$fasta\_clusterPicks.fas $fasta\_clusterPicks_$cp3.fas 
mv $fasta\_clusterPicks.nwk.figTree $fasta\_clusterPicks_$cp3.nwk.figTree

cp_no=$(grep ">" $fasta\_clusterPicks_$cp3.fas | grep -v "$ref" |wc -l)
cluster_no=$(grep ">" $fasta\_clusterPicks_$cp3.fas | cut -f 1 -d "_" | sort | uniq |wc -l)
single_count=$(cat $fasta\_clusterPicks_$cp3\_list.txt | awk '{if ($2==-1) print $1}'|grep -v "$ref" | grep -v "SequenceName" |wc -l)

echo ""$cp_no" clustered sequences were detected in the "$fasta"_clusterPicks_"$cp3".fas file in "$cluster_no" clusters. 
There were "$single_count" singletons. 
Clustering results can be viewed in the "$fasta"_clusterPicks_"$cp3".nwk.figTree file.  This tree file shows you in different colours how the clusters were picked for the cluster filter analysis. 
If the cluster size looks like it should be altered, you can rerun and adjust the clustering variable in the setup." >> cross-contamination-report-$$

rename 's/__/_/g' *

#Next step is to create a tab-del file with essential elements based on the fasta names- for trouble-shooting, check out the tab delimited input file

echo "Defining variables"
#Firsly pull out the names

rm -rf names*

#Get rid of reference names and pull out names of samples

grep ">" $fasta\_clusterPicks_$cp3.fas | cut -f 2 -d ">" | sort | uniq > names_$cp3

names_count=$(cat names_$cp3 | wc -l)
echo ""$names_count" sequences occurred within clusters and were then extracted from the "$fasta"_clusterPicks_"$cp3".fas file" >> cross-contamination-report-$$

cut -f 1 -d "_" names_$cp3 | sed 's/Clust//g'> clusters_$cp3

#This is specific to this study - i.e. the names begin with EBB - and then you list the 8 characters after EBB (the 11 adds 3 for EBB also). Needs to be adapted to other setups/filenames

awk '/EBB/{print substr($0, index($0, "EBB"), 11)}' names_$cp3 | sed 's/_/-/g'> patient_$cp3

#List the number of mapped reads for each file (it is appended at the end of every sam filename)
cat names_$cp3 | rev | cut -f 1 -d "_" | rev > mapped_$cp3

echo "Making tab delimited files"
rm -rf seq_no*
seq_no1=$(cat names_$cp3 | wc -l) 
for i in `seq $seq_no1`; do echo "$cp3">>seq_no1; done

paste patient_$cp3 mapped_$cp3 clusters_$cp3 | awk -v FIL0="$read_threshold" '{if ($2 >= FIL0) {print $1,$2,$3,$4="PASS"} else { print $1,$2,$3,$4="FAIL"}}' | grep "PASS" | awk '{print $1,$2,$3,$4}'> Table1_Total_Mapping_Reads_Filter-$cp3-$accession-$virus

#Count how many sequences removed due to read threshold
paste patient_$cp3 mapped_$cp3 clusters_$cp3 | awk -v FIL0="$read_threshold" '{if ($2 >= FIL0) {print $1,$2,$3,$4="PASS"} else { print $1,$2,$3,$4="FAIL"}}' | grep "FAIL" | awk '{print $1,$2,$3,$4}' > filter0_fail
f0=$(cat filter0_fail | wc -l)

echo "
FILTER 1. The Total Mapping Reads Filter removes any sequences  where the total number of reads is less than your pre-defined amount, a read threshold of "$read_threshold" reads" >> cross-contamination-report-$$

#Record of full input table - including clusters
tdf_no=$(cat Table1_Total_Mapping_Reads_Filter-$cp3-$accession-$virus | wc -l)
echo "You can check the tab delimited file which is called Tab_delimited_input_file_"$virus" to see the input file for sorting later. There are "$tdf_no" lines in the file." >> cross-contamination-report-$$

echo "FILTER 2 - sorting singletons..."
#FILTER 2 (singletons)- samples that form a singleton cluster are designated as -1 on cluster picker
#These are considered to be good but probably best to use a conservative (i.e. larger) cluster size

cat $fasta\_clusterPicks_$cp3\_list.txt | awk '{if ($2==-1) print $1}' | grep -v "SequenceName" |sed 's/^/>/g'> singles
cat singles | rev | cut -f 1 -d "_" | rev > singles_2
cat singles | rev | cut -f 2,3,4,5,6 -d "_" | rev > singles_1
paste singles_1 singles_2 > singles_tab

#Need to remove those below the read threshold
cat singles_tab | awk -v FIL0="$read_threshold" '{if ($2 >= FIL0) {print $1,$2,$3="PASS"} else { print $1,$2,$3="FAIL"}}' | grep "PASS" | awk '{print $1,$2}' | sed 's/ /_/g' > singletons_filter2_pass
rm -rf $virus\_singletons_filter2_pass_$cp3\.fasta
while read a; do grep -A 1 "$a" $allconfile.one ; done < singletons_filter2_pass >>$virus\_singletons_filters_1-2_pass_$accession-$cp3\.fasta
f2=$(grep ">" $virus\_singletons_filters_1-2_pass_$accession-$cp3\.fasta | wc -l)
echo "
"$f2" sequences passed to the singletons file to be added at the end of filtering to the final file - these also passed the read threshold of "$read_threshold" reads" >> cross-contamination-report-$$

#FILTER 3 - MAX READ NUMBER IN CLUSTER
#Probably the most important one
echo "Running FILTER 3 - max read numbers in clusters"
#htf stands for hcv table full
cat Table1_Total_Mapping_Reads_Filter-$cp3-$accession-$virus | sort  -V -k 3 -V -k 2  | sed 's/ /\t/g'> full_table_sorted_$cp3
rm -rf pnt_id2_reads_run*_cluster*

while read a b c d; do echo "$a"",""$b"",""$c" >>pnt_reads_cluster_$cp3; done < full_table_sorted_$cp3

echo "Calculating read numbers and percentage in clusters"
#Calculate the highest number of reads for each cluster and make a table

# Sort the file by cluster (numerically) and then by readcount (numerically in reverse)
sort -t, -k3,3n -k2,2nr pnt_reads_cluster_$cp3 > pnt_reads_cluster_sorted_$cp3

# Use awk to calculate the maximum read count for each cluster and append it as a 4th column
awk -F, '
{
    if (NR == 1 || $3 != prevCluster) {
        maxReadCount = $2; # Reset maxReadCount for the new cluster
    }
    if ($2 > maxReadCount) {
        maxReadCount = $2; # Update maxReadCount if the current readcount is higher
    }
    print $0","maxReadCount; # Print the line with maxReadCount appended
    prevCluster = $3; # Update prevCluster for the next iteration
}' pnt_reads_cluster_sorted_$cp3 | sed 's/,/\t/g' > pnt_reads_cluster_sorted_max_$cp3

# Clean up the intermediate sorted file
rm -rf pnt_reads_cluster_sorted_$cp3

echo "
Filter 3 calculates the number of reads for each sequence within each cluster. Then it calculates the maximum number of reads for one sequence within every cluster. This is used as a denominator to calculate a threshold for removing reads within the same cluster under the likely level for cross-contamination. 
The clustering pairwise distance in this run was "$cp3" and the threshold was "$filter_reads"%" >> cross-contamination-report-$$

#The join or awk command adds in the max reads to the original table to allow the % calculation to be carried out - the files need to be sorted first

#This calculation gives the % of the highest number of reads in the cluster and puts it in column 5
awk -v OFS='\t' '{$5 = ($2/$4)*100}1' pnt_reads_cluster_sorted_max_$cp3 > Table_with_max_read_counts_$cp3
mrc=$(cat Table_with_max_read_counts_$cp3 | wc -l)
echo "
The input table for calculating the maximum number of reads in each cluster can be found in the Table_with_max_read_counts_"$cp3" file which has "$mrc" rows. If you think there may have been an error, you can take a look at this file more closely to trouble-shoot. " >> cross-contamination-report-$$
#rm -rf filter3_reads* f3_reads*
echo "Filtering with a read cut-off of " $filter_reads"%"
#Need to tell awk about the $filter_reads variable
rm -rf f3_reads_$cp3
cat Table_with_max_read_counts_$cp3 |  sed 's/\t/ /g' |cut -f 1,2,3,4,5 -d " " | awk -v FIL3="$filter_reads" '{if ($5 >= FIL3) { print $1,$4,$2,$3,$6="PASS"; } else { print $1,$4,$2,$3,$6="FAIL";  } }' | grep "PASS"| sed 's/_/\t/g'| sed 's/-/_/g'| sed 's/ /\t/g' >>f3_reads_$cp3

f3=$(cat f3_reads_$cp3 | wc -l)

while read a b c d e; do grep -A 1 "$a" $allconfile.one ; done < f3_reads_$cp3 >>$virus\_filters_1-3_pass_$accession-$cp3\.fasta

echo "
Filter 3 resulted in "$f3" sequences passing filtration. 
These can be found in the $virus\_singletons_filters_3_pass_$accession-$cp3\.fasta file" >> cross-contamination-report-$$

#Final results
cat $virus\_filters_1-3_pass_$accession-$cp3\.fasta $virus\_singletons_filters_1-2_pass_$accession-$cp3\.fasta > $virus\_filters_1-2-3_pass_$accession-$cp3\.fasta

echo "Final data can be found in the "$virus"_filters_1-2-3_pass_"$accession-$cp3".fasta file"
echo "Final data can be found in the "$virus"_filters_1-2-3_pass_"$accession-$cp3".fasta file" >>cross-contamination-report-$$

rm -rf filter* cluster0* multi* hcv_tab_full* full_analysis_table* filter5_reads* No_geno* Good_geno* Bad_geno* cluster_run* good* bad_* accession* singletons_filter2_pass subgeno* id2*  seq_no* accession* f4* run*_* clusters_* geno_mis* names* mapping_geno phylo_geno geno_check geno_mismatches* gc1 singles* htf_* patient_cluster_* pnt_cluster* newconfile pnt* patient* mapped_* *samstats *clusterPicks*list.txt *clusterPicks*nwk *clusterPicks*log.txt *clusterPicks*.fas *_aligned full_table_sorted* f3_reads* *.one Table1_Total_Mapping* Table_with_max_read*
