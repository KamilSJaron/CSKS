# generate k-mer spectrum according to genome and population parameters
echo "########################################"
echo "COALESCENT SIMULATION FOR K-MER SPRECTRA"

# makes sure the script won't proceeds if any of the commands returns a non-zero exit status (i.e. an error)
set -e

# activate conda environment (run 00setUpEnv.sh to set up)
conda activate CSKS

# make a temp dir
td=$(mktemp -d)

# simulation parameters
gs=1000000 # Size in nt of each chromosome. There are 10 chromosomes, so the overall haploid genome size is 'gs'*10!
ploy=4 # ploidy level
theta=0.005 # population-scaled mutation rate (= heterozygosity in a random-mating population)
#gc=0.5 #  genome GC-content (ignored for now)
#div=0.05 # genome divergence (allotetraploids only)
kmax=5000 # the max coverage that is counted in the kmer database / histogram (there is no reason we calculate a higher coverage than are the coverages expored to the histogram)

# run msprime (this produces the genomes)
python code/makeGenomes.py $gs $ploy $theta $td

# make reads from genomes
python code/makeReads.py $td

# run kmc
echo "########################################"
echo "KMC DB"
# parameter -m is specifying RAM used in GB, you might want to change it to something higher if you operate on a big server, the same with -t (number of threads)
kmc -k21 -t3 -m4 -cs"$kmax" -fm $td/reads.fa $td/db21 $td
echo "DB done."
echo ""
echo "########################################"
echo "MAKING SPECTRUM"
# run kmc_tools to make spectrum
kmc_tools transform  $td/db21 histogram  G$gs'P'$ploy'T'$theta.hist -cx"$kmax"
echo "removing zero lines"
awk '{ if( $2 != 0 ){ print $0 } }' G$gs'P'$ploy'T'$theta.hist > G$gs'P'$ploy'T'$theta.hist.no0 && rm G$gs'P'$ploy'T'$theta.hist
echo "Done."
echo ""

echo "########################################"
echo "Removing temp files..."
rm -rf $td # comment out for debug!
echo "Histogram written to G"$gs"P"$ploy"T"$theta".hist.no0"
echo "PIPELINE DONE."
