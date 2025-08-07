This repository is intended to accompany the manuscript "Small phytoplankton community composition cycles annually with a coastal bloom" 

 
The analysis presented in that work is carried out in a variety of formats (R, bash, and matlab). Here we outline and document all the steps in one place. 


## Part 1: Going from raw NCBI sequences to ASVs 

	# MVCO_process_4batches.R
	    This script starts from fastq data downloaded from our NCBI projects and processes them with dada2. It compares the asvs to the PR2 database. 

	    It produces the following output files: 
            	Outputs4batches_all.Rdata, 
            	taxa4batch.txt
		taxa4_boot.txt
		taxa4tab.txt
		mergtab4_trans.txt

## Part 2: Sort and do some statistical analyses in MATLAB

  	# Metadata file is MVCO_metadata_AllSeqeunceSamples.csv

	# Contingencies for plotting in MatLab
		brewermap and viridis

	# Part2_sorting.m  
	     This script loads mergtab4_trans.txt and metadata and carries out initial statistical analysis.

	     The sections of the script are as follows: 
		# Load and organize initial data
		# Cut to Small Phytoplankton and save 'Taxa_Included_in_SmallPhytoplankton.csv'
		# Report numbers of ASVs and look at small phyto grouped at Division (Figure S1)
		# Look at chlorophyta grouped at Class (Figure S1 part 2) 
		# Calculate Shannon Diversity Indices at Genus level 
		# Supplementary Figure S2
		# Barplots for Figure 1
		# PCA at Species level assignments - Figure 2
		# Load MVCO Flow Cytometry concentration data 
		# Compare seasonality of concentration to seasonality of PCA - Figure S7
		# Look at annual cycles as climatologies - Figure 3
		# Supplementary Figure S4 - correlations between PC1 and env. 
		# Test temperature vs. day of year as predictors of PC1 
		# Generate time series of Phaeocystis data - Figure S6

## Part 3: SparCC 

	# Part3_CreateSparCCInputs.m 
	      This script simply takes outputs from Part 1 and converts them into csv that is useful for SparCC. There is certainly a way to do this without going through MATLAB, but I am more comfortable managing tables in MATLAB. 

	      Outputs are chrono_species_gt80_forR.csv and chrono_species_gt80_taxa.csv 

	# sparcc_analysis.R
	      This scipt uploads data from csv files made above, and metadata csv. Then runs 
		SparCC correlation analysis to produce an adjacency matrix which we cut down to 
		consider picophytoplankton chlorophytes only. We export the adjacency matrix as 
		two new csvs: 

		Outputs are sparCCtable.csv and sparCCtaxa.csv
			as well as sparccout.RData for those who prefer R objects. 

## Part 4: Topic Model Inference
	
   	# 4. 1: Format inputs for rost-cli 
	    # Part4p1_PrepWordlist.m 
	     	This script is used to convert our tables of counts into "wordlists" which are formatted as a row for each sample and repeated instances of integers for each sequence assigned to a taxon. 
	
	    	This script calls the function Barplot_to_Wordlist.m 

	    	Outputs: 4Batch_chlorophytes.csv and 4Batch_SmallPhytos.csv

	    # Next we need to remove the commas at the end of each row of the csv (because the rows in the matlab cell all had to be of the same length, but ROST will interpret blank ,,, strings as 0s). You can do this however you would like, but I did in a terminal window with the following commands. 

		sed ‘s/,,*,//’ <4Batch_chlorophytes.csv >4Batch_chlorophytes_nocommas.csv
		sed ‘s/,,*,//’ <4Batch_SmallPhytos.csv >4Batch_SmallPhytos_nocommas.csv

   	# 4.2 Execute Realtime Online Spatiotemporal Topic Modeling (ROST)
	    All scripts needed for ROST are available and documented at https://gitlab.com/warplab/rost-cli

	    Once packages are installed etc. We apply the method to our data using the following: 

	   ./topics.refine.t --in.words=../4Batch_chlorophytes_no_commas.csv --iter=100 --alpha=0.1 --beta=0.5 --vocabsize=69 -K 2

	    The vocabsize parameter is the number of species included in the wordlist. If you are applying this method to a different set of data, the vocabsize will be the height of the word_dict variable or 1 greater than the highest number in the wordlist table (excluding the timestamps in the first column) since the "words" begin with number 0.  

	    This command will create files 
   		perplexity.csv
		perplexity.iter.csv
		topicmodel.csv
		topics.csv
		topics.maxlikelihood.csv
		topics.position.csv

   	# 4.3 Return to MATLAB for analysis of the topic model outputs. 
	     # Part4p3_AnalyzeTopicResults.m 
	 	This script reads the results of the topic model assignments and generates barplot figures for prevalence of each topic within each sample and top taxa represented within each topic. 

	     The sections of the script are as follows: 
		# Load data and sort metadata chronologically 
		# Make time series of co-occurrence community proportions - Figure 5 
		# Make barplots of co-occurrence community species - Figure S4 




