#**endseq**

endseq has three (currently only 2) modules. 
  
    usage: endseq [-h] [-v] {analysis,plots} ...
    
    optional arguments:
    -h, --help          show this help message and exit
    -v, --version     Installed endseq version
  
  [sub-commands]:
    {analysis,plots}
      analysis        endseq Alignments and A-Tail counting
      plots          	  Plots for endseq analysis


###**Analysis**

Analysis is the primary module that processes:

    reads in paired-end fastq files, 
    checks reads for adapter, 
    counts tail lengths, 
    aligns to reference genome with bowtie (outputs .sam file in output directory),
    assigns alignments to reference transcriptome,
    and outputs data_table
    that is necessary for plots and stats modules


    usage:
    endseq analysis [-h] [-q] -1 READ1 -2 READ2 -i INDEX -G GFF -o
                         OUTPUTFOLDER [-a ADAPTER] [-s R1LEN] [-e ERROR]
                         [-t THRESHOLD] [-w WINDOW] [-S STRING] [-d DISTANCE]

Required arguments:  

      -1			Fastq file read 1  
      -2      Fastq file read 2   
      -i			Bowtie index of reference genome (*.ebwt)  
      -G			GFF of reference transcriptome  
      -o			Output folder name  

Additional arguments:

      -s			Read 1 length used for alignment to genome, default=25
      -a			3’end adapter sequence default=GATGGTGCCTACAGTTTTT
      -e			Error rate - counting strategy 1, default=10
      -w			Window size - counting strategy 2, default=10
      -t			A tail threshold - counting strategy 2, default=4
      -s			String size - counting strategy 3, default =’TTTTT’
      -d			Max distance between strings - counting strategy 3, default =20

###**Plots**

      usage:
      endseq plots [-h] [-q] -t [TABLES [TABLE1,TABLE2,...,TABLEN]] -l
                          [LABELS [LABEL1,LABEL2,...,LABELN]]
                          [-p {KDE,CDF,HIST,Heatmap,Genewise,ScatterMatrix,SingleGene}]
                          [-m {counts,window,strings}] [-a MAX_LENGTH] [-s BINSIZE]
                          [-r COUNTS] [-c CONTROL] [-i GENEID]

Required arguments:

      -t /--tables 		Comma separated list of tables Table1,Table2,...,TableN
      -l/--label		Comma separated list of labels Label1,Label2,...,LabelN
      -p/--plottype		Choice: 
                            KDE,
                            CDF,
                            HIST,
                            Heatmap,
                            Genewise,
                            ScatterMatrix,
                            SingleGene

Optional arguments:

      -m/--metric		Choice:
                          counts (strategy1)
                          window (strategy2)
                          default=strings (strategy3) 

      -a/--max_length	   Max length of A-tail, default=200 
      -s/--binsize  		 Binsize, default=1 - change to 5 for Heatmap
      -r/--readsthreshold	Minimum read counts, default=30
      -c/--controlsample	Control Sample label, default=first label from labels
      -i/--geneid		Gene ID for Single Gene plots

###**Stats**
