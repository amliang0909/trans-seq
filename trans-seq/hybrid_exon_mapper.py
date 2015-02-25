'''Mohan T. Bolisetty hybrid_exon_mapper.py version 0.2 07/25/2013

TO FIX: 05/18/2013
1.BOWTIE EXON INDEX HEADERS HAVE COORDINATES WHICH REQUIRES STRING SPLITTING
2.REPETITIVE FILE READ AND ALIGNMENTS CAN BE REPLACED WITH APPROPRIATE FUNCTIONS
3.ITIERATE OVER FILES IN A FOLDER INSTEAD OF REPLACING FILE ASSIGNMENTS - should use a wrapper
4.LOLA/MOD REQURIES MULTIPLE LINE REPLACEMENTS

-----------------------------------------------------------------------------------------------------------

CHANGES 07/2013 ver0.2

1.ADDED sys.argv TO TAKE INPUTS FOR R1,R2 AND OUTPUT FROM COMMAND LINE
2.ADDED A SETTINGS .csv FILE WHICH CONTAINS path/to/index.ebwt, PRIMER SEQUENCES, ALIGNMENT LENGTHS FOR CE and VE
3.ADDED AN EXTRA OUTPUT COLUMN TO REPORT AMBIGIOUS ASSIGNMENTS
4.FIXED NO.4 05/18/2013 - via SETTINGS FILE

USAGE: python hybrid_exon_mapper.py <path/to/settings_file> <path/to/samplesheet> <input_folder_path> <output_folder_with_path>
------------------------------------------------------------------------------------------------------------
'''

''' THE FOLLOWING CODE READS MISEQ-R1,R2.

MISEQ-R1 IS ASSIGNED TO A dmel or dsim COMMON EXON USING LOCAL ALIGNMENT

MISEQ-R2 IS ASSIGNED TO A dmel or dsim VARIABLE EXON USING BOWTIE - REQUIRES
BOWTIE INDEXES .ebwt

THE ASSIGNMENTS ARE THEN PARSED TO ASSIGN VALUES M/M, M/S, S/M, S/S TO DETERMINE
HYBRID mRNAs FROM CIS mRNAs

-------------------------------------------------------------------------------------------------------------
'''

import sys
import os
import math
import fileinput
##import editdist
import shutil
import time
import csv

args_in = list(sys.argv)
settings_file=args_in[1]
samplesheet=args_in[2]
input_folder=args_in[3]
output_folder=args_in[4]

if os.path.exists(output_folder) == 'TRUE':
    print 'Output folder exists, please use a different folder name'
    sys.exit(0)
else:
    os.mkdir(output_folder, 0755)

''' READ LOLA/MOD SETTINGS FILE, .CSV FORMAT; CONTAINS INDEX LOCATIONS, PRIMER SEQUENCES, ALIGNMENT LENGTHS'''

settings_in = []
for line in fileinput.input([settings_file]):
    line=line.replace('\n','')
    settings_in.append(line.split(','))

dmel_index = settings_in[4][1]
dsim_index = settings_in[5][1]
mel_CE=settings_in[0][1]
sim_CE=settings_in[1][1]
primer=settings_in[2][1]
rev_primer=settings_in[3][1]
coordinates_file=settings_in[6][1]
max_CE = int(settings_in[7][1])
max_VE = int(settings_in[8][1])

''' READ SAMPLESHEET'''
samples_in = []
for line in fileinput.input([samplesheet]):
    line=line.replace('\n','')
    samples_in.append(line.split(','))

def local_alignment(primer,read):
    align_list = [[],[]]
    for i in range(len(read)-len(primer)):
        align_list[0].append(i)
        align_list[1].append(sum(map(lambda (x, y): 0 if x == y else 1, zip(primer, read[i:i+len(primer)]))))
    return align_list[0][align_list[1].index(min(align_list[1]))], min(align_list[1])

def read_out(out):
    a = out.readlines()
    for line in a:
        return line

''' READ IN EXON ASSIGNMENTS FROM .bed FILE, EXON ASSIGNMENTS IN coordinates[0] R1/R2 [1]=Mel/Mel [2]=Mel/Sim [3]= Sim/Mel [4]=Sim/Sim'''    

''' PRIMER SEQUENCES, EXPECTED COMMON EXON SEQUENCES AND PATH TO INDEXES'''

'''MISEQ-R1 INPUT AND COMMON EXON ASSIGNMENTS'''

''' DICTIONARY TO HOLD COMMON EXON ASSIGNMENTS'''
for i in range(len(samples_in)):
    R1=input_folder+samples_in[i][1]+'_'+samples_in[i][0]+'_L001'+'_R1'+'_001.fastq'
    R2=input_folder+samples_in[i][1]+'_'+samples_in[i][0]+'_L001'+'_R2'+'_001.fastq'
    print R1, R2
    output=output_folder+'/'+samples_in[i][8]

    coordinates=[[],[],[],[],[],[]]
    for line in fileinput.input([coordinates_file]):
        a = line.split('\t')
        coordinates[0].append(a[3])
        coordinates[1].append(0)
        coordinates[2].append(0)
        coordinates[3].append(0)
        coordinates[4].append(0)
        coordinates[5].append(0)

    mel_counter=0
    sim_counter=0
    mis_counter=0


    CE_cluster_id={}  
    lc_r1=0
    for line in fileinput.input([R1]):
        if lc_r1==0:
            cluster_idr1 = line.replace('\n','')
            cluster_idr1 = cluster_idr1[0:(cluster_idr1.index(' 1:N'))]
            lc_r1+=1
        elif lc_r1==1:
            seq = line.replace('\n','')
            lc_r1+=1
        elif lc_r1==2:
            symb=line.replace('\n','')
            lc_r1+=1
        elif lc_r1==3:
            qual=line.replace('\n','')
            lc_r1=0
            if primer in seq:
                trim5 = seq[(seq.index(primer)+len(primer)):(len(seq)-1)]
                if trim5[0:max_CE] in mel_CE[0:max_CE]:
                    CE_cluster_id[cluster_idr1]='dmel'
                    mel_counter=mel_counter+1
                elif trim5[0:max_CE] in sim_CE[0:max_CE]:
                    CE_cluster_id[cluster_idr1]='dsim'
                    sim_counter=sim_counter+1
                else:
                    CE_cluster_id[cluster_idr1]='misprimed'
                    mis_counter=mis_counter+1
            else:
                pass
    print mel_counter,'\t',sim_counter,'\t',mis_counter
    fileinput.close()

    '''MISEQ-R2 INPUT AND tmp/tmp.fastq READS OUTPUT FOR BOWTIE ALIGNMENT'''

    lc_r2=0
    os.mkdir('tmp', 0755)
    tmp_file = open('tmp/tmp_reads.fq','w')
    for line in fileinput.input([R2]):
        if lc_r2==0:
            cluster_idr2 = line.replace('\n','')
            cluster_idr2 = cluster_idr2[0:(cluster_idr2.index(' 2:N'))]        
            lc_r2+=1
        elif lc_r2==1:
            seq = line.replace('\n','')
            lc_r2+=1
        elif lc_r2==2:
            symb=line.replace('\n','')
            lc_r2+=1
        elif lc_r2==3:
            qual=line.replace('\n','')
            lc_r2=0
            if rev_primer in seq:
                trim3=seq[(seq.index(rev_primer)+len(rev_primer)):(len(seq)-1)]
                line1 = cluster_idr2
                line2 = trim3[4:max_VE]
                line3 = symb
                line4 = qual[(seq.index(rev_primer)+len(rev_primer)):(len(seq)-1)]
                line='%s\n%s\n%s\n%s\n' %(line1,line2,line3,line4[4:max_VE])
                tmp_file.writelines(line)
    tmp_file.close()
    fileinput.close()

    '''tmp/tmp.fastq ALIGNED TO dmel and dsim EXON INDEXES WITH DEFAULT BOWTIE SETTINGS'''

    cmd = 'bowtie --quiet %s %s %s --un %s' %(dmel_index,'tmp/tmp_reads.fq','tmp/tmp_align_dmel','tmp/tmp_unalign_dmel')
    process=os.popen(cmd)
    process.close()
    cmd = 'bowtie --quiet %s %s %s --un %s' %(dsim_index,'tmp/tmp_reads.fq','tmp/tmp_align_dsim','tmp/tmp_unalign_dsim')
    process=os.popen(cmd)
    process.close()

    '''DICTIONARY TO HOLD BOWTIE dmel: EXON ALIGNMENTS'''
    dmel_align_exon={} 
    '''DICTIONARY TO HOLD BOWTIE dmel: MISMATCH VALUES'''
    dmel_align_mismatch={}

    for line in fileinput.input(['tmp/tmp_align_dmel']):
        a=line.replace('\n','')
        b=a.split('\t')
        c=b[2].split(':')
        dmel_align_exon[('@'+b[0])]=c[0]
        try:
            dmel_align_mismatch[('@'+b[0])]=b[6]
        except IndexError:
            dmel_align_mismatch[('@'+b[0])]=''
    fileinput.close()        

    '''DICTIONARY TO HOLD BOWTIE dsim: EXON ALIGNMENTS'''
    dsim_align_exon={}
    '''DICTIONARY TO HOLD BOWTIE dsim: MISMATCH VALUES'''
    dsim_align_mismatch={} 

    for line in fileinput.input(['tmp/tmp_align_dsim']):
        a=line.replace('\n','')
        b=a.split('\t')
        c=b[2].split(':')
        dsim_align_exon[('@'+b[0])]=c[0]
        try:
            dsim_align_mismatch[('@'+b[0])]=b[6]
        except IndexError:
            dsim_align_mismatch[b[0]]=''
    fileinput.close()

    ''' DELETE tmp/ FOLDER. COMMENT OUR FOR TROUBLESHOOTING WILL CALL OS.ERROR IF tmp/ FOLDER EXISTS'''
    shutil.rmtree('tmp/') 

    ''' FINAL ASSIGNMENTS R1/R2 APPENDED TO COORDINATES'''
    same_counter=0
    for k in CE_cluster_id:
        if CE_cluster_id[k]=='dmel':
            try:
                exon_mel=dmel_align_exon[k]
            except KeyError:
                exon_mel='none'
            try:
                mismatch_mel=dmel_align_mismatch[k]
            except KeyError:
                mismatch_mel='none'
            try:
                exon_sim=dsim_align_exon[k]
            except KeyError:
                exon_sim='none'
            try:
                mismatch_sim=dsim_align_mismatch[k]
            except KeyError:
                mismatch_sim='none'

            if exon_mel == exon_sim and exon_mel!='none' and exon_sim!='none':
                if len(mismatch_mel) < len(mismatch_sim):
                    coordinates[1][coordinates[0].index(exon_mel)] = coordinates[1][coordinates[0].index(exon_mel)]+1
                elif len(mismatch_mel) > len(mismatch_sim):
                    coordinates[2][coordinates[0].index(exon_mel)] = coordinates[2][coordinates[0].index(exon_mel)]+1
                elif len(mismatch_mel) == len(mismatch_sim):
                    same_counter+=1
                    coordinates[5][coordinates[0].index(exon_mel)] = coordinates[5][coordinates[0].index(exon_mel)]+1
            elif exon_mel != exon_sim and exon_mel!='none':
                coordinates[1][coordinates[0].index(exon_mel)] = coordinates[1][coordinates[0].index(exon_mel)]+1
            elif exon_sim!=exon_mel and exon_sim!='none':
                coordinates[2][coordinates[0].index(exon_sim)] = coordinates[2][coordinates[0].index(exon_sim)]+1
        elif CE_cluster_id[k]=='dsim':
            try:
                exon_mel=dmel_align_exon[k]
            except KeyError:
                exon_mel='none'
            try:
                mismatch_mel=dmel_align_mismatch[k]
            except KeyError:
                mismatch_mel='none'
            try:
                exon_sim=dsim_align_exon[k]
            except KeyError:
                exon_sim='none'
            try:
                mismatch_sim=dsim_align_mismatch[k]
            except KeyError:
                mismatch_sim='none'
            
            if exon_mel == exon_sim and exon_mel!='none' and exon_sim!='none':
                if len(mismatch_mel) < len(mismatch_sim):
                    coordinates[3][coordinates[0].index(exon_mel)] = coordinates[3][coordinates[0].index(exon_mel)]+1
                elif len(mismatch_mel) > len(mismatch_sim):
                    coordinates[4][coordinates[0].index(exon_mel)] = coordinates[4][coordinates[0].index(exon_mel)]+1
                elif len(mismatch_mel) == len(mismatch_sim):
                    same_counter+=1
                    coordinates[5][coordinates[0].index(exon_mel)] = coordinates[5][coordinates[0].index(exon_mel)]+1                
            elif exon_mel != exon_sim and exon_mel!='none':
                coordinates[3][coordinates[0].index(exon_mel)] = coordinates[3][coordinates[0].index(exon_mel)]+1
            elif exon_sim!=exon_mel and exon_sim!='none':
                coordinates[4][coordinates[0].index(exon_sim)] = coordinates[4][coordinates[0].index(exon_sim)]+1


    ##'''OUTPUT COORDINATES LIST'''
    ##
    ##for i in range(len(coordinates[0])):
    ##	print coordinates[0][i], '\t', coordinates[1][i], '\t', coordinates[2][i], '\t', coordinates[3][i], '\t', coordinates[4][i]

    '''OUTPUT COORDINATES LIST TO FILE'''
    output_file = open(output,'w')
    line= 'Exon\tM/M\tM/S\tS\M\tS\S\tAmbigious\n'
    output_file.writelines(line)
    for i in range(len(coordinates[0])):
            line='%s\t%s\t%s\t%s\t%s\t%s\n' %(coordinates[0][i], coordinates[1][i], coordinates[2][i], coordinates[3][i], coordinates[4][i], coordinates[5][i])
            output_file.writelines(line)
    output_file.close()
