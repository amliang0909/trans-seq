'''version 0.1.0 Freeze 12/05/14'''

import subprocess
import os
import shlex
from operator import *
try:
    import HTSeq
except:
    pass
import HTSeq
import itertools
import re

class mate:
    def __init__(self,r1,r2):
        self.r1=r1
        self.r2=r2
        self.r=itertools.izip(HTSeq.FastqReader(r1),HTSeq.FastqReader(r2))
    
    def __iter__(self):
        return self
    
    def next(self):
        return self.r.next()

class counting_strategies:
    @staticmethod
    def counter(seq,error_rate):
        start,error_count=0,0
        c=zip(seq[start:len(seq)-1])
        for i in range(len(c)):
            if error_count <=error_rate:
                if c[i][0] =='T':
                    start+=1
                else:
                    error_count+=1
            else:
                break
        return start

    @staticmethod
    def percent(seq, window, threshold):
        windows=[]
        for i in range(len(seq)):
            windows.append(list(seq[i:i+window]).count('T'))
        try:
            if min(windows) <= threshold:
                for i in range(len(windows)):
                    if windows[i] <=threshold:
                        return i
            else:
                return len(windows)
        except ValueError:
            return 0

    @staticmethod
    def strings(seq, length, distance,adapter):
        string='T'*length
        if seq.count(string)>1:
            t_list=[m.end() for m in re.finditer(('(?=%s)'%(string)),seq)]
            t_list.append(len(seq))
            delta_list=[t_list[i+1]-t_list[i] for i in range(len(t_list)-1)]
            if max(delta_list)<=distance:
                return delta_list[-1]+adapter.count('T')
            else:
                for i in range(len(delta_list)):
                    if delta_list[i]>=distance:
                        return t_list[i]+adapter.count('T')
        else:
            try:
                return seq.index(string)+adapter.count('T')
            except ValueError:
                return adapter.count('T')

def out(odir):
    if odir[-1] == '/':
        pass
    else:
        odir=odir+'/'
    return odir

def check_bowtie():
    try:
        proc=subprocess.Popen(['bowtie','--version'],stdout=subprocess.PIPE)
        val=proc.communicate()[0]
    except OSError:
        errmsg='Error: bowtie not found on this system'
        print errmsg
        sys.exit(1)
        
        
    
def analyze(r1,r2,adapter,alignlen,error,wsize,wthreshold,slen,sdist,output_dir):
    out_file=out(output_dir)+'r1_reads_for_alignment.fastq'
    writeout=open(out_file,'w')
    tr,tad=0,0
    for paired in mate(r1,r2):
        tr+=1
        if adapter in paired[1].seq:
            tad+=1
            ta=counting_strategies.counter(paired[1].seq[paired[1].seq.index(adapter)+len(adapter):],error)
            tb=counting_strategies.percent(paired[1].seq[paired[1].seq.index(adapter)+len(adapter):],wsize,wthreshold)
            tc=counting_strategies.strings(paired[1].seq[paired[1].seq.index(adapter)+len(adapter):],slen,sdist,adapter)
            line='@%s_%s_%s_%s\n%s\n+\n%s\n' %(paired[1].name.split()[0],
                                               str(ta),
                                               str(tb),
                                               str(tc),
                                               paired[0].seq[4:4+alignlen],
                                               paired[0].qualstr[4:4+alignlen])
            writeout.writelines(line)
    writeout.close()
    return tr,tad,out_file

def tbin(read,start,end,bins):
    return int(round((float(read-start)/float(end-start))*bins))

def align(fastq,index,output_dir):
    sam_file=out(output_dir)+'aligned.sam'
    cmd='bowtie -m 1 %s %s -S %s' %(index,fastq,sam_file)
    cmd=shlex.split(cmd)
    line_logs='#bowtie cmd:\t%s\n\n' %(cmd)

    proc=subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    for line in proc.stderr:
        line_logs=line_logs+'%s' %(line)

    return line_logs,sam_file

def logs(output_dir,*args):
    log_file=out(output_dir)+'logs'
    log=open(log_file,'w')
    line='Total reads:\t%s\nReads with adapter:\t%s\n\n\n%s\n\n\nAssigned reads:%s' %(args[0],args[1],args[2],args[3])
    log.writelines(line)
    try:
        for line in args[4:]:
            log.writeline(line)
    except:
        pass
    log.close()
        
def samfile(sam):
    samreads={}
    with open(sam,'rU') as sam:
        for line in sam:
            if '@HD\tV' in line or '@SQ\t' in line or '@PG\t' in line:
                pass
            else:
                line=line.replace('\n','').split('\t')
                line[3]=int(line[3])
                line[4]=int(line[4])
                if line[1] == '*':
                    pass
                elif line[1] == '0':
                    line[1]='+'
                    try:
                        samreads[line[2]].append((line))
                    except KeyError:
                        samreads[line[2]]=[(line)]
                elif line[1] == '16':
                    line[1]= '-'
                    try:
                        samreads[line[2]].append((line))
                    except KeyError:
                        samreads[line[2]]=[(line)]

    for chrom in samreads.keys():
        samreads[chrom]=sorted(samreads[chrom],key=itemgetter(3,1))

    return samreads


def assign(transcripts,isoforms,samreads,data):
    assigned_reads=0
    chroms=list(set(transcripts.keys()) & set(samreads.keys()))
    for chrom in chroms:
        ts,reads=transcripts[chrom],samreads[chrom]
        counter=0
        for alignment in range(len(reads)):
            while counter < len(ts):
                if reads[alignment][3] >= ts[counter][2] and reads[alignment][3] <= ts[counter][3]:
                    if reads[alignment][1] == ts[counter][0]:
                        alpha=reads[alignment][0].split('_')
                        assigned_reads+=1
                        data[ts[counter][1]][0].append(tbin(reads[alignment][3],ts[counter][2],ts[counter][3],20))
                        data[ts[counter][1]][1].append(int(alpha[1]))
                        data[ts[counter][1]][2].append(int(alpha[2]))
                        data[ts[counter][1]][3].append(int(alpha[3]))
                        break
                    else:
                        break
                elif reads[alignment][3] >= ts[counter][3]:
                    counter+=1
                elif reads[alignment][3] <= ts[counter][2]:
                    break
    return assigned_reads,data

def gff(gff):
    gff_transcripts,data,isoforms,transcripts={},{},{},{}
    gff_file = HTSeq.GFF_Reader(gff, end_included=True )

    for feature in gff_file:
       if feature.type == "exon":
          transcript_id = feature.attr['transcript_id']
          if transcript_id not in gff_transcripts:
             gff_transcripts[ transcript_id ] = list()
          gff_transcripts[ transcript_id ].append( feature )
          
    for transcript_id in sorted(gff_transcripts):
        start,end,exons_start,exons_end,spliced_length=[],[],[],[],0
        for exon in gff_transcripts[transcript_id]:
            strand,chrom=exon.iv.strand,exon.iv.chrom
            start.append(exon.iv.start)
            end.append(exon.iv.end)
            spliced_length+= abs(int(exon.iv.start) - int(exon.iv.end))
            exons_start.append(int(exon.iv.start))
            exons_end.append(int(exon.iv.end))
        try:
            transcripts[chrom].append((strand,transcript_id,min(start),max(end),spliced_length,exons_start,exons_end))
            data[transcript_id]=[[],[],[],[],[],[]]
        except KeyError:
            transcripts[chrom]=[(strand,transcript_id,min(start),max(end),spliced_length,exons_start,exons_end)]
            data[transcript_id]=[[],[],[],[],[],[]]
        
        if '-' in transcript_id:
            b=transcript_id.split('-')
        elif '.' in transcript_id:
            b=transcript_id.split('.')
        else:
            b=transcript_id

        try:
            isoforms[b[0]].append((strand,transcript_id,min(start),max(end),spliced_length,exons_start,exons_end))
        except KeyError:
            isoforms[b[0]]=[(strand,transcript_id,min(start),max(end),spliced_length,exons_start,exons_end)]

    for chrom in transcripts.keys():
        transcripts[chrom]=sorted(transcripts[chrom],key=itemgetter(2))
    
    return transcripts,isoforms,data

def table(data,output_dir):
    out_file=out(output_dir)+'aligned_data_table'
    line='>Transcript_Name\tNo.Aligned_Reads\tBins\tA_length_counts\tA_length_window\tA_length_strings\n'
    outfile=open(out_file,'w')
    outfile.writelines(line)
    for k,v in data.iteritems():
        if len(v[0]) >=1:
            line='%s\t%s\t%s\t%s\t%s\t%s\n' %(k,
                                              len(v[0]),
                                              str(v[0]).strip('[]'),
                                              str(v[1]).strip('[]'),
                                              str(v[2]).strip('[]'),
                                              str(v[3]).strip('[]'))
            outfile.writelines(line)
    outfile.close()

def run(parser,options):
    os.mkdir(options.outputfolder)
    tr,tad,fastq=analyze(options.read1,
                         options.read2,
                         options.adapter,
                         options.r1len,
                         options.error,
                         options.window,
                         options.threshold,
                         options.string,
                         options.distance,
                         options.outputfolder)
    line_logs,sam=align(fastq,
                            options.index,
                            options.outputfolder)
    transcripts,isoforms,data=gff(options.GFF)
    samreads=samfile(sam)
    logs_assigned,assigned_data=assign(transcripts,isoforms,samreads,data)
    logs(options.outputfolder,tr,tad,line_logs,logs_assigned)
    table(assigned_data,options.outputfolder)

    print 'Alignments Completed, please check output folder %s'  %(options.outputfolder)
