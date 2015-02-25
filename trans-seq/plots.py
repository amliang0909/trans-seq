''' FREEZE version 0.1.0 - test 12.04.14'''

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import scipy
import scipy.cluster.hierarchy as sch
from operator import *
pd.core.config.option_context('mode.use_inf_as_null',True)

try:
    import seaborn as sns  
    rcdefsns = plot.rcParams.copy()
except:
    pd.options.display.mpl_style = 'default'
    plot.ioff()

def binning(maxlength,size):    
    bins=np.arange(5,maxlength,size)
    return bins

def mhist(values,*args):
    return np.histogram(values,args[0],normed=True)[0]

def kde(values,*args):
    try:
        density=scipy.stats.gaussian_kde(values)
        ys=density(args[0])
        return ys
    except (ValueError,np.linalg.LinAlgError):
        return np.nan

def cdf(values,*args):
    try:
        a=scipy.stats.cumfreq(values,len(args[0]))
        return a[0]/max(a[0])
    except (ValueError,np.linalg.LinAlgError):
        return np.nan

def cumulative_plots(tables,plottype,metric):
    fig = plot.figure(frameon=False,facecolor='white')
    color=colors(tables)
    for i in range(len(tables)):
        plot.plot(tables[i][1][plottype].sum()/len(tables[i][1][plottype]),color=color[i],label=str(tables[i][0]))
    plot.xlabel('Tail Length (%s)' %(metric))
    plot.ylabel('Frequency')
    plot.legend()
    return fig
    
def table(filein,metric,counts,maxlength,size):
    table = pd.io.parsers.read_table(filein,index_col=0)
    metric='A_length_'+str(metric)
    table[metric]=table[metric].str.split(',').apply(np.asarray,dtype=int)
    
    index=table['No.Aligned_Reads'].where(table['No.Aligned_Reads']>=counts).dropna().index
    table=table.ix[index]
    table['MEDIAN']=table[metric].apply(np.median)
    table['KDE']=table[metric].apply(kde,args=(binning(maxlength,size),))
    table['HIST']=table[metric].apply(mhist,args=(binning(maxlength,size),))
    table['CDF']=table[metric].apply(cdf,args=(binning(maxlength,size),))    
    return table
  
def colors(tables):
    try:
        import brewer2mpl
        allcolors1=brewer2mpl.get_map('Set1', 'Qualitative',7).mpl_colors
        allcolors2=brewer2mpl.get_map('Set2', 'Qualitative',6).mpl_colors
        allcolors=allcolors1+allcolors2
    except:
        allcolors=['b','g','r','c','m','y','k','w']

    while len(tables)> len(allcolors):
        allcolors+=allcolors

    return allcolors[:len(tables)]

def common_index(tables):
    index=tables[0][1].index
    for i in range(1,len(tables)):
        index=index.intersection(tables[i][1].index)
    return index

def genewise(tables,fig):
    color=colors(tables)
    index=common_index(tables)
    for i in index:
        fig, ((ax1,ax2),(ax3,ax4)) = plot.subplots(nrows=2, ncols=2)
        fig.suptitle(str(i))
        reads=[[],[]]
        for q in range(len(tables)):
            ax1.plot(tables[q][1].ix[i].KDE,color=color[q],label=tables[q][0])
            ax2.plot(tables[q][1].ix[i].CDF,color=color[q],label=tables[q][0])
            reads[0].append(tables[q][0])
            reads[1].append(tables[q][1].ix[i]['No.Aligned_Reads'])
        ax1.legend()
        ax1.set_ylabel('Frequency')
        ax1.set_xlabel('Tails')
        ax2.legend()
        ax2.set_ylabel('Frequency')
        ax2.set_xlabel('Tails')
        ind=np.arange(len(reads[0]))+0.5
        ax3.bar(ind,reads[1],color='k')
        ax3.set_xticks(ind+0.5)
        ax3.set_xticklabels(reads[0])
        ax3.set_ylabel('Reads')
        fig_pdf.savefig()
        plot.close(fig)

def singlegene(tables,gene):
    color=colors(tables)
    fig, ((ax1,ax2),(ax3,ax4)) = plot.subplots(nrows=2, ncols=2)
    fig.suptitle(str(gene))
    reads=[[],[]]
    for q in range(len(tables)):
        ax1.plot(tables[q][1].ix[gene].KDE,color=color[q],label=tables[q][0])
        ax2.plot(tables[q][1].ix[gene].CDF,color=color[q],label=tables[q][0])
        reads[0].append(tables[q][0])
        reads[1].append(tables[q][1].ix[gene]['No.Aligned_Reads'])
    ax1.legend()
##        ax1.set_xlim([0,80])
    ax1.set_ylabel('Frequency')
    ax1.set_xlabel('Tails')
    ax2.legend()
    ax2.set_ylabel('Frequency')
    ax2.set_xlabel('Tails')
    ind=np.arange(len(reads[0]))+0.5
    ax3.bar(ind,reads[1],color='k')
    ax3.set_xticks(ind+0.5)
    ax3.set_xticklabels(reads[0])
    ax3.set_ylabel('Reads')
    return fig        

def heatmap(tables):
    fig, axes = plot.subplots(nrows=1, ncols=len(tables))
    index=common_index(tables)
    if len(tables) > 1:
        for data, ax in zip(tables, axes.flat):
            datap=np.row_stack(data[1].ix[index]['HIST'].values)
            heatmap = ax.matshow(datap, aspect='auto', cmap=plot.cm.Greens,vmin=0,vmax=0.05)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(str(data[0]))
            cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
            fig.colorbar(heatmap, cax=cax)
    else:
        ax=axes
        data=tables[0]
        datap=np.row_stack(data[1].ix[index]['HIST'].values)
        heatmap = ax.matshow(datap, aspect='auto', cmap=plot.cm.Greens,vmin=0,vmax=0.05)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(str(data[0]))
        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        fig.colorbar(heatmap, cax=cax)
    return fig
        

def wrap(txt, width=8):
    '''helper function to wrap text for long labels'''
    import textwrap
    return '\n'.join(textwrap.wrap(txt, width))

def scattermatrix(tables):
    fig = plot.figure(frameon=False,facecolor='white')
    index=common_index(tables)
    data=pd.DataFrame(index=index)
    for i in tables:
        data[i[0]]=i[1].ix[index]['MEDIAN']
    axs=pd.scatter_matrix(data, alpha=0.2, figsize=(8,8), diagonal='none', marker='.',)
    
    for ax in axs[:,0]:
        ax.grid('off', axis='both')
        ax.set_ylabel(wrap(ax.get_ylabel()), rotation=0, va='center', labelpad=30)
        ax.set_yticks([])
    for ax in axs[-1,:]:
        ax.grid('off', axis='both')
        ax.set_xlabel(wrap(ax.get_xlabel()), rotation=90)
        ax.set_xticks([])
    return fig        


def ksstat(control,experimental):
    return scipy.stats.ks_2samp(control,experimental)

def ks(tables,metric,control):
    metric='A_length_'+str(metric)
    index=common_index(tables)
    data=pd.DataFrame(index=index)
    ks_data=pd.DataFrame(index=index)
    for i in tables:
        data[i[0]]=i[1].ix[index][[metric],['HIST']]

    for column in data.columns:
        ksstat='KS-Stat'+str(control)+str(column)
        pval='Pvalue'+str(control)+str(column)
        data[ksstat],data[pval]=zip(*data.apply(lambda x:ksstat(x[control],x[column]),axis=1))

    return data        
        
def heatmap_ks(tables,metric):
    data=ks(tables,metric)

def control_sample(labels,control):
    if control == '':
        return labels[0]
    else:
        if control in labels:
            return control
        else:
            print 'Control Label is not present in Labels revert to sample %s' %(options.label[0])
            return label[0]    

##def cluster(tables,sample):
##    'test'
##


def run(parser,options):
    data=[]

    for i in range(len(options.tables[0])):
        data.append((options.labels[0][i],
                       table(options.tables[0][i],
                             'strings',
                             options.counts,
                             options.max_length,
                             options.binsize)))

    fig_pdf=PdfPages('%s.pdf' %(options.plottype))
    
    if options.plottype=='KDE' or options.plottype=='CDF' or options.plottype=='HIST':
        fig=cumulative_plots(data,options.plottype,options.metric)
        fig_pdf.savefig()
        plot.close(fig)
        fig_pdf.close()        
    elif options.plottype=='Heatmap':
        fig=heatmap(data)
        fig_pdf.savefig()
        plot.close(fig)
        fig_pdf.close()        
    elif options.plottype=='Genewise':
        genewise(data)
        fig_pdf.close()
    elif options.plottype=='ScatterMatrix':        
        fig=scattermatrix(data)
        fig_pdf.savefig()
        plot.close(fig)        
        fig_pdf.close()
    elif options.plottype=='SingleGene':        
        fig=singlegene(data,options.geneid)
        fig_pdf.savefig()
        plot.close(fig)        
        fig_pdf.close()
        
