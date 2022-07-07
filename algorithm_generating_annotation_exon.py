#!/usr/bin/env python
'''
Algorithm for generating annotation information structure for each exon (Python
script to parse exons from ChrAll_knownGene.txt file)

USAGE:
    python Algorithm_generating_annotation_exon.py ChrAll_knownGene.txt
'''

import sys, os
import pandas as pd, numpy as np

def args_check(numArgs):
    '''Checks if in proper number of arguments are passed; gives instructions on
    proper use'''
    if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
        print("\nPlease supply ONLY ONE file to parse.\n\nCorrect Syntax is:\n    python " + sys.argv[0] + " ChrAll_knownGene.txt\n\nThe output file will be the input file name with \".exons\" appended.\n")
        exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

# Uses the function defined above to check if the number of arguments is correct
args_check(2)

#number of bases on exon boundary
EXONBOUNDARY_OFFSET=1
UP_FLANK=2000
DOWN_FLANK=500

inFile = sys.argv[1]  # Stores file one for input checking; note that sys.argv[0] is this script file name
inFilet=os.path.splitext(inFile)[1]
try:
    print("Reading your data file "+inFile+" ...")
    if inFilet=='.xls' or inFilet=='.xlsx':
        a = pd.read_excel(inFile)
    else:
        a = pd.read_csv(inFile, delimiter="\t")
except IOError:  #handling Exceptions
    print("Failed to open " + inFile)
    exit(1)

a.insert(0, 'bin', 'NOINFO')
a.rename(columns={'#name': 'name','exonStarts':'exStart','exonEnds':'exEnd'},
         inplace=True)
a['exStart'] = a['exStart'].str.rstrip(',')
a['exEnd'] = a['exEnd'].str.rstrip(',')
b=a.set_index(['bin','name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonCount','proteinID','alignID']).apply(lambda x: x.str.split(',').explode())
b['ExonNumber'] = b.groupby(b.index).cumcount()+1
b['ExonNumber-'] = - b.groupby(b.index).cumcount()
c=b.reset_index()
c.loc[c.strand == '-','ExonNumber'] = c.exonCount + c['ExonNumber-']
c.drop(columns=['proteinID', 'alignID', 'ExonNumber-'],inplace=True)
c=c.astype({'exStart': 'int32','exEnd': 'int32'})
c['intronStart']=c['exEnd'].shift(1) + EXONBOUNDARY_OFFSET
c['intronEnd']=c['exStart'] - EXONBOUNDARY_OFFSET
c.loc[c.txStart == c.exStart,["intronStart","intronEnd"]] = np.nan

c['utr5Start']=np.nan
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart),'utr5Start']=c.exStart
c.loc[(c.strand=='-')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart>c.cdsEnd),
      'utr5Start']=c.exStart
c.loc[(c.strand=='-')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart<=c.cdsEnd),
      'utr5Start']=c.cdsEnd + 1
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd),
      'utr5Start']=c.cdsEnd + 1

c['utr5End']=np.nan
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd<c.cdsStart),
      'utr5End']=c.exEnd
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart),
      'utr5End']=c.cdsStart - 1
c.loc[(c.strand=='-')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd),
      'utr5End']=c.exEnd
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd),
      'utr5End']=c.exEnd

c['utr3Start']=np.nan
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart),'utr3Start']=c.exStart
c.loc[(c.strand=='+')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart>c.cdsEnd),
      'utr3Start']=c.exStart
c.loc[(c.strand=='+')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart<=c.cdsEnd),
      'utr3Start']=c.cdsEnd + 1
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd),
      'utr3Start']=c.cdsEnd + 1

c['utr3End']=np.nan
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd<c.cdsStart),
      'utr3End']=c.exEnd
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart),
      'utr3End']=c.cdsStart - 1
c.loc[(c.strand=='+')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd),
      'utr3End']=c.exEnd
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd),
      'utr3End']=c.exEnd

c['upstreamStart']=np.nan
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exStart==c.txStart),
      'upstreamStart']=c.txStart - UP_FLANK
c.loc[(c.strand=='-')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exEnd==c.txEnd),
      'upstreamStart']=c.txEnd + 1
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd)&(c.exEnd==c.txEnd),
      'upstreamStart']=c.txEnd + 1

c['upstreamEnd']=np.nan
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exStart==c.txStart),
      'upstreamEnd']=c.txStart - 1
c.loc[(c.strand=='-')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exEnd==c.txEnd),
      'upstreamEnd']=c.txEnd + UP_FLANK
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd)&(c.exEnd==c.txEnd),
      'upstreamEnd']=c.txEnd + UP_FLANK

c['downstreamStart']=np.nan
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exStart==c.txStart),
      'downstreamStart']=c.txStart - DOWN_FLANK
c.loc[(c.strand=='+')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exEnd==c.txEnd),
      'downstreamStart']=c.txEnd + 1
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd)&(c.exEnd==c.txEnd),
      'downstreamStart']=c.txEnd + 1

c['downstreamEnd']=np.nan
c.loc[(c.strand=='-')&(c.exStart<=c.cdsStart)&(c.exStart==c.txStart),
      'downstreamEnd']=c.txStart - 1
c.loc[(c.strand=='+')&(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exEnd==c.txEnd),
      'downstreamEnd']=c.txEnd + DOWN_FLANK
c.loc[(c.strand=='+')&(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd)&(c.exEnd==c.txEnd),
      'downstreamEnd']=c.txEnd + DOWN_FLANK

c['cdsStart_']=np.nan
c.loc[(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart<=c.cdsEnd),
      'cdsStart_']=c.exStart
c.loc[(c.exStart>c.cdsStart)&(c.exEnd<=c.cdsEnd),'cdsStart_']=c.exStart
c.loc[(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart),'cdsStart_']=c.cdsStart

c['cdsEnd_']=np.nan
c.loc[(c.exStart>c.cdsStart)&(c.exEnd>c.cdsEnd)&(c.exStart<=c.cdsEnd),
      'cdsEnd_']=c.cdsEnd
c.loc[(c.exStart>c.cdsStart)&(c.exEnd<=c.cdsEnd),'cdsEnd_']=c.exEnd
c.loc[(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd<c.cdsEnd),
      'cdsEnd_']=c.exEnd
c.loc[(c.exStart<=c.cdsStart)&(c.exEnd>=c.cdsStart)&(c.exEnd>=c.cdsEnd),
      'cdsEnd_']=c.cdsEnd

c.drop(columns=['cdsStart', 'cdsEnd', 'strand', 'exonCount'],inplace=True)
c.rename(columns={'cdsStart_': 'cdsStart','cdsEnd_':'cdsEnd'},
         inplace=True)
c=c[['bin','name','chrom','txStart','txEnd','cdsStart','cdsEnd','ExonNumber',
     'exStart','exEnd','intronStart','intronEnd','utr5Start','utr5End','utr3Start',
     'utr3End','upstreamStart','upstreamEnd','downstreamStart','downstreamEnd']]
c=c.astype({'cdsStart': 'Int32','cdsEnd': 'Int32','intronStart': 'Int32',
            'intronEnd': 'Int32','utr5Start': 'Int32','utr5End': 'Int32',
            'utr3Start': 'Int32','utr3End': 'Int32','upstreamStart': 'Int32',
            'upstreamEnd': 'Int32','downstreamStart': 'Int32','downstreamEnd':'Int32',})

outFile = inFile + ".exons"
try:
    print("Saving data to file "+outFile+" ...")
    c.to_csv(outFile,index=False,sep='\t',na_rep='NA')
except IOError:
    print("Failed to create " + outFile)
    exit(1)
