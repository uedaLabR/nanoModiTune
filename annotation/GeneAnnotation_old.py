import pandas as pd
from itertools import product
import numpy as np
seqprepast = 5


import pysam
class GenomeAnnotator:

    def __init__(self, gtf_file,stringtie_gtf,neibormargin=7):

        self.data = pd.read_csv(gtf_file, delimiter='\t')
        data = self.data
        self.neibormargin = neibormargin
        didx = 0
        self.ididx = {}
        for id in  data['#hg38.knownGene.name']:

            self.ididx[id] = didx
            didx+=1

        # for genes
        self.chromosome_data = {}
        self.max_tx_length = {}
        # from matrix which include stringtie2 result as well
        self.chromosome_data_st = {}
        self.neighbor_seq = None

        strands = ['+','-']
        for chrom,strand in product(data['hg38.knownGene.chrom'].unique(),strands):

            key = chrom+":"+str(strand)
            chrom_data = data[data['hg38.knownGene.chrom'] == chrom]
            chrom_data = chrom_data[chrom_data['hg38.knownGene.strand'] == strand]
            chrom_data = chrom_data.sort_values(by='hg38.knownGene.txStart').reset_index(drop=True)
            chrom_data['txLength'] = chrom_data['hg38.knownGene.txEnd'] - chrom_data['hg38.knownGene.txStart']
            max_tx_length = chrom_data['txLength'].max()
            if not np.isnan(max_tx_length):
                self.chromosome_data[key] = chrom_data
                self.max_tx_length[key] = max_tx_length


        #
        strands = [1, -1]
        self.exonsdata = pd.read_csv(convertMatrix, delimiter='\t', header=None)
        new_headers = ['ID', 'chrom', 'strand', 'exons']
        self.exonsdata.columns = new_headers
        self.exonsdata['start'] = self.exonsdata.iloc[:, -1].apply(lambda x: int(x.split(',')[0].split(':')[0]))

        for chrom, strand in product(self.exonsdata['chrom'].unique(), strands):

            key = chrom + ":" + str(strand)
            chrom_data_st = self.exonsdata[self.exonsdata['chrom'] == chrom]
            chrom_data_st = chrom_data_st[chrom_data_st['strand'] == strand]
            chrom_data_st = chrom_data_st.sort_values(by='start').reset_index(drop=True)
            self.chromosome_data_st[key] = chrom_data_st

        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def reverse_complement(self,seq):

        reverse_complement_seq = "".join(self.complement.get(base, base) for base in reversed(seq))
        return reverse_complement_seq

    def getSeq(self, chrom, strand, neighborseqInterval,record):

        allseq = ""
        for s, e in neighborseqInterval:
            # seq = self.aligner.seq(chrom, s, e)
            # seq = self.fastafile.fetch(chrom,s,e).upper()
            seq = str(record.seq[s:e]).upper()
            if seq is None:
                continue
            allseq = allseq + seq

        if not strand:
            allseq = self.reverse_complement(allseq)
        return allseq


    def getSeqSimple(self, chrom, strand, s, e,record):

        # seq = self.aligner.seq(chrom, s, e)
        #seq = self.fastafile.fetch(chrom, s, e).upper()
        seq = str(record.seq[s:e]).upper()
        if seq is None:
            return ""
        if not strand:
            seq = self.reverse_complement(seq)
        return seq

    def anno(self,row,bstrand,position,margin, neighbor_seq,record):

        spliceJunction = False
        gene_id ,gene_symbol = row['#hg38.knownGene.name'],row['hg38.kgXref.geneSymbol']
        txStart, txEnd = row['hg38.knownGene.txStart'], row['hg38.knownGene.txEnd']
        chrom = row['hg38.knownGene.chrom']
        if position < txStart - margin:
            return None
        if position > txEnd + margin:
            return None


        streamlabel1 = 'upstream'
        utrlabel1 = '5pUTR'
        streamlabel2 = 'downstream'
        utrlabel2 = '3pUTR'
        if not bstrand:
            streamlabel1 = 'downstream'
            streamlabel2 = 'upstream'
            utrlabel1 = '3pUTR'
            utrlabel2 = '5pUTR'

        if position < txStart:
            dist = txStart - position
            return (streamlabel1, gene_id,gene_symbol,dist,neighbor_seq,spliceJunction)

        if position > txEnd:
            dist = position - txEnd
            return (streamlabel2, gene_id, gene_symbol,dist,neighbor_seq,spliceJunction)

        cdsStart, cdsEnd = row['hg38.knownGene.cdsStart'], row['hg38.knownGene.cdsEnd']
        cdsflg = "CDS"
        if cdsStart == cdsEnd:
            cdsflg = "noncoding"
            utrlabel1 =  "noncoding"
            utrlabel2 = "noncoding"


        if position < cdsStart:
            dist = abs(cdsStart - position)
            return (utrlabel1, gene_id, gene_symbol, dist,neighbor_seq,spliceJunction)

        if position > cdsEnd:
            dist = abs(position - cdsEnd)
            return (utrlabel2, gene_id, gene_symbol, dist,neighbor_seq,spliceJunction)

        exonStarts = list(map(int, row['hg38.knownGene.exonStarts'].strip(',').split(',')))
        exonEnds = list(map(int, row['hg38.knownGene.exonEnds'].strip(',').split(',')))
        dist = 0
        prevexonend = 0
        idx = -1
        for exonStart, exonEnd in zip(exonStarts, exonEnds):

            idx+=1
            if exonStart <= position <= exonEnd:

                dist += (position-exonStart)
                if bstrand:
                    dist = dist - (cdsStart-txStart)
                else:
                    dist = dist - (txEnd-cdsEnd)

                d1 = position-exonStart
                if d1 < self.neibormargin:
                    diff = self.neibormargin-d1
                    itvl = [(prevexonend-diff-1, prevexonend), (exonStart,  position+self.neibormargin-diff)]
                    neighbor_seq = self.getSeq(chrom,bstrand,itvl,record)
                    spliceJunction = True


                d2 = exonEnd-position
                if d2 < self.neibormargin:

                    diff = self.neibormargin - d2
                    nidx = idx+1
                    if nidx < len(exonStarts):
                        nextexonStart = exonStarts[nidx]
                        itvl = [(position-self.neibormargin-1, exonEnd), (nextexonStart,  nextexonStart+diff)]
                        neighbor_seq = self.getSeq(chrom,bstrand,itvl,record)
                        spliceJunction = True

                return (cdsflg, gene_id, gene_symbol, dist,neighbor_seq,spliceJunction)

            elif position < exonStart:

                if bstrand:
                    dist = dist - (cdsStart-txStart)
                else:
                    dist = dist - (txEnd-cdsEnd)

                return ("INTRON", gene_id, gene_symbol, dist,neighbor_seq,spliceJunction)

            else:

                dist += (exonEnd-exonStart)
            prevexonend = exonEnd

        return "Error"


    def select_tuple(self,tuples):

        priority_order = ['CDS', 'INTRON', '3pUTR', '5pUTR','noncoding','novel_ncRNA','genomic', 'upstream', 'downstream','genomic']
        for category in priority_order:
            filtered_tuples = [tup for tup in tuples if tup[0] == category]
            if filtered_tuples:
                return min(filtered_tuples, key=lambda tup: abs(tup[3]))

        return None

    def inExon(self,position,exons):

        intervals = exons.split(',')
        for interval in intervals:
            start, end = map(int, interval.split(':'))
            if start <= position <= end:
                return True

        return False

    def getNeighborSeq(self,chrom,bstrand,pos,record):

        neighbor_seq = self.getSeqSimple(chrom, bstrand, pos - self.neibormargin - 1, pos + self.neibormargin, record)
        return neighbor_seq

    def annotate_genomic_position(self, chrom,bstrand,position,geneID,record):

        margin = 1000
        annotations = []
        retanno = None

        neighbor_seq = self.getSeqSimple(chrom,bstrand,position-self.neibormargin-1,position+self.neibormargin,record)
        self.neighbor_seq = neighbor_seq

        if bstrand:
            strand = '+'
        else:
            strand = '-'

        annotation_set=False
        if len(geneID) > 10:
            if geneID in self.ididx:

                idx = self.ididx[geneID]
                row = self.data.iloc[idx]
                annotation = self.anno(row, bstrand, position, margin, neighbor_seq,record)
                if annotation is not None:
                    print('anno',annotation, annotation is None,len(annotation))
                    retanno = annotation
                    annotation_set = True
                    neighbor_seq = annotation[4]


        key = chrom + ":" + str(strand)
        if key not in self.chromosome_data:
            return neighbor_seq, False, None

        chrom_data = self.chromosome_data[key]
        max_tx_length = self.max_tx_length[key]

        idx = chrom_data['hg38.knownGene.txStart'].searchsorted(position-margin-max_tx_length, side='right') - 1
        idx2 = chrom_data['hg38.knownGene.txStart'].searchsorted(position+margin + max_tx_length, side='right')
        if idx < 0 and idx2 < 0:
            return neighbor_seq, False, None
        elif idx < 0 and idx2 >= 0:
            idx = 0

        intronExsist = False
        exonExist = False
        for n in range(idx,idx2):

            row = chrom_data.iloc[n]
            start = row['hg38.knownGene.txStart']
            if position+margin < start:
                break
            annotation = self.anno(row,bstrand,position,margin, neighbor_seq,record)
            if annotation is not None:
                if annotation[0] == "INTRON":
                    intronExsist = True
                if "UTR" in annotation[0] or annotation[0] == "CDS":
                    exonExist = True

                if annotation_set == False:
                    annotations.append(annotation)
                neighbor_seq = annotation[4]

        if retanno is None:
            retanno = self.select_tuple(annotations)


        possibleASRegion = intronExsist

        #check introm from matrix info
        if bstrand:
            strand = '1'
        else:
            strand = '-1'

        key = chrom + ":" + str(strand)
        chrom_data_st = self.chromosome_data_st[key]

        idx = chrom_data_st['start'].searchsorted(position - margin-max_tx_length, side='right') - 1
        idx2 = chrom_data_st['start'].searchsorted(position + margin + max_tx_length, side='right')
        if idx2 >= idx and idx >=0 and idx2 >=0:
        #
            intron = True
            for n in range(idx,idx2):

                row = chrom_data_st.iloc[n]
                start = row['start']
                exons = row['exons']
                if self.inExon(position,exons):
                    intron = False
                    break
                if position+margin < start:
                    break

            if intron:
                possibleASRegion = True

        return neighbor_seq,possibleASRegion,retanno






# file_path = '/mnt/share/ueda/RNA004/pipeline/gencode_v44.tsv'
# ref = "/share/reference/hg38.fa"
# convertMatrix = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/bam_pass/U87wt.matrix"
# annotator = GenomeAnnotator(file_path,ref,convertMatrix)
# chrom = 'chr1'
# position = 10399634
# annotations = annotator.annotate_genomic_position(chrom,True,position,'ENST00000477958.5')
# for annotation in annotations:
#     print(annotation)
#
# position = 9749873
# annotations = annotator.annotate_genomic_position(chrom,False,position,'ENST00000464286.1')
# for annotation in annotations:
#     print(annotation)
