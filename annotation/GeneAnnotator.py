import pandas as pd
from itertools import product
import numpy as np
seqprepast = 5
from pybedtools import BedTool

class IntervalList:

    def __init__(self,neibormargin):
        self.intervals = []
        self.genes = []
        self.merged = None
        self.neibormargin = neibormargin

    def addInterval(self,gene):


        s = gene.start - self.neibormargin
        if s < 1:
            s  =1
        e = gene.end + self.neibormargin
        if e < self.neibormargin:
            e  = self.neibormargin
        new_gene = (gene.chrom, s, e, gene.strand)
        self.intervals.append(new_gene)
        self.genes.append(gene)

    def resolveinterval(self):

        genes = BedTool(self.intervals)
        sorted = genes.sort()
        self.merged = sorted.merge(d=0)
        # print(self.merge)


class GeneHolder:

    def __init__(self,neibormargin=20):

        self.plusgenes = {}
        self.minusgenes = {}
        self.neibormargin = neibormargin

        self.plusMaingenes = {}
        self.minusMaingenes = {}

    def addGene(self,gene):

        strand = gene.strand
        chrom = gene.chrom
        if strand == '+':
            if chrom not in self.plusgenes:
                self.plusgenes[chrom] = IntervalList(self.neibormargin)
            self.plusgenes[chrom].addInterval(gene)
        else:
            if chrom not in self.minusgenes:
                self.minusgenes[chrom] = IntervalList(self.neibormargin)
            self.minusgenes[chrom].addInterval(gene)

    def addMainGene(self,gene):

        strand = gene.strand
        chrom = gene.chrom
        if strand == '+':
            if chrom not in self.plusMaingenes:
                self.plusMaingenes[chrom] = []
            self.plusMaingenes[chrom].append(gene)
        else:
            if chrom not in self.minusMaingenes:
                self.minusMaingenes[chrom] = []
            self.minusMaingenes[chrom].append(gene)

    def addGenes(self, genes):

        if genes is not None:
            for gene in genes:

                if len(gene)>2:
                    if (gene[2]=="gene"):
                        self.addMainGene(gene)

                self.addGene(gene)

    def resolveinterval(self):

        for key,value in self.plusgenes.items():

            intervallist = self.plusgenes[key]
            intervallist.resolveinterval()

        for key, value in self.minusgenes.items():
            intervallist = self.minusgenes[key]
            intervallist.resolveinterval()

    def getIntervalsByKey(self):

        retdict = {}
        for key,intervallist in self.plusgenes.items():

            retlist = []
            for iv in intervallist.merged:
                retlist.append(iv)

            retdict[key+"-"+"p"] = retlist

        for key, intervallist in self.minusgenes.items():

            retlist = []
            for iv in intervallist.merged:
                retlist.append(iv)
            retdict[key + "-" + "m"] = retlist

        return retdict



from pybedtools import BedTool
import pysam

def extract_gene_names(feature):
    attributes = feature.attrs
    if 'gene_name' in attributes:
        return attributes['gene_name']
    return None

def get_level(gene):
    level_str = gene.attrs.get('level')
    if level_str is not None:
        try:
            return int(level_str)
        except ValueError:
            return float('inf')
    else:
        return float('inf')

class GenomeAnnotator:


    def __init__(self, ref,gtf_file,stringtie_gtf,neibormargin=20):

            print(gtf_file)
            print(stringtie_gtf)

            self.genes = BedTool(gtf_file)
            if stringtie_gtf!=None and len(stringtie_gtf)>0:
                self.stringtie_genes = BedTool(stringtie_gtf)


            def is_gene_or_transcript(feature):
                return feature[2] == 'gene' or feature[2] == 'transcript'

            def is_UTR(feature):
                return "UTR" in feature[2]

            genes = self.genes.filter(is_gene_or_transcript)
            # print(genes)
            if stringtie_gtf:
                stgenes = self.stringtie_genes.filter(is_gene_or_transcript)

            utrs = self.genes.filter(is_UTR)

            self.gh = GeneHolder(neibormargin)
            self.gh.addGenes(genes)
            if stringtie_gtf:
                self.gh.addGenes(stgenes)
            self.gh.resolveinterval()

            self.gh_utr = GeneHolder(neibormargin)
            self.gh_utr.addGenes(utrs)
            self.gh_utr.resolveinterval()
            print("finish")


    def getIntervalsByKey(self):

        return self.gh.getIntervalsByKey()

    def annotate_is_UTR(self,chrom, strand, gpos):

        if strand:
            il = self.gh_utr.plusgenes[chrom]
        else:
            il = self.gh_utr.minusgenes[chrom]

        if il is None:
            return False

        il = il.genes
        matching_gene = [
            interval for interval in il if interval.start <= gpos < interval.end
        ]
        if len(matching_gene) > 0:

            return True

        return False



    def annotate_genomic_position(self,chrom, strand, gpos):

        if strand:
            il = self.gh.plusMaingenes[chrom]
        else:
            il = self.gh.minusMaingenes[chrom]


        matching_gene = [
            interval for interval in il if interval.start <= gpos < interval.end
        ]
        if len(matching_gene) > 0:

            if len(matching_gene)==1:
                best_gene = matching_gene[0]
            else:
                overlapping_genes_sorted = sorted(matching_gene, key=get_level)
                best_gene = overlapping_genes_sorted[0]

            print(best_gene)
            gene_id = best_gene.attrs.get('ID', '')
            gene_type = best_gene.attrs.get('gene_type', '')
            gene_name = best_gene.attrs.get('gene_name', '')

            return "ID="+gene_id+",gene_name="+gene_name+",gene_type="+gene_type

        return None


ref = "/mnt/ssdnas07/pipeline/rna_v08/source/mm10.fa"
gtf_file = '/mnt/ssdnas07/pipeline/rna_v08/source/gencode.vM25.annotation.gff3'
# stringtie_gtf = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/Adipocyte_2/Adipocyte_2/Adipocyte_2out.gtf"
stringtie_gtf = None

annotator = GenomeAnnotator(ref,gtf_file,stringtie_gtf)
chrom="chr6"
strand=True
gpos = 145216699
145220737
annoret = annotator.annotate_genomic_position(chrom, strand, gpos)
print(annoret)
isUTR = annotator.annotate_is_UTR(chrom, strand, gpos)
print(isUTR)
