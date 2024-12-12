import click

@click.group()
def cmd():
    pass

# from utils.GtfToFasta import process
# @cmd.command()
# @click.option('-gtf', '--gtf')
# @click.option('-gtf_stringtie', '--gtf_stringtie')
# @click.option('-ref', '--ref')
# @click.option('-outmatrix', '--outmatrix')
# @click.option('-outfasta', '--outfasta')
# def gtfToFasta(gtf,gtf_stringtie,ref,outmatrix,outfasta):
#
#     process(gtf,gtf_stringtie,ref,outmatrix,outfasta)
#
#
# from mapping.BamConvert import run_convert
# @cmd.command()
# @click.option('-inbam', '--inbam')
# @click.option('-outbam', '--outbam')
# @click.option('-ref', '--ref')
# @click.option('-convertmatrix', '--convertmatrix')
# def bamConvert(inbam, outbam, ref, convertmatrix):
#
#     run_convert(inbam, outbam, ref, convertmatrix)

from recalibration.Recalibration import run_recalib
@cmd.command()
@click.option('-inbam', '--inbam', required=True)
@click.option('-outbam', '--outbam', required=True)
@click.option('-ref', '--ref', required=True)
@click.option('-recalib_db', '--recalib_db', required=True)
@click.option('-out_stats', '--out_stats', required=True)
def recalib(inbam, outbam, ref, recalib_db, out_stats):

    print("recalib")
    run_recalib(inbam, outbam, ref, recalib_db, out_stats)


from pileup.PileUPAll import pileup_all
@cmd.command()
@click.option('-prop', '--property', required=True)
@click.option('-recalib_stats','--recalib_stats', required=True)
@click.option('-inbam', '--inbam', required=True)
@click.option('-outdir', '--outdir', required=True)
@click.option('-ref', '--ref', required=True)
@click.option('-gtf_file', '--gtf_file', required=True)
@click.option('-stringtie_gtf', '--stringtie_gtf')
@click.option('-ncore', '--ncore', required=True)
def pileup(property,recalib_stats,inbam,outdir, ref, gtf_file, stringtie_gtf,ncore):

    # pileup_all(yamlf, recalib_stats, bamfile_name, outdir, gtf_file, ref, stringtie_gtf, ncore=8):
    pileup_all(property, recalib_stats, inbam, outdir, ref, gtf_file,  stringtie_gtf, ncore=int(ncore))


from filter.FurtherFilter import classification
@cmd.command()
@click.option('-invcf', '--invcf')
@click.option('-outvcf', '--outvcf')
@click.option('-nn_wight', '--nn_wight')
@click.option('-knowndir', '--knowndir')
@click.option('-genome', '--genome')
def filter(invcf,outvcf,nn_wight,knowndir,genome):

    classification(invcf, outvcf, nn_wight,knowndir,genome)

from summary.SummaryAll import summary
@cmd.command()
@click.option('-vcf', '--vcf')
@click.option('-gtf', '--gtf')
@click.option('-out', '--out')
def summaryall(vcf,gtf,out):

   summary(vcf, gtf, out)


from filter.AttentionClassfication import trainNN
@cmd.command()
@click.option('-m6Apath', '--m6Apath')
@click.option('-m5Cpath', '--m5Cpath')
@click.option('-psudepath', '--psudepath')
@click.option('-fp_ivtpath', '--fp_ivtpath')
@click.option('-ref', '--ref')
@click.option('-weightpath', '--weightpath')
def train(m6Apath,m5cpath,psudepath,fp_ivtpath,ref,weightpath):

   trainNN(m6Apath,m5cpath,psudepath,fp_ivtpath,ref,weightpath)

if __name__ == "__main__":
    cmd()

