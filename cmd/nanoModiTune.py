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
@click.option('-inbam', '--inbam')
@click.option('-outbam', '--outbam')
@click.option('-ref', '--ref')
@click.option('-recalib_db', '--recalib_db')
@click.option('-out_stats', '--out_stats')
def recalib(inbam, outbam, ref, recalib_db, out_stats):

    run_recalib(inbam, outbam, ref, recalib_db, out_stats)


from pileup.PileUPAll import pileup_all
@cmd.command()
@click.option('-prop', '--property')
@click.option('-recalib_stats','--recalib_stats')
@click.option('-inbam', '--inbam')
@click.option('-outdir', '--outdir')
@click.option('-ref', '--ref')
@click.option('-gtf_file', '--gtf_file')
@click.option('-stringtie_gtf', '--stringtie_gtf')
@click.option('-ncore', '--ncore')
def pileup(property,recalib_stats,inbam,outdir, ref, gtf_file, stringtie_gtf,ncore):

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

from summary.SummaryAll_old import summary
@cmd.command()
@click.option('-vcf', '--vcf')
@click.option('-tpm', '--tpm')
@click.option('-genetsv', '--genetsv')
@click.option('-out', '--out')
def summaryall(vcf,tpm,genetsv,out):

   summary(vcf, tpm, genetsv, out)



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


