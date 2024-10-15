import click

@click.group()
def cmd():
    pass

from utils.GtfToFasta import process
@cmd.command()
@click.option('-gtf', '--gtf')
@click.option('-gtf_stringtie', '--gtf_stringtie')
@click.option('-ref', '--ref')
@click.option('-outmatrix', '--outmatrix')
@click.option('-outfasta', '--outfasta')
def gtfToFasta(gtf,gtf_stringtie,ref,outmatrix,outfasta):

    process(gtf,gtf_stringtie,ref,outmatrix,outfasta)


from mapping.BamConvert import run_convert
@cmd.command()
@click.option('-inbam', '--inbam')
@click.option('-outbam', '--outbam')
@click.option('-ref', '--ref')
@click.option('-convertmatrix', '--convertmatrix')
def bamConvert(inbam, outbam, ref, convertmatrix):

    run_convert(inbam, outbam, ref, convertmatrix)

from recalibration.Recalibration import run_recalib
@cmd.command()
@click.option('-inbam', '--inbam')
@click.option('-outbam', '--outbam')
@click.option('-ref', '--ref')
@click.option('-recalib_db', '--recalib_db')
@click.option('-out_stats', '--out_stats')
def recalib(inbam, outbam, ref, recalib_db, out_stats):

    run_recalib(inbam, outbam, ref, recalib_db, out_stats)


from pileup.PileUpAll import pileup_all
@cmd.command()
@click.option('-prop', '--property')
@click.option('-inbam', '--inbam')
@click.option('-convertmatrix', '--convertmatrix')
@click.option('-outdir', '--outdir')
@click.option('-gtf_file', '--gtf_file')
@click.option('-ref', '--ref')
@click.option('-dbsnpdir', '--dbsnpdir')
@click.option('-errormatrix_m6a', '--errormatrix_m6a')
@click.option('-errormatrix_all', '--errormatrix_all')
@click.option('-ncore', '--ncore')
# yamlf,bamfile_name,convertMatrix,outdir,gtf_file, ref,dbSNPdir,errorMatrix_m6A,errorMatrix_All
def pileup(property,inbam,convertmatrix,outdir,gtf_file, ref,dbsnpdir,errormatrix_m6a,errormatrix_all,ncore):

    pileup_all(property,inbam,convertmatrix,outdir,gtf_file, ref,dbsnpdir,errormatrix_m6a,errormatrix_all,ncore=int(ncore))



from filter.FurtherFilter import classification
@cmd.command()
@click.option('-invcf', '--invcf')
@click.option('-outvcf', '--outvcf')
@click.option('-nn_wight', '--nn_wight')
@click.option('-m5cpath', '--m5cpath')
@click.option('-psudepath', '--psudepath')
def filter(invcf,outvcf,nn_wight):

    classification(invcf, outvcf, nn_wight)


