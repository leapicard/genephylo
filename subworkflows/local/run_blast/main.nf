// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

// BLAST - make a subworkflow?
include { BLAST_UPDATEBLASTDB } from '../modules/nf-core/blast/updateblastdb/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
//include { BLAST_TBLASTN } from '../modules/nf-core/blast/tblastn/main'
include { BLAST_BLASTDBCMD } from '../modules/nf-core/blast/blastdbcmd/main'


workflow RUN_BLAST {

    take:
    // TODO nf-core: edit input (take) channels
    ch_samplesheet // channel: [ val(meta), [ bam ] ]
    ch_blastdb // channel; [ id:'mito', single_end:false ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    BLAST_UPDATEBLASTDB ( ch_blastdb )
    ch_versions = ch_versions.mix(BLAST_UPDATEBLASTDB.out.versions.first())

    BLAST_BLASTN ( ch_samplesheet, BLAST_UPDATEBLASTDB.out )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
