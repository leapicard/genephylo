/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// BLAST database
include { BLAST_UPDATEBLASTDB } from '../modules/nf-core/blast/updateblastdb/main'
include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb/main'

// retrieve sequences
include { BLAST_BLASTN } from '../modules/nf-core/blast/blastn/main'
include { BLAST_TBLASTN } from '../modules/nf-core/blast/tblastn/main'
include { BLAST_EXTRACT } from '../modules/local/blast_extract'
include { BLAST_BLASTDBCMD } from '../modules/nf-core/blast/blastdbcmd/main'
include { SEQKIT_RMDUP } from '../modules/nf-core/seqkit/rmdup/main'
include { BLAST_FILTER } from '../modules/local/blast_filter'

// ALIGN AND TREE
include { MAFFT_ALIGN } from '../modules/nf-core/mafft/align/main'
include { IQTREE } from '../modules/nf-core/iqtree/main'
include { FASTTREE } from '../modules/nf-core/fasttree/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genephylo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENEPHYLO {

		take:
		ch_samplesheet // channel: samplesheet read in from --input

		main:

		ch_versions = channel.empty()

		//
		// SUBWORKFLOW: build BLAST database for subsequent analyses (3 options)
		//

		ch_blastdb_in = channel.empty()
		ch_blastdb = channel.empty()

		// 1. download from ncbi databases (update)
		if (params.blastdb_option == 'update') {

			ch_blastdb_in = channel.of( tuple( [ id: "BLASTDB" ], params.blastdb_update ))

			BLAST_UPDATEBLASTDB( ch_blastdb_in )
			
			// decompression step added through modules.config

			ch_blastdb = BLAST_UPDATEBLASTDB.out.db
			ch_versions = ch_versions.mix(BLAST_UPDATEBLASTDB.out.versions)

		}
		// 2. build custom databases from sequences (build)
		else if (params.blastdb_option == 'build') {

			channel
				.fromPath(params.blastdb_build, checkIfExists: true)
				.map { dbpath -> tuple( [ id: "BLASTDB" ], file(dbpath) ) }
				.set { ch_blastdb_in }

			BLAST_MAKEBLASTDB( ch_blastdb_in )

			ch_blastdb = BLAST_MAKEBLASTDB.out.db
			ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

		}
		// 3. use existing database (current)
		else if (params.blastdb_option == 'current') {

			ch_blastdb = channel.of(tuple([ id: "BLASTDB" ], file(params.blastdb_dir) ))

		}

		//
		// SUBWORKFLOW: run BLAST on target sequences (2 options)
		//
		
		ch_blast_in = ch_samplesheet.map { meta, fasta, tree -> tuple(meta, file(fasta)) }
		ch_blastdb_in = ch_blastdb.first()
		ch_blast_out = channel.empty()

		// 1. nt sequence: blastn
		if ( params.blast_type == 'nt' ) {

			BLAST_BLASTN( ch_blast_in, ch_blastdb_in )

			ch_blast_out = BLAST_BLASTN.out.txt
			ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

		}
		// 2. aa sequence: tblastn
		else if ( params.blast_type == 'aa' ) {

			BLAST_TBLASTN( ch_blast_in, ch_blastdb_in )

			ch_blast_out = BLAST_TBLASTN.out.txt
			ch_versions = ch_versions.mix(BLAST_TBLASTN.out.versions)

		}

		// get the first column of the blast results file
		BLAST_EXTRACT( ch_blast_out )

		ch_accessions_out = BLAST_EXTRACT.out.accessions
		ch_versions = ch_versions.mix(BLAST_EXTRACT.out.versions)

		ch_extract_in = ch_accessions_out.map { meta, batch_file ->
			tuple(meta, null, batch_file)
		}
		
		BLAST_BLASTDBCMD(ch_extract_in, ch_blastdb_in)

		ch_extract_out = BLAST_BLASTDBCMD.out.fasta
		ch_versions = ch_versions.mix(BLAST_BLASTDBCMD.out.versions)

		SEQKIT_RMDUP(ch_extract_out)

		ch_rmdup = SEQKIT_RMDUP.out.fastx
		ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

		BLAST_FILTER(ch_rmdup, ch_blast_out)

		ch_aln_in = BLAST_FILTER.out.fasta
		ch_versions = ch_versions.mix(BLAST_FILTER.out.versions)

		// SUBWORKFLOW: phylo
		//

		// build initial alignment
		MAFFT_ALIGN (
			ch_aln_in, [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ], [ [:], [] ], false
			)

		ch_mafft_out = MAFFT_ALIGN.out.fas
		ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions)

		// build phylogenetic tree
		FASTTREE ( ch_mafft_out )

		ch_fasttree_out = FASTTREE.out.phylogeny
		ch_versions = ch_versions.mix(FASTTREE.out.versions)

		ch_iqtree_in = ch_mafft_out
			.join(ch_fasttree_out)
			.map { meta, alignment, tree -> tuple(meta, alignment, tree) }

		IQTREE (
			ch_iqtree_in, [], [], [], [], [], [], [], [], [], [], [], []
			)
		
		ch_iqtree_out = IQTREE.out.phylogeny
		ch_versions = ch_versions.mix(IQTREE.out.versions)

		// if ( params.tree_tool == "iqtree" ) {
			
		// }
		// else if ( params.tree_tool == "fasttree" ) {
		// 	ch_fasttree_in = ch_mafft_out.map { meta, alignment -> alignment }
		// }

		softwareVersionsToYAML(ch_versions)
				.collectFile(
						storeDir: "${params.outdir}/pipeline_info",
						name: 'nf_core_'  +  'genephylo_software_'  + 'mqc_'  + 'versions.yml',
						sort: true,
						newLine: true
				).set { ch_collated_versions }


		emit:
		versions       = ch_versions.toList()                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
