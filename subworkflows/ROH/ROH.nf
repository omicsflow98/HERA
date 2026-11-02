#!/usr/bin/env nextflow

include { overlap } from '../../processes/overlap.nf'
include { combine_overlap } from '../../processes/overlap.nf'
include { merge_files } from '../../processes/mergefiles.nf'
include { splitROH } from '../../processes/splitROH.nf'
include { fixROH } from '../../processes/fixROH.nf'
include { snpislands } from '../../processes/snpislands.nf'
include { topsnp } from '../../processes/topsnp.nf'

workflow ROH {

	take:
	collected_tsvs
	bim_file
	bed_file
	state
	genes_enabled
	filetype

	main:

	chromosomes= file( "${params.settings.chromsfile}" )
	overlap_script = file("${projectDir}/scripts/overlap_roh.py")
	split_script = file("${projectDir}/scripts/splitROH.R")
	fix_script = file("${projectDir}/scripts/fixROH.R")
	topsnp_script = file("${projectDir}/scripts/topsnp.py")
	overlap_version = channel.empty()

	island_info = channel.value(
		tuple(
			params.island.type,
			params.island.threshold
		)
	)

	genome = params.settings.genomelength

	if (genes_enabled) {
			collected_tsvs
			| flatten()
			| set { individual_tsvs }
		overlap_map = overlap(overlap_script, individual_tsvs, bed_file)
		overlap_version = overlap_map.OL_version.collect().map{ it[0] }
		overlap_map.overlap_bed
			| collect
			| set { all_tsvs }
		combine_overlap_map = combine_overlap(all_tsvs, chromosomes)

	} else {
		combine_overlap_map = combine_overlap(collected_tsvs, chromosomes)
	}

	if (filetype != "rohan") {
		islands = snpislands(combine_overlap_map.ROH_file, bim_file)
		island_version = islands.snpislands_version
		topislands = topsnp(state, topsnp_script, islands.snpcount, island_info)
		topsnps_version = topislands.topsnps_version
		finalislands = topislands.topsnps
	} else {
		finalislands = channel.value(file(chromosomes))
		topsnps_version = channel.empty()
		island_version = channel.empty()
	}

	splitROH_map = splitROH(split_script, combine_overlap_map.ROH_file, state, genome)
	indiv_info = splitROH_map.indiv_file

	splitROH_map.split_tsv
		| flatten()
		| set { individual_ROH }

	fixROH_map = fixROH(fix_script, individual_ROH)

	fixROH_version = fixROH_map.fixROH_version.collect().map{ it[0] }

	fixROH_map.fixed_tsv
		| collect
		| set { fixed_tsv }

	merge_map = merge_files(fixed_tsv, state)

	ROH_results = merge_map.tsv_file

	version_info = fixROH_version
		| combine(combine_overlap_map.combineOL_version)
		| combine(splitROH_map.splitROH_version)
		| combine(topsnps_version)
		| combine(island_version)

	all_versions = version_info.mix(overlap_version)

	emit:
	ROH_results
	finalislands
	all_versions
	indiv_info
}
