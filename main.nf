#!/usr/bin/env nextflow

include { ROH as ROHOM } from './subworkflows/ROH/ROH.nf'
include { ROH as ROHET } from './subworkflows/ROH/ROH.nf'
include { quarto } from './processes/quarto.nf'
include { plinkhet } from './processes/plink.nf'
include { plink as plinkhom } from './processes/plink.nf'
include { plink as plinkhetbim } from './processes/plink.nf'
include { detectruns } from './processes/detectruns.nf'
include { annotationtobed } from './processes/bedops.nf'
include { split_plink } from './processes/split_plink.nf'
include { all_versions } from './processes/MergeVersions.nf'
include { vcftools } from './processes/vcftools.nf'
include { rohan } from './processes/fixrohan.nf'

workflow {
	
    //Pipeline version
    pipeline_version = "1.1.0"

    // Read in the scripts and files needed for the workflow
    dr_script = file("${projectDir}/scripts/detectRUNS.r")
    rohan_script = file("${projectDir}/scripts/fixrohan.py")
    chromosomes= file( "${params.settings.chromsfile}" )
    split_results = file("${projectDir}/scripts/split_plink.py")
    version_script = file("${projectDir}/scripts/version_merge.py")
    subpop_tsv = params.settings.subpopulations? channel.fromPath(params.subpopulation.subpop_file):channel.of([])

    // Initialize empty channels for outputs
    ROHOM_out = channel.empty()
    ROHET_out = channel.empty()
    ROHOM_indiv = channel.empty()
    ROHET_indiv = channel.empty()
    ROHOM_islands = channel.empty()
    ROHET_islands = channel.empty()
    all_versions_ch = channel.empty()

    // Read in parameters from the config file
    dr_params = channel.value(params.detectruns.parameters)
    plink_params = channel.value(params.plink.parameters)
    if (params.inputfiles.filetype != 'rohan') {
        island_info = channel.value(
            tuple (
            params.island.type,
            params.island.minSNP,
            params.island.indiv_count
            )
        )
    } else {
        island_info = channel.value (
            tuple (
                'fill',
                0,
                0
            )
        )
    }


    // Read in files needed to generate final report
    quarto_files = channel.of(
    tuple(
        file("${projectDir}/scripts/Quarto/report.qmd"),
        file("${projectDir}/scripts/Quarto/test.R"),
        file("${projectDir}/scripts/Quarto/subpopulations.R"),
        file("${projectDir}/scripts/Quarto/style.scss")
        )
    )

    // Read in gene file if user wants to include gene information in the final report
    if (params.settings.genes) {
        def anno = file("${params.genes.file}", checkIfExists: true)
        def annoName = anno.name.toLowerCase()
        def annoFmt  = annoName.endsWith('.gtf')                              ? 'gtf'
             : (annoName.endsWith('.gff') || annoName.endsWith('.gff3')) ? 'gff'
             : annoName.endsWith('.bed')                              ? 'bed'
             : null
        if ( !annoFmt ) { error "Annotation file must end in .gtf, .gff, .gff3 or .bed: ${anno.name}" }

        if (annoFmt == "bed") {
            bed_file = channel.value(anno)
        } else {
            bed_file = annotationtobed( channel.value([annoFmt, anno]) ).bedfile
        }        
    } else {
        bed_file = channel.value([])
    }

    if (params.inputfiles.filetype != 'vcf' && params.inputfiles.filetype != 'plink' && params.inputfiles.filetype != 'rohan') {
        error "Input file type must be either vcf, plink, or rohan"
    }

    // Read in input files and determine if the user is starting with a VCF, PLINK files, or ROHan file. Extract VCF info if starting with a VCF file
    if (params.inputfiles.filetype == 'vcf') {
        infile = channel.value(
            tuple(
                params.inputfiles.filetype,
                file(params.inputfiles.file),
                [],
                [],
                []
            )
        )
        vcf_info = vcftools(file(params.inputfiles.file))
    } else if (params.inputfiles.filetype == 'plink') {
        prefix = params.inputfiles.file
        def bed = file("${prefix}.bed")
        def bim = file("${prefix}.bim")
        def fam = file("${prefix}.fam")
            
        infile = channel.value(
            tuple(
                'plink',
                [],
                bed,
                bim,
                fam
            )
        )
        vcf_info = [ vcf_info: file("${prefix}.fam") ]
    } else {
        rohan_files = channel.fromPath("${params.inputfiles.file}/*.hmmrohl", checkIfExists: true).collect()
        vcf_info = [ vcf_info: file("${params.inputfiles.file}") ]
    }

    //Perform homozygosity analysis if user has selected this option in the config file
    if (params.settings.homozygosity) {
        homstate = "HOM"

        if (params.inputfiles.filetype != "rohan") {
            plinkhom_map = plinkhom(infile, plink_params, "")
            bim_file = plinkhom_map.bim_file
        } else {
            plinkhom_map = rohan(rohan_script, rohan_files)
            bim_file = [ bim_file: file("${params.inputfiles.file}") ]
        }

        split_plink_map = split_plink(split_results, plinkhom_map.roh_file, chromosomes)

        ROHOM_map = ROHOM(split_plink_map.separated_plink, bim_file, bed_file, homstate, params.settings.genes, params.inputfiles.filetype)
        ROHOM_out = ROHOM_map.ROH_results
        ROHOM_version = ROHOM_map.all_versions
        ROHOM_indiv = ROHOM_map.indiv_info
        ROHOM_islands = ROHOM_map.finalislands
    }

    //Perform homozygosity analysis if user has selected this option in the config file
    if (params.settings.heterozygosity) {
        
        chroms = channel.fromList(chromosomes.readLines())
        hetstate = "HET"
        plinkhet_map = plinkhet(infile, chroms)
        plinkhet_version = plinkhet_map.version.collect().map{ it[0] }

        if (params.inputfiles.start_vcf) {
            make_bim = plinkhetbim(infile, plink_params, "")
        } else {
            make_bim = [ bim_file: file("${prefix}.bim") ]
        }

        detectruns_map = detectruns(dr_script, plinkhet_map.dr_files, dr_params, plinkhet_map.chromosome)

        detectruns_map.dr_tsv
            | collect
            | set {all_runs}

        detectruns_version = detectruns_map.dr_version.collect().map{ it[0] }

        ROHET_map = ROHET(all_runs, make_bim.bim_file, bed_file, hetstate, params.settings.genes, params.inputfiles.filetype)
        ROHET_out = ROHET_map.ROH_results
        ROHET_version = ROHET_map.all_versions
        ROHET_indiv = ROHET_map.indiv_info
        ROHET_islands = ROHET_map.finalislands
    }

    //Collect all output information from subworkflows
    ROH_all = ROHOM_out.mix(ROHET_out).collect()
    ROH_island = ROHOM_islands.mix(ROHET_islands).collect()
    indiv_all = ROHOM_indiv.mix(ROHET_indiv).collect()

    //Use all information to generate final report using QUARTO
    quarto_map = quarto(ROH_all,
    ROH_island, 
    indiv_all, 
    vcf_info.vcf_info, 
    quarto_files, 
    params.settings.homozygosity, 
    params.settings.heterozygosity, 
    params.settings.subpopulations, 
    params.inputfiles.filetype,
    params.settings.genes,
    subpop_tsv, 
    bed_file,
    island_info )

    //populate version information channel
    all_versions_ch = quarto_map.quarto_version

    if (params.settings.homozygosity) {
        all_versions_ch = all_versions_ch
        | combine(ROHOM_version)
        | combine(plinkhom_map.version)
        | combine(split_plink_map.splitplink_version)
    }
    if (params.settings.heterozygosity) {
        all_versions_ch = all_versions_ch
        | combine(ROHET_version)
        | combine(plinkhet_map.version)
        | combine(detectruns_version)
    }

    //output version information to a text file
    all_versions(version_script, all_versions_ch, pipeline_version)
}
