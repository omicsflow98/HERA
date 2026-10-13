process plink {
        ext version: '1.0.0'
        
        label 'plink'

        publishDir "${params.main.outdir}/output/plink"

        input:
        tuple val(type), path(vcf), path(bed), path(bim), path(fam)
        val plink_map
        val chrom

        output:
        path("plink.*"), emit: plink_files
        path("plink.bim"), emit: bim_file
        path("plink.hom"), emit: roh_file
        val(chrom), emit: chromosome
        tuple path("plink.map"), path("plink.ped"), emit: dr_files
        path("${task.process}_versions.yml"), emit: plink_version

        script:

        def cleanname = task.process.split(':')[-1]
        def command
        def make_bed
        def stage = ""

        if (type == 'VCF') {
          command = "--vcf ${vcf}"
          make_bed = true
        } else {
          command = "--bfile plink"
          make_bed = false

          def link = { f, ext ->
            f.toString() == "plink.${ext}" ? "" : "ln -s ${f} plink.${ext}"
          }

          stage = [
            link(bed, 'bed'),
            link(bim, 'bim'),
            link(fam, 'fam')
            ].findAll { it }.join(' && ')
        }

        """
        ${stage}

        plink \
        ${command} \
        --homozyg \
        --chr-set 35 no-xy no-mt \
        ${ make_bed ? "--make-bed" : "" } \
        ${ plink_map.snp ? "--homozyg-snp ${plink_map.snp}" : "" } \
        ${ plink_map.kb ? "--homozyg-kb ${plink_map.kb}" : "" } \
        ${ plink_map.density ? "--homozyg-density ${plink_map.density}" : "" } \
        ${ plink_map.gap ? "--homozyg-gap ${plink_map.gap}" : "" } \
        ${ plink_map.windowsnp ? "--homozyg-window-snp ${plink_map.windowsnp}" : "" } \
        ${ (plink_map?.containsKey('WindowHet') && plink_map.WindowHet != null) ? "--homozyg-window-het ${plink_map.WindowHet}" : "" } \
        ${ (plink_map?.containsKey('WindowMissing') && plink_map.WindowMissing != null) ? "--homozyg-window-missing ${plink_map.WindowMissing}" : "" } \
        ${ plink_map.WindowThreshold ? "--homozyg-window-threshold ${plink_map.WindowThreshold}" : "" } \
        ${ (plink_map?.containsKey('HomozygHet') && plink_map.HomozygHet != null) ? "--homozyg-het ${plink_map.HomozygHet}" : "" } \
        ${ plink_map.nosex ? "--allow-no-sex" : "" } \
        --allow-extra-chr \
        ${ chrom ? "--chr ${chrom}" : "" } \
        --recode \
        --threads 2 \
        --out plink

        cut -f1 plink.bim | uniq > chr-map.txt

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          plink: \$(plink --version | awk '{print \$2}')
        EOF
        
        """
}

process plinkhet {
        ext version: '1.0.0'
        
        label 'plink'

        publishDir "${params.main.outdir}/output/plinkhet"

        input:
        tuple val(type), path(vcf), path(bed), path(bim), path(fam)
        val chrom

        output:
        tuple path("plink.map"), path("plink.ped"), emit: dr_files
        path("${task.process}_versions.yml"), emit: plink_version
        val(chrom), emit: chromosome

        script:
        def prefix = bim.toString().replace('.bim', '')
        def cleanname = task.process.split(':')[-1]

        """
        plink \
        --bfile ${prefix} \
        --chr ${chrom} \
        --chr-set 35 \
        --recode \
        --allow-extra-chr \
        --out plink

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          plink: \$(plink --version | awk '{print \$2}')
        EOF
        """
}