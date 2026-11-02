process quarto {
        ext version: '1.0.3'

        label 'quarto'

        publishDir = [
          [
            path: "${params.main.outdir}/output/report",
            mode: 'copy',
            pattern: '*.zip'
          ],
          [
            path: "${params.main.outdir}/output/report/figures",
            mode: 'copy',
            pattern: '*.png'
          ],
          [
            path: "${params.main.outdir}/output/report/tables",
            mode: 'copy',
            pattern: '*.out.tsv'
          ]
        ]

        input:
        path(ROHall)
        path(ROHisland)
        path(indiv_all)
        path(vcf_info)
        tuple path(quartofile), path(plotsfile), path(subpopfile), path(stylesfile)
        val(homstate)
        val(hetstate)
        val(substate)
        val(filetype)
        val(genestate)
        path(subpop)
        path(bedfile)
        tuple val(islandtype), val(minsnp), val(indiv_count)

        output:
        path("*.zip")
        path("*.out.tsv")
        path("*.png")
        path("${task.process}_versions.yml"), emit: quarto_version

        script:

        def cleanname = task.process.split(':')[-1]
        
               
        """
        export HOME=\${PWD}

        if [[ "${genestate}" == "true" ]]; then
          mv ${bedfile} tempfilename
          mv tempfilename genes.bed

        else
          touch gene_fill.txt
        fi

        ${ subpop ? "mv ${subpop} subfile.tsv" : "" }

        quarto render ${quartofile} \\
        -P hom:${homstate} \\
        -P het:${hetstate} \\
        -P substate:${substate} \\
        -P filetype:${filetype} \\
        -P genestate:${genestate} \\
        -P type:${islandtype} \\
        -P minsnp:${minsnp} \\
        -P indiv:${indiv_count}

        zip -r report.zip report.html report_files

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          Quarto: \$(quarto --version)
	      EOF

        """
}
