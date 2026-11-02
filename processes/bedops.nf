process annotationtobed {

        ext version: '1.0.9'
	label 'bedops'

        publishDir "${params.main.outdir}/output/bedfile", mode: 'copy', overwrite: true

        input:
        tuple val(format), path(annofile)

        output:
        path("genes.bed"), emit: bedfile

        script:

        """
        convert2bed -i ${format} \
        < ${annofile} \
        > temp.bed

        awk '\$8 ~ /gene/' temp.bed > genes.bed
        
        """
}
