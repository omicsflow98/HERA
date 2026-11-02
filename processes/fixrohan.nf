process rohan {

        ext version: '1.0.0'
	label 'fixrohan'

        input:
        path rohan_script
        path rohan_files

        output:
        path("rohan.hom"), emit: roh_file
        path("${task.process}_versions.yml"), emit: version

        script:

        def cleanname = task.process.split(':')[-1]
        
        """
        python3 ${rohan_script} \
        --rohan ${rohan_files} \
        --version ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml
        
        """
}