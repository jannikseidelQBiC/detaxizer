process SUMMARIZER {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
        path(tosummarize)
        

    output:
        path("summary.tsv"), emit: summary
        //path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env python
    import glob
    files_kraken2 = glob.glob('*.kraken2_summary.tsv')
    files_blastn = glob.glob('*.blastn_summary.tsv')

    with open("summary.tsv", 'w') as f:
        for file in files_kraken2:
            f.write(file + '\\n')
        for file in files_blastn:
            f.write(file + '\\n')
    """
}