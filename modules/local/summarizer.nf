process SUMMARIZER {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"
    input:
        tuple val(meta),path(tosummarize)
        

    output:
        path("summary.tsv"), emit: summary
        //path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env python
    tosummarize_list = "${tosummarize}".split(' ')
    for entry in tosummarize_list:
        print(entry)
    with open("summary.tsv", 'w') as f:
        f.write("test")
    """
}