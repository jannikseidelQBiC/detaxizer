process SUMMARY_KRAKEN2 {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"
    input:
        tuple val(meta),path(kraken2),path(isolatedkraken2)
        

    output:
        tuple val(meta), path("*.kraken2_summary.tsv"), emit: summary
        //path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    lines_in_kraken2 = 0
    lines_in_isolatedkraken2 = 0
    
    with open("${kraken2}", 'r') as f:
        for line in f:
            if line != "\\n":
                lines_in_kraken2 += 1
    
    with open("${isolatedkraken2}", 'r') as f:
        for line in f:
            if line != "\\n":
                lines_in_isolatedkraken2 += 1
    
    data = {
        "kraken2": [lines_in_kraken2],
        "isolatedkraken2": [lines_in_isolatedkraken2]
        }
    
    df = pd.DataFrame(data)

    df.index = ["${meta.id}"]

    df.to_csv("${meta.id}.kraken2_summary.tsv",sep='\\t')
    """
}