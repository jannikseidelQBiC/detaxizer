process ISOLATE_IDS_FROM_KRAKEN2_TO_BLASTN {
    input:
        tuple val(meta), path(kraken2results)

    output:
        tuple val(meta), path('*classified.txt'), emit: classified
        path "versions.yml", emit: versions

    script:
    """
    grep '${params.tax2filter}' $kraken2results > ${meta.id}.classified.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | grep -oP 'grep \\(GNU grep\\) \\K\\d+(\\.\\d+)*')
    END_VERSIONS
    """
}