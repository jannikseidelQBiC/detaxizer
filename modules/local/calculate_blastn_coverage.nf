process CALCULATE_BLASTN_COVERAGE {

    container = 'https://depot.galaxyproject.org/singularity/python%3A3.10.4'

    input:
        tuple val(meta), path(blast_output)

    output:
        tuple val(meta), path('*coverage.txt'), emit: classified
        path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python

    def calculate_coverage(query_length, alignment_length):
        return (alignment_length / query_length) * 100

    with open('$blast_output', 'r') as f, open('${meta.id}.coverage.txt', 'w') as out:
        for line in f:
            parts = line.strip().split('\\t')
            query_id = parts[0]
            identity = float(parts[2])
            alignment_length = int(parts[3])
            start = int(parts[6])
            end = int(parts[7])
            query_length = end - start + 1 # Assuming the provided data covers the entire query
            coverage = calculate_coverage(query_length, alignment_length)
            if identity >= $params.blast_similarity and coverage >= $params.blast_coverage:
                out.write(f"{query_id}\\tIdentity: {identity:.2f}%\\tCoverage: {coverage:.2f}%\\n")

    import subprocess
    def get_version():
        version_output = subprocess.getoutput('python --version')
        version = version_output.split()[1]
        return version

    with open('versions.yml', 'w') as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f'    python: {get_version()}\\n')
    """

}