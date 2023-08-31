process PARSE_KRAKEN2REPORT {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"
    
    input:
        tuple val(meta), path(kraken2report)
    
    output:
        path "versions.yml", emit: versions
        //path "*.txt", emit: txt
        path "*.json", emit: json
    
    script:
    """
    #!/usr/bin/env python
import json
def read_in_kraken2report(report):
    kraken_list = []
    with open(report, 'r') as f:
        for line in f:
            if line == "":
                continue
            line = line.strip('\\n').split('\\t')
            kraken_list.append(line)
    return kraken_list

def kraken_taxonomy2hierarchy(kraken_list):
    tax_dict = {}
    stack = []
    for line in kraken_list:
        level = (len(line[5])-len(line[5].strip()))/2
        node = line[5].strip()
        id = line[4]

        while len(stack) > level:
            stack.pop()
        
        if level == 0 or not stack:
            tax_dict[node] = {"id": id}
            stack = [tax_dict[node]]
        else:
            current_level = stack[-1]
            current_level[node] = {"id": id}
            stack.append(current_level[node])
    with open('${meta.id}.json', 'w') as f:
        json.dump(tax_dict,f)
    return tax_dict
        
def get_all_keys(nested_dict, parent_key='', sep='\\t'):
    keys = []

    for k, v in nested_dict.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        # If the current key is 'id', it's a leaf node and we add the parent_key to keys.
        if k == 'id':
            keys.append(parent_key)
        # If the value v is a dictionary, we need to traverse further.
        elif isinstance(v, dict):
            keys.extend(get_all_keys(v, new_key, sep=sep))
        
    return keys


def remove_incomplete_taxa(list_kraken):
    list_kraken.sort()
    to_remove = set()
    for i in range(len(list_kraken) - 1):
        # If the current item (plus a comma, to prevent removing "A,B,C" when there's "A,B,C,D") 
        # is a prefix of the next one, mark it for removal.
        if list_kraken[i+1].startswith(list_kraken[i] + ","):
            to_remove.add(i)

    # Create a new list without the marked indexes
    result = [item for idx, item in enumerate(list_kraken) if idx not in to_remove]

    return result

def generate_dict_for_lookup(removed_list):
    lookup_dict = {}
    for entry in removed_list:
        new_list = entry.split('\\t')
        while new_list:
            new_entry = new_list.pop(0)
            if new_entry in lookup_dict.keys():
                lookup_dict[new_entry].update(new_list)
            else:
                lookup_dict[new_entry] = set(new_list)
    for key in lookup_dict.keys():
        lookup_dict[key] = sorted(list(lookup_dict[key]))
    sorted_keys = sorted(lookup_dict.keys())
    sorted_dict = {k: lookup_dict[k] for k in sorted_keys}
    return sorted_dict

kraken_list = read_in_kraken2report('$kraken2report')
kraken_dict = kraken_taxonomy2hierarchy(kraken_list)
list_kraken = get_all_keys(kraken_dict)
result = remove_incomplete_taxa(list_kraken)
result = generate_dict_for_lookup(result)
result_to_filter = result["$params.tax2filter"] + ["$params.tax2filter"]
with open('test.json', "w") as f:
    json.dump(result_to_filter, f)
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