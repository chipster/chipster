def read_input_definitions():
    input_names = {}
    with open('chipster-inputs.tsv') as inputs:
        for line in inputs:
            # remove \n
            line = line[:-1]
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            input_name = columns[0]
            dataset_name = columns[1]
            input_names[input_name] = dataset_name
    
    return input_names

def write_output_definitions(output_names):
    with open('chipster-outputs.tsv', 'w') as outputs:
        for output_name in output_names:
            dataset_name = output_names[output_name]
            outputs.write(output_name + '\t' + dataset_name + '\n')

def remove_postfix(string, postfix):
    if string.endswith(postfix):
        return string[:-len(postfix)]
    return string

