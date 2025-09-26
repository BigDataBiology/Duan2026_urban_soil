import os
from fasta import fasta_iter

folder_path = '~/soil/13_compare/marker/marker_spire'
output_file = 'spire_cog0495.fa'

with open(output_file, 'wt') as out:
    for filename in os.listdir(folder_path):
        if filename.endswith('fetchMGs.fna'):
            file_path = os.path.join(folder_path, filename)
            for h,seq in fasta_iter(file_path):
                if 'COG0495' in h:
                    head = f"{filename}#{h}"
                    out.write(f'>{head}\n{seq}\n')