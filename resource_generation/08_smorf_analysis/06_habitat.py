
def map_habitat(infile1,infile2,outfile):
    habitat_dict = {}
    with open(infile1,'rt') as f:
        for line in f:
            smorf,habitat = line.strip().split('\t')
            smorf = smorf.replace('_contig','')
            habitat_dict[smorf] = habitat
    with open(outfile,'wt') as out:
        with open(infile2,'rt') as f:
            for line in f:
                habitat_set = set()
                linelist = line.strip().split('\t')
                smorfs = linelist[3].split(',')
                for item in smorfs:
                    parts = item.split('#')
                    header = f'{parts[0]} # {parts[1]}'
                    habitats = habitat_dict[header].split(',')
                    for h in habitats:
                        habitat_set.add(h)
                length = len(habitat_set)
                final = ','.join(sorted(list(habitat_set)))
                out.write(f'{linelist[0]}\t{final}\t{length}\n')

infile1 = 'all_habitat.out.smorfs_newname.tsv'
infile2 = 'cluster_count.tsv'
outfile = 'cluster_habitat.tsv'
map_habitat(infile1,infile2,outfile)

def general(infile,outfile):
    higher_level = {
            'fermented food' : 'anthropogenic',
            'activated sludge' : 'anthropogenic',
            'wastewater' : 'anthropogenic',
            'built environment' : 'anthropogenic',
            'anthropogenic': 'anthropogenic',
            'groundwater' : 'aquatic',
            'river associated' : 'aquatic',
            'lake associated' : 'aquatic',
            'water associated' : 'aquatic',
            'marine' : 'aquatic',
            'pond associated' : 'aquatic',
            'plant associated' : 'soil/plant',
            'soil' : 'soil/plant',
            'bird gut' : 'other animal',
            'chicken gut' : 'other animal',
            'cattle rumen' : 'other animal',
            'bee gut' : 'other animal',
            'dog associated' : 'other animal',
            'cattle associated' : 'other animal',
            'insect gut' : 'other animal',
            'crustacean associated' : 'other animal',
            'planarian associated' : 'other animal',
            'sponge associated' : 'other animal',
            'goat rumen' : 'other animal',
            'crustacean gut' : 'other animal',
            'annelidae associated' : 'other animal',
            'bird skin' : 'other animal',
            'beatle gut' : 'other animal',
            'termite gut' : 'other animal',
            'fish gut' : 'other animal',
            'tunicate associated' : 'other animal',
            'mussel associated' : 'other animal',
            'mollusc associated' : 'other animal',
            'ship worm associated' : 'other animal',
            'wasp gut' : 'other animal',
            'insect associated' : 'other animal',
            'coral associated' : 'other animal',
            'turtle gut' : 'other animal',
            'human urogenital tract' : 'other human',
            'human associated' : 'other human',
            'human respiratory tract' : 'other human',
            'human skin' : 'other human',
            'human digestive tract' : 'other human',
            'human saliva' : 'other human',
            'human mouth' : 'other human',
            'human gut' : 'human gut',
            'isolate' : 'isolate',
            'dog gut' : 'mammal gut',
            'cat gut' : 'mammal gut',   
            'rat gut' : 'mammal gut',
            'cattle gut' : 'mammal gut',
            'deer gut' : 'mammal gut',
            'mouse gut' : 'mammal gut',
            'primate gut' : 'mammal gut',
            'pig gut' : 'mammal gut',
            'bear gut' : 'mammal gut',
            'bat gut' : 'mammal gut',
            'goat gut' : 'mammal gut',
            'rodent gut' : 'mammal gut',
            'fisher gut' : 'mammal gut',
            'coyote gut' : 'mammal gut',
            'rabbit gut' : 'mammal gut',
            'horse gut' : 'mammal gut',
            'guinea pig gut' : 'mammal gut',
            'dolphin gut' : 'mammal gut',
            'whale gut' : 'mammal gut'
            }
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                cluster,habitat,n = line.strip().split('\t')
                habitat_list = habitat.split(',')
                general = set()
                for item in habitat_list:
                    if item in higher_level.keys() and item !='isolate':
                        general.add(higher_level[item])
                    else:
                        general.add('other')
                final = ','.join(sorted(list(general)))
                out.write(f'{line.strip()}\t{final}\t{len(general)}\n')

infile = 'cluster_habitat.tsv'
outfile = 'cluster_habitat_general.tsv'
general(infile,outfile)


def cal_sample(infile,outfile):
    with open(outfile,'wt') as out:
        with open(infile,'rt') as f:
            for line in f:
                sample_set = set()
                linelist = line.strip().split('\t')
                smorfs = linelist[3].split(',')
                for item in smorfs:
                    sample = item.split('_')[0]
                    sample_set.add(sample)
                length = len(sample_set)
                out.write(f'{linelist[0]}\t{length}\n')

infile = 'cluster_count.tsv'
outfile = 'cluster_sample.tsv'
cal_sample(infile,outfile)

def cal_fraction(infile,outfile):
    habitat_dict = {}
    with open(infile,'rt') as f:
        for line in f:
            linelist = line.strip().split('\t')
            habitatlist = linelist[3].split(',')
            for item in habitatlist:
                if item not in habitat_dict.keys():
                    habitat_dict[item] = 1
                else:
                    habitat_dict[item] += 1
    with open(outfile,'wt') as out:
        for key,value in habitat_dict.items():
            out.write(f'{key}\t{value}\n')

infile = 'cluster_habitat_general.tsv'
outfile = 'cluster_habitat_general_fraction_general.tsv'
cal_fraction(infile,outfile)