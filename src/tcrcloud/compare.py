import json
from collections import defaultdict

import tcrcloud.colours
import tcrcloud.format

# Import default colours based on the V gene
TRAV = tcrcloud.colours.TRAV
TRBV = tcrcloud.colours.TRBV
TRGV = tcrcloud.colours.TRGV
TRDV = tcrcloud.colours.TRDV


def format(args):
    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_cloud(samples_df)
    samples = formatted_samples.groupby(['chain', 'repertoire_id'])
    keys = [key for key, _ in samples]
    new = {}
    for j in keys:
        df = samples.get_group(j)

        if len(df) > 1:
            # create a dict associating the junction to the gene
            # with highest counts
            family = df[['junction_aa', 'v_call']].drop_duplicates(
                subset='junction_aa',
                keep="first",
                inplace=False).set_index('junction_aa').squeeze().to_dict()
        else:
            family = {df['junction_aa'].iloc[0]: df['v_call'].iloc[0]}

        for i in family:
            new[i] = eval(family.get(i)[:4]).get(family.get(i), [])
    return(new)


def compare(args):
    args.rearrangements = args.file1
    list1 = format(args)
    args.rearrangements = args.file2
    list2 = format(args)
    set1 = set(list1)
    set2 = set(list2)
    newlist1 = set1.difference(set2)
    newlist2 = set2.difference(set1)
    newlist3 = set1.intersection(set2)
    exclusive1 = {}
    exclusive2 = {}
    common = {}
    for i in newlist1:
        exclusive1[i] = list1.get(i)
    for i in newlist2:
        exclusive2[i] = list2.get(i)
    for i in newlist3:
        common[i] = list1.get(i)

    exclusive1_inverted = defaultdict(list)
    {exclusive1_inverted[v].append(k) for k, v in exclusive1.items()}
    exclusive1 = dict(exclusive1_inverted)

    exclusive2_inverted = defaultdict(list)
    {exclusive2_inverted[v].append(k) for k, v in exclusive2.items()}
    exclusive2 = dict(exclusive2_inverted)

    common_inverted = defaultdict(list)
    {common_inverted[v].append(k) for k, v in common.items()}
    common = dict(common_inverted)

    with open(args.file1[:-4] + '.json', 'w') as fileout:
        final = json.dumps(exclusive1)
        print(final, file=fileout)
        print("colours saved as " + args.file1[:-4] + '.json')
    with open(args.file2[:-4] + '.json', 'w') as fileout:
        final = json.dumps(exclusive2)
        print(final, file=fileout)
        print("colours saved as " + args.file2[:-4] + '.json')
    with open(args.file1[:-4] + "+" + args.file2[:-4]
              + '.json', 'w') as fileout:
        final = json.dumps(common)
        print(final, file=fileout)
        print("colours saved as " + args.file1[:-4] + "+"
              + args.file2[:-4] + '.json')
