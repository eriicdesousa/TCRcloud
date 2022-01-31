#import glob, airr, matplotlib, json, math
import json
import math

import airr
import matplotlib
matplotlib.use('Agg')
from wordcloud import (WordCloud)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
import numpy as np
import pandas as pd


# This colours the wordclouds
class SimpleGroupedColorFunc(object):
    """Create a color function object which assigns EXACT colors
    to certain words based on the color to words mapping

    Parameters
    ----------
    color_to_words : dict(str -> list(str))
        A dictionary that maps a color to the list of words.

    default_color : str
        Color that will be assigned to a word that's not a member
        of any value from color_to_words.
    """

    def __init__(self, color_to_words, default_color):
        self.word_to_color = {word: color
                              for (color, words) in color_to_words.items()
                              for word in words}

        self.default_color = default_color

    def __call__(self, word, **kwargs):
        return self.word_to_color.get(word, self.default_color)

# Import default colours based on the V gene
TRAV = json.load(open("./dict/TRAV.json"))
TRBV = json.load(open("./dict/TRBV.json"))
TRGV = json.load(open("./dict/TRGV.json"))
TRDV = json.load(open("./dict/TRDV.json"))
    

def format_data(args):
    
    reader = airr.read_rearrangement(args.rearrangements)
    empty_list=[]

    # keep only the Junction, Vgene and Repertoire ID columns
    keys = ['junction_aa', 'v_call', 'repertoire_id']

    for row in reader:
        empty_list.append({x:row[x] for x in keys})

    df = pd.DataFrame(empty_list)

    # replace cells without junction with Nan
    df['junction_aa'].replace('', np.nan, inplace=True)
    
    # delete lines with Nan
    df.dropna(subset=['junction_aa'], inplace=True)
    
    # delete lines with an X on the junction
    df = df[~df.junction_aa.str.contains("X")]
    
    # Series to select only one Vgene when there are multiple in the column
    temp_v_call  = df.v_call.str.split(",", n = 0, expand = True)
    
    # replace the v_call column with the new one
    df['v_call'] = temp_v_call[0]

    # format the df to aggregate by junction and count them
    aggregate=df.pivot_table(index=['junction_aa','v_call'], 
                                aggfunc='size').reset_index()
    
    # replace the old information in the df
    df['junction_aa']= aggregate['junction_aa']
    df['v_call']= aggregate['v_call']
    df['counts']= aggregate[0]
    
    # the size of the df changed so removing extra lines
    df.dropna(subset=['junction_aa'], inplace=True)
    
    # sort the df by highest numbers of counts
    df=df.sort_values(by='counts', ascending=False)
    
    # remove allele information from v_call and keep only the gene information
    df['v_call'] = df.apply(lambda x: x['v_call'][:-3], axis = 1)
    df['chain'] = df.apply(lambda x: x['v_call'][2], axis = 1)
        
    return df



def wordcloud(args):
    
    samples_df=format_data(args)
    samples = samples_df.groupby(['chain','repertoire_id'])
    keys = [key for key, _ in samples]
    
    for j in keys:
        df=samples.get_group(j)
        
        if len(df) > 1:
            # create a dict associating the junction to the gene
            # with highest counts
            family=df[['junction_aa','v_call']].drop_duplicates(
                subset='junction_aa', 
                keep="first", 
                inplace=False).set_index('junction_aa').squeeze().to_dict()
            # create a dict associating the junction with a count number
            text=df[['junction_aa',
                'counts']].groupby('junction_aa').sum().squeeze().to_dict()
        else:
            family={df['junction_aa'].iloc[0]:df['v_call'].iloc[0]}
            text={df['junction_aa'].iloc[0]:df['counts'].iloc[0]}

        # create the wordcloud
        
        wordcloud = WordCloud(
                        width=1000, 
                        height=1000,
                        background_color="white",
                        relative_scaling=1,
                        prefer_horizontal=1.0,
                        max_words=len(df)).generate_from_frequencies(text)
        color_to_words = {}
        for i in family:
            color_to_words.setdefault(
                eval(family.get(i)[:4]).get(family.get(i)),[]).append(i)
        default_color = 'grey'
        grouped_color_func = SimpleGroupedColorFunc(color_to_words, default_color)
        wordcloud.recolor(color_func=grouped_color_func)

        plt.figure(dpi=300.0)
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.xticks([])
        plt.yticks([])

        colours_for_legend = {}
        for i in family:
            tempdict=eval(family.get(i)[:4]).get(family.get(i))
            colours_for_legend[family.get(i)]=tempdict

        sorted_legend=sorted(colours_for_legend)
        patchList = []
        for key in sorted_legend:
            data_key = mpatches.Patch(color=colours_for_legend[key], label=key)
            patchList.append(data_key)
        outputname = args.rearrangements[:-4]+"_"+j[1]+"_"+j[0]+".png"
        plt.legend(handles=patchList,
                    bbox_to_anchor=(1.55, 1.0),
                    loc='upper right',
                    ncol=2,
                    prop={'size': 6})
        plt.savefig(outputname,dpi=300,bbox_inches='tight')