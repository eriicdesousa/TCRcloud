import json
import sys

import matplotlib

from wordcloud import (WordCloud)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import tcrcloud.colours
import tcrcloud.format


matplotlib.use('Agg')

# Import default colours based on the V gene
TRAV = tcrcloud.colours.TRAV
TRBV = tcrcloud.colours.TRBV
TRGV = tcrcloud.colours.TRGV
TRDV = tcrcloud.colours.TRDV


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


def wordcloud(args):

    samples_df = tcrcloud.format.format_data(args)
    formatted_samples = tcrcloud.format.format_cloud(samples_df)
    samples = formatted_samples.groupby(['chain', 'repertoire_id'])
    keys = [key for key, _ in samples]

    for j in keys:
        df = samples.get_group(j)

        if len(df) > 1:
            # create a dict associating the junction to the gene
            # with highest counts
            family = df[['junction_aa', 'v_call']].drop_duplicates(
                subset='junction_aa',
                keep="first",
                inplace=False).set_index('junction_aa').squeeze().to_dict()
            # create a dict associating the junction with a count number
            text = df[['junction_aa',
                       'counts']].groupby('junction_aa'
                                          ).sum().squeeze().to_dict()
        else:
            family = {df['junction_aa'].iloc[0]: df['v_call'].iloc[0]}
            text = {df['junction_aa'].iloc[0]: df['counts'].iloc[0]}

        # create the wordcloud

        wordcloud = WordCloud(width=1000,
                              height=1000,
                              background_color="white",
                              relative_scaling=0.7,
                              prefer_horizontal=1.0,
                              max_words=len(df)
                              ).generate_from_frequencies(text)
        color_to_words = {}
        if args.colours is not None:
            try:
                with open(args.colours) as json_file:
                    color_to_words = json.load(json_file)
            except FileNotFoundError:
                sys.stderr.write("TCRcloud error: " + args.colours
                                 + " doesn't seem to exist\n")
                exit()
            except json.decoder.JSONDecodeError:
                sys.stderr.write("TCRcloud error: " + args.colours
                                 + " doesn't seem properly formatted. Check \
https://github.com/oldguyeric/TCRcloud for more information\n")
                exit()
        else:
            for i in family:
                color_to_words.setdefault(
                    eval(family.get(i)[:4]).get(family.get(i)), []).append(i)

        default_color = 'grey'
        try:
            grouped_color_func = SimpleGroupedColorFunc(color_to_words,
                                                        default_color)
        except TypeError:
            sys.stderr.write("TCRcloud error: " + args.colours
                             + " doesn't seem properly formatted. Check \
https://github.com/oldguyeric/TCRcloud for more information\n")
            exit()
        wordcloud.recolor(color_func=grouped_color_func)

        plt.figure(dpi=300.0)
        plt.imshow(wordcloud, interpolation="bilinear")
        plt.xticks([])
        plt.yticks([])
        colours_for_legend = {}
        if args.colours is None:
            for i in family:
                tempdict = eval(family.get(i)[:4]).get(family.get(i))
                colours_for_legend[family.get(i)] = tempdict

            sorted_legend = sorted(colours_for_legend)
            patchList = []
            for key in sorted_legend:
                data_key = mpatches.Patch(color=colours_for_legend[key],
                                          label=key)
                patchList.append(data_key)

            plt.legend(handles=patchList,
                       bbox_to_anchor=(0.5, -0.01),
                       loc='upper center',
                       ncol=4,
                       prop={'size': 6})
        outputname = args.rearrangements[:-4] + "_" + j[1] + "_" \
                                              + j[0] + ".png"
        plt.tight_layout()
        plt.savefig(outputname, dpi=300, bbox_inches='tight')
        print("Word cloud saved as " + outputname)
