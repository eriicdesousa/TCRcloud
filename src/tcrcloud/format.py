import pandas as pd

import airr


def format_data(args):

    with open(args.rearrangements) as f:
        first_line = f.readline()
        if "duplicate_count" in first_line:
            keys = ["junction_aa", "v_call", "j_call", "junction",
                    "repertoire_id", "duplicate_count", "productive"]
        else:
            keys = ["junction_aa", "v_call", "j_call", "junction",
                    "repertoire_id", "productive"]

    airr.validate_rearrangement(args.rearrangements, True)
    reader = airr.read_rearrangement(args.rearrangements)
    print("Preparing data")
    empty_list = []
    # keep only part of the data
    for row in reader:
        productive = row.get("productive")
        CDR3 = row.get("junction_aa")
        try:
            v_call = row.get("v_call")[2]
        except IndexError:
            v_call = ""
        try:
            j_call = row.get("j_call")[2]
        except IndexError:
            j_call = ""
        if productive is True:
            if CDR3 != "":
                if v_call != "":
                    if j_call != "":
                        if "X" or "\*" not in CDR3:
                            if CDR3[0] == "C":
                                if CDR3[-1] == "F" or CDR3[-1] == "W":
                                    if v_call == j_call:
                                        empty_list.append({x: row[x] for x in keys})

    df = pd.DataFrame(empty_list)

    # keep only one first V gene when there are multiple in the column
    df["v_call"] = df.v_call.str.split(",", n=1, expand=True)[0]

    # remove allele information from v_call and keep only the gene information
    if "*" in str(df.v_call.head(1)):
        df["v_call"] = df.apply(lambda x: x["v_call"][:-3], axis=1)

    # create column with chain information
    df["chain"] = df.apply(lambda x: x["v_call"][2], axis=1)

    return df


def format_convergence(df):
    # format the df to aggregate by junction
    aggregate = df.pivot_table(index=["junction_aa", "junction",
                                      "repertoire_id", "chain"], aggfunc="size"
                               ).reset_index()
    del df
    for_convergence = aggregate.pivot_table(index=["junction_aa",
                                                   "repertoire_id", "chain"],
                                            aggfunc="size").reset_index()

    for_convergence.rename(columns={0: "counts"}, inplace=True)
    for_convergence = for_convergence.sort_values(by="counts", ascending=False)
    return for_convergence


def format_metrics(df):
    if "duplicate_count" in df.columns:
        aggregate = df.loc[:, ("junction_aa", "repertoire_id",
                               "chain", "duplicate_count")]
        del df
        aggregate.rename(columns={"duplicate_count": "counts"}, inplace=True)
        aggregate = aggregate.groupby(["junction_aa",
                                       "repertoire_id", "chain"]
                                      ).sum().reset_index()
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    else:
        aggregate = df.pivot_table(index=["junction_aa", "repertoire_id",
                                          "chain"],
                                   aggfunc="size").reset_index()
        del df
        aggregate.rename(columns={0: "counts"}, inplace=True)
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    return aggregate


def format_cloud(df):
    if "duplicate_count" in df.columns:
        aggregate = df.loc[:, ("junction_aa", "v_call", "repertoire_id",
                               "chain", "duplicate_count")]
        del df
        aggregate.rename(columns={"duplicate_count": "counts"}, inplace=True)
        aggregate = aggregate.groupby(["junction_aa", "v_call",
                                       "repertoire_id", "chain"]
                                      ).sum().reset_index()
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    else:
        aggregate = df.pivot_table(index=["junction_aa", "v_call",
                                          "repertoire_id", "chain"],
                                   aggfunc="size").reset_index()
        del df
        aggregate.rename(columns={0: "counts"}, inplace=True)
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    return aggregate


def format_barplot(df):
    if "duplicate_count" in df.columns:
        aggregate = df.loc[:, ("junction_aa", "v_call", "repertoire_id",
                               "chain", "duplicate_count")]
        del df
        aggregate.rename(columns={"duplicate_count": "counts"}, inplace=True)
        aggregate = aggregate.groupby(["junction_aa", "v_call",
                                       "repertoire_id", "chain"]
                                      ).sum().reset_index()
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    else:
        aggregate = df.pivot_table(index=["junction_aa", "v_call",
                                          "repertoire_id", "chain"],
                                   aggfunc="size").reset_index()
        del df
        aggregate.rename(columns={0: "counts"}, inplace=True)
        aggregate = aggregate.sort_values(by="counts", ascending=False)
    aggregate["CDR3_length"] = aggregate["junction_aa"].str.len()
    return aggregate


# print("""
#                   _____ ____ ____      _                 _
#                  |_   _/ ___|  _ \ ___| | ___  _   _  __| |
#                    | || |   | |_) / __| |/ _ \| | | |/ _` |
#                    | || |___|  _ < (__| | (_) | |_| | (_| |
#                    |_| \____|_| \_\___|_|\___/ \__,_|\__,_|
#    """)
