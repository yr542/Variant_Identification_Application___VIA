import pandas as pd
import numpy as np


# filter the dataFrame (df) by the minimum allele depth (ad) in a particular
# column (name)
def filter_AD(df, name, ad):
    strings=np.array(df[name])
    ADindices=[s.split(":").index("AD") for s in df["FORMAT"]]
    ADs=[strings[i].split(":")[ADindices[i]] for i in range(len(strings))]
    ADs=[int(ad.split(",")[1])/max(int(ad.split(",")[0]),1) for ad in ADs]
    df["AD"]=ADs
    df=df[df["AD"]>ad]
    #print(len(df))
    return df

# filter the dataFrame (df) by minimum depth in a particular column (name)
# if inplace is set to any integer other than 1, it will be filtered into a new data frame
# by default, the function filters in place when the inplace arg is left out of the function call
def filter_DP(df, name, dp, inplace=1):
    strings = np.array(df[name])
    DPindices=[s.split(":").index("DP") for s in df["FORMAT"]]
    DPs=[int(strings[i].split(":")[DPindices[i]]) for i in range(len(strings))]
    df["DP"]=DPs
    if inplace == 1:
        df=df[df["DP"] >= dp].copy()
        #print(len(df))
        return df
    else:
        dfcopy = df[df["DP"] >= dp].copy()
        return dfcopy

# filter the dataFrame (df) by the maximum number of occurences (cap) of a
# particular zygosity (zyg), e.g. "0/1", in a range of columns
# [namestart,nameend]
def filter_occurences(df, zyg, namestart, nameend, cap):
    freqs=[]
    for i in range(df.columns.get_loc(namestart), df.columns.get_loc(nameend)+1):
        mask=df.iloc[:,i].str.contains(zyg)
        freqs.append(mask)
    freqs=sum(freqs)
    indices=[]
    for key in freqs.keys():
        if freqs[key]>cap:
            indices.append(key)
    df.drop(indices,inplace=True)
    #print(len(df))
    return df

# filter the dataFrame (df) by the maximum population allele frequency (cap)
def filter_AF(df, cap):
    AF_columns=["AF_popmax","PopFreqMax","GME_AF","Kaviar_AF","abraom_freq"]
    AF_columns = AF_columns + [col + ".1" for col in AF_columns]

    for col in AF_columns:
        if col in df.columns:
            df.loc[df[col]==".",col]="-1"
            df[col]=df[col].astype(float)
            df=df[df[col]<=cap].copy()
    #print(len(df))
    return df

# filter the dataFrame (df) for the zygosity (zyg), e.g. "0/1", in a particular
# column (name)
def filter_zyg(df, name, zyg):
    if name in df.columns:
        df=df[df[name].str.contains(zyg)]
    return df

# filter the dataFrame (df) to exclude a certain zygosity (zyg) in a particular
# column (name)
def exclude_zyg(df, name, zyg):
    if name in df.columns:
        df = df[~df[name].str.contains(zyg)]
    return df

# filter out variants that are "Benign" or "Likely benign"
def filter_benign(df):
    df=df[(df["CLNSIG"].str.contains("enign")==False)]
    return df

# filter the dataFrame (df) by MAXIMUM depth in a particular column (name)
#if inplace is 1, it filters df in place; if option is not 1, it filters into a new data frame
def filter_DP_Max(df, name, dp, inplace=1):
    strings = np.array(df[name])
    DPindices = [s.split(":").index("DP") for s in df["FORMAT"]]
    DPs = [int(strings[i].split(":")[DPindices[i]]) for i in range(len(strings))]
    df["DP"] = DPs
    if inplace == 1:
        df = df[df["DP"] < dp]
        return df
    else:
        dfcopy = df[df["DP"] < dp].copy()
        return dfcopy

# filter the dataFramd (df) if you only want to keep the rows in which the gene is
# located in the X chromosome
def filter_chr(df, name, chrom):
    if name in df.columns:
        df=df[df[name].str.contains(chrom)]
    return df
