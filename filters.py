# filter the dataFrame (df) by the minimum allele depth (ad) in a particular
# column (name)
def filter_ADs(df, name, ad):
    strings=np.array(df[name])
    ADindices=[s.split(":").index("AD") for s in df["FORMAT"]]
    ADs=[strings[i].split(":")[ADindices[i]] for i in range(len(strings))]
    ADs=[int(ad.split(",")[1])/max(int(ad.split(",")[0]),1) for ad in ADs]
    df["AD"]=ADs
    df=df[df["AD"]>ad]
    print(len(df))
    return df

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
    print(len(df))
    return df

# filter the dataFrame (df) by the maximum population allele frequency (cap)
def filter_AF(df, cap):
    AF_columns=["AF_popmax","PopFreqMax","GME_AF","Kaviar_AF","abraom_freq"]
    AF_columns = AF_columns + [col + ".1" for col in AF_columns]

    for col in AF_columns:
            if col in df.columns:
                    df.loc[df[col]==".",col]="-1"
                    df[col]=df[col].astype(float)
                    df=df[(df[col]<=cap)]
    print(len(df))
    return df

# filter the dataFrame (df) for the zygosity (zyg), e.g. "0/1", in a particular
# column (name)
def filter_zyg(df, name, zyg):
    if name in df.columns:
        df=df[df[name].str.contains(zyg)]
    return df

# filter out variants that are "Benign" or "Likely benign"
def filter_benign(df):
    df=df[(df["CLNSIG"].str.contains("enign")==False)]
    return df

# filters the dataFrame (df) by the maximum population allele frequency (cap),
# but puts the candidate variants in a new data frame without changing the old one
# and returns the new data frame (efficiency measure to maintain old data frame for
# multiple models)
def filter_AF_into_new_DataFrame(df, cap):
    newdf = pd.DataFrame()
    firstcol = 1
    AF_columns = ["AF_popmax", "PopFreqMax", "GME_AF", "Kaviar_AF", "abraom_freq"]
    AF_columns = AF_columns + [col + ".1" for col in AF_columns]

    for col in AF_columns:
        if firstcol == 1 and col in df.columns:
            df.loc[df[col] == ".", col] = "-1"
            df[col] = df[col].astype(float)
            newdf = df[(df[col] <= cap)]
            firstcol = 0
        elif col in df.columns:
            newdf.loc[newdf[col] == ".", col] = "-1"
            newdf[col] = newdf[col].astype(float)
            newdf = newdf[(newdf[col] <= cap)]

    return newdf
