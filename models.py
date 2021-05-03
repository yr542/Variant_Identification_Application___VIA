# This file is for functions that apply the filters to the families based on
# the different inheritance models (ar, xl, xldn, ch, addn, ad)
# where:
# ar: autosomal recessive
# xl: x-linked
# xldn: x-linked de novo
# ch: compount heterozygous
# addn: autosomal dominant de novo
# ad: autosomal dominant

# The functions should take in the variant data frame and a Family object,
# and output a data frame of possible variants.

from filters import *
import pandas as pd

#add_columns adds three columns to the dataFrame df containing info to be outputted for
#candidate variants for a specific Family object fam and the relevant inheritance
#model number modelno
def add_columns(df, fam, modelno):
    df.insert(0, "inh model", modelno)
    df.insert(1, "family", fam.ID)
    df.insert(2, "sample", fam.child)

#ad_model takes in a data frame and Family object and returns a new data frame
#containing candidate variants
def ad_model(df, fam):
    min_allelic_depth = 1650 * 6
    numAffected = 0
    newdf = filter_AF_into_new_DataFrame(df, .0005)
    if fam.mother_phen == "Affected":
        numAffected += 1
        filter_zyg(newdf, fam.mother, "0/1")
        filter_ADs(newdf, fam.mother, min_allelic_depth)
    if fam.father_phen == "Affected":
        numAffected += 1
        filter_zyg(newdf, fam.father, "0/1")
        filter_ADs(newdf, fam.father, min_allelic_depth)
    if fam.child != "":
        numAffected += 1
        filter_zyg(newdf, fam.child, "0/1")
        filter_ADs(newdf, fam.child, min_allelic_depth)
    for sib in fam.siblings:
        if sib.phen == "Affected":
            numAffected += 1
            filter_zyg(newdf, sib.ID, "0/1")
            filter_ADs(newdf, sib.ID, min_allelic_depth)
    if numAffected <= 1:
        return pd.DataFrame() #returns an empty Data Frame if nothing should be output for this model
    else:
        add_columns(newdf, fam, 4)
        return newdf
