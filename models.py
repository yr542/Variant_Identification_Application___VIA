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

import filters.py

#add_columns adds three columns to the dataFrame df containing info to be outputted for
#candidate variants for a specific Family object fam and the relevant inheritance
#model number modelno
def add_columns(df, fam, modelno):
    df.insert(0, "inh model", modelno)
    df.insert(1, "family", fam.ID)
    df.insert(2, "sample", 'sampleno') #sampleno will be replaced with value(s) when location is known

#ad_model takes in a data frame and Family object and returns a new data frame
#containing candidate variants
def ad_model(df, fam):
    numAffected = 0
    newdf = pd.DataFrame.copy()  #deep copy of data frame so other models can use the same one (likely inefficient, stopgap for now)
    filter_AF(newdf, .0005)
    if fam.mother_phen == "Affected":
        numAffected += 1
        filter_zyg(newdf, fam.mother, "0/1")
        filter_ADs(newdf, fam.mother, 'min allelic depth') #min allelic depth will be replaced with a number when known here and below
    if fam.father_phen == "Affected":
        numAffected += 1
        filter_zyg(newdf, fam.father, "0/1")
        filter_ADs(newdf, fam.father, 'min allelic depth')
    if fam.child != "":
        numAffected += 1
        filter_zyg(newdf, fam.child, "0/1")
        filter_ADs(newdf, fam.child, 'min allelic depth')
    for sib in fam.siblings:
        if sib.phen == "Affected"
            numAffected += 1
            filter_zyg(newdf, sib.ID, "0/1")
            filter_ADs(newdf, sib.ID, 'min allelic depth')
    if numAffected <= 1:
        return pd.DataFrame() #returns an empty Data Frame if nothing should be output for this model
    else:
        add_columns(newdf, fam, 4)
        return newdf
  #to ask in meeting: should it filter out variants where a person is 0/1 and not affected? Seems like yes but not in specs
