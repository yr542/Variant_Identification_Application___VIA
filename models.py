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
    df.insert(2, "sample", fam.child.ID)

#ad_model takes in a data frame and Family object and returns a new data frame
#containing candidate variants
def ad_model(df, fam):
    min_allelic_depth = 6
    numAffected = 0
    newdf = filter_AF_into_new_DataFrame(df, .0005)
    for person in fam.people:
        if person.phen == "Affected":
            numAffected += 1
            newdf = filter_zyg(newdf, person.ID, "0/1")
            newdf = filter_DP(newdf, person.ID, min_allelic_depth)
        else:
            newdf = filter_zyg(newdf, person.ID, "0/0")
    if numAffected <= 1:
        return pd.DataFrame()  # returns an empty Data Frame if nothing should be output for this model
    else:
        add_columns(newdf, fam, 4)
        return newdf

# de_novo_model takes a dataframe (the cleaned data) and a family object
# return value: a new dataframe with all possible de novo candidate
# genes
def de_novo_model(df, fam):

    # re-filter for MAF
    revised_df = filter_AF_into_new_DataFrame(df, .0005)
    
    # keep track of number of individuals we are identifying variants
    # for
    num_affected = 0

    # If either mother or father is affected, no de novo, so return
    # empty data frame
    if fam.mother.phen == "Affected" or fam.father.phen == "Affected":
        return pd.DataFrame()
   
    # filter child for all 0/1
    if fam.child.ID != "":
        num_affected += 1
        revised_df = filter_zyg(revised_df, fam.child.ID, "0/1")
        
	# filter to make sure DP is at least 6x
        revised_df = filter_DP(revised_df, fam.child.ID, 6)

        # check that no unaffected siblings are 0/1
        de_novo_check_siblings(revised_df, fam)

    # filter parents for 0/0
    revised_df = filter_zyg(revised_df, fam.father.ID, "0/0")
    revised_df = filter_zyg(revised_df, fam.mother.ID, "0/0")

    # filter siblings to identify more candidate genes
    for sib in fam.siblings:
        if sib.phen == "Affected":
            num_affected += 1
            revised_df = filter_zyg(revised_df, fam.sib.ID, "0/1")
            revised_df = filter_DP(revised_df, fam.sib.ID, 6)
            de_novo_check_siblings(revised_df, fam)
            revised_df = pd.concat([revised_df, sib_df])
    
    if num_affected:

        # add on the columns with family info
        add_columns(revised_df, fam, 3)
        return revised_df

    # if no affected individuals, return empty data frame
    return pd.DataFrame()

# Checks to make sure that no sibling is unaffected and also 0/1
def de_novo_check_siblings(df, fam):

    # keep a list of rows to remove
    to_remove = []
    for i, row in df.iterrows():
        errors = 0
        for sib in fam.siblings:
            
	    # add to list to remove if sibling is Unaffected yet also
	    # 0/1
            if sib.phen != "Affected" and row[sib.ID] == "0/1":
                errors += 1
            if errors:
                to_remove.append(i)
    
    # remove all bad rows
    df.drop(index = to_remove)


def cmpd_het_model(df, fam):
   
    # keep track of individuals we are identifying variants for
    num_affected = 0

    # filter child for having 0/1 in 2 variants of the same gene

    # first create newdf to include all instances of child 0/1
    if fam.child != "":
        num_affected += 1
        newdf = filter_zyg(df, fam.child, "0/1")

        # use Gene.refGene column to create new column in df with
        # gene names (deals with semicolon issue)
        
        partitioned_gene=0
        newdf['Gene'] = ''
        for i in range(0, len(newdf)):
            partitioned_gene=newdf["Gene.refGene"][i].partition(";")
            newdf['Gene'][i]= partitioned_gene[0]
        
        # create new df with only duplicate genes for child with 0/1
        # and delete the gene column we created
        finaldf = newdf[newdf.duplicated(subset=['Gene'], keep=False)]
        del finaldf['Gene']
    
    # add on the columns with family info
    if num_affected:
        add_columns(finaldf, fam, 2)
        return finaldf

    # if no affected individuals, return empty data frame
    return pd.DataFrame()
