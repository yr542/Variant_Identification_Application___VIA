# This file is for functions that apply the filters to the families based on
# the different inheritance models (ar, xl, xldn, ch, addn, ad)
# where:
# ar: autosomal recessive (model #1)
# xl: x-linked (model #1)
# xldn: x-linked de novo (model #1)
# ch: compound heterozygous (model #2)
# addn: autosomal dominant de novo (model #3)
# ad: autosomal dominant (model #4)

# The functions should take in the variant data frame and a Family object,
# and output a data frame of possible variants.

from family import Family
from filters import *
import pandas as pd

#add_columns adds three columns to the dataFrame df containing info to be outputted for
#candidate variants for a specific Family object fam and the relevant inheritance
#model number modelno
def add_columns(df, fam, modelno):
    df.insert(0, "inh model", modelno)
    df.insert(1, "family", fam.ID)
    df.insert(2, "sample", fam.child.ID)

# ad_model takes in a data frame and Family object and returns a new data frame
# containing candidate variants
def ad_model(df, fam):
    min_allelic_depth = 6  # will filter for 6x coverage minimum for at least one affected individ
    numAffected = 0
    newdf = df.copy()
    newdf = filter_AF(newdf, .0005)  # filters all AF cols for entries <= .0005
    dpdf = pd.DataFrame()
    for person in fam.people:
        if person.phen == "Affected":
            numAffected += 1
            newdf = filter_zyg(newdf, person.ID, "0/1")  # filters for 0/1 entries for affected individs
            if numAffected == 1: #entire if-else works to filter DP appropriately
                dpdf = filter_DP_Max(newdf, person.ID, min_allelic_depth, 0) #saves low DP variants for future use
                newdf = filter_DP(newdf, person.ID, min_allelic_depth) #filters for min DP for first affected individ
            else: #makes sure if later affected individs have high enough DPs those variants are included in the output
                  #even if the first affected individ had a low DP
                dpdf2 = filter_DP(dpdf, person.ID, min_allelic_depth, 0)
                newdf.append(dpdf2, ignore_index=True)
                dpdf = filter_DP_Max(dpdf, person.ID, min_allelic_depth) #removes kept variants from dpdf
        else:
            newdf = filter_zyg(newdf, person.ID, "0/0")  # filters for 0/0 entries for unaffected individs

    # returns an empty Data Frame if nothing should be output for this model (<= 1 affected individs)
    if numAffected <= 1:
        return pd.DataFrame()
    else:
        add_columns(newdf, fam, 4)  # adds on columns with family info
        return newdf

# de_novo_model takes a dataframe (the cleaned data) and a family object
# return value: a new dataframe with all possible de novo candidate
# genes
def de_novo_model(df, fam):

    # re-filter for MAF
    revised_df = df.copy()
    revised_df = filter_AF(revised_df, .0005)
    
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

    # filter child for having 0/1 in >=2 variants of the same gene
	
    # create newdf to include all instances of child 0/1
    if fam.child.ID != "":
        num_affected += 1
        newdf = df.copy()
        newdf = filter_zyg(newdf, fam.child.ID, "0/1")
	
        # use Gene.refGene column to create new column "Gene' with
        # gene names (deals with semicolon issue in some genes)
        partitioned_gene=0
        newdf.loc[:,'Gene'] = ''
        newdf=newdf.dropna(subset=["Gene.refGene"])
        newdf = newdf.reset_index(drop=True)
        newdf['Gene'] = newdf["Gene.refGene"].copy().str.partition(";")[0]
        
        # create new df where genes must meet the following criteria: 
	# there must be at least 1 0/1 variant in mother that is not in father
	# and at least 1 0/1 variant in father that is not in mother
        finaldf = newdf[newdf.duplicated(subset=['Gene'], keep=False)] 

        if(fam.father.ID != "" and fam.mother.ID != ""):
            genes = finaldf["Gene"].unique()
            for gene in genes:
                genedf = finaldf[finaldf["Gene"]==gene]
                mom = sum(genedf[fam.mother.ID].str.contains("0/1") &
                          genedf[fam.father.ID].str.contains("0/0"))
                dad = sum(genedf[fam.mother.ID].str.contains("0/0") &
                          genedf[fam.father.ID].str.contains("0/1"))
                if(mom==0 or dad==0):
                    finaldf = finaldf[finaldf["Gene"]!=gene]
        # delete the gene column we created
        del finaldf['Gene']

        
    # add on the columns with family info
    if num_affected:
        add_columns(finaldf, fam, 2)
        return finaldf





def xl_model(df, fam):
    newdf = df.copy()
    min_allelic_depth = 0.5
    numAffected = 0
    name = 'Chr'
    x_df = filter_chr(newdf, name, "chrX")
    if fam.mother.phen == "Affected" or fam.father.phen == "Affected":
        return pd.DataFrame()
    # filter child for all 0/1
    if fam.child.ID != "":
        if fam.child.sex == 'Male':
            if fam.mother.ID == '0/1':
                numAffected += 1
                x_df = filter_zyg(x_df, fam.child.ID, "1/1")
    add_columns(x_df, fam, 1)
    return(x_df)

def xldn_model(df, fam):
    newdf = df.copy()
    min_allelic_depth = 0.5
    numAffected = 0
    name = 'Chr'
    x_df = filter_chr(newdf, name, "chrX")
    if fam.mother.phen == "Affected" or fam.father.phen == "Affected":
        return pd.DataFrame()
    # filter child for all 0/1
    if fam.child.ID != "":
        if fam.child.sex == 'Male':
            if fam.mother.ID == '0/0':
                numAffected += 1
                x_df = filter_zyg(x_df, fam.child.ID, "1/1")
    add_columns(x_df, fam, 1)
    return(x_df)
    

# ar_model takes the data frame and Family object. Returns: a new data frame containing
# all possible autosomal recessive candidate genes
def ar_model(df, fam):
    min_allelic_depth = 6  # will filter for 6x coverage minimum for at least one affected individ
    numAffected = 0
    newdf = df.copy()
    newdf = filter_AF(newdf, .005) 
    dpdf = pd.DataFrame()
    for person in fam.people:
        if person.phen == "Affected":
            numAffected += 1
            newdf = filter_zyg(newdf, person.ID, "1/1") 
            if numAffected == 1: 
                dpdf = filter_DP_Max(newdf, person.ID, min_allelic_depth, 0) #saves low DP variants for future use
                newdf = filter_DP(newdf, person.ID, min_allelic_depth)
    # returns an empty Data Frame if nothing should be output for this model (<= 1 affected individs)
    if numAffected < 1:
        return pd.DataFrame()
    else:
        add_columns(newdf, fam, 1)  # adds on columns with family info
        return newdf
    
