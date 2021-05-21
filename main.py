import argparse
import pandas as pd
from family import *
from models import *

if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    argp.add_argument('-p', '--pedfile', default="Test_Ped.txt")
    argp.add_argument('-d', '--data', default="Test_cleaned.txt")
    argp.add_argument('-o', '--output', default="filtered.csv")
    argp.add_argument('-f', '--family', default="")

    args = argp.parse_args()

    pedDf = pd.read_csv(args.pedfile,sep='\t')

    families = {}
    for i in range(0, len(pedDf)):
        #print(i)
        Fam_ID = pedDf["Family_ID"][i]
        if Fam_ID not in families.keys():
            fam = Family(Fam_ID)
            families[Fam_ID]=fam
        else:
            fam = families.get(Fam_ID)

        ID = pedDf["individual_ID"][i]
        status = pedDf["Status"][i]
        sex = pedDf["Sex"][i]
        phen = pedDf["Phenotype"][i]
        newperson = Person(ID, sex, phen)
        fam.people.append(newperson)
        if status == "Father":
            fam.father = newperson
        elif status == "Mother":
            fam.mother = newperson
        elif status == "Child":
            fam.child = newperson
        elif status == "Sibling":
            fam.siblings.append(newperson)


    df = pd.read_csv(args.data, sep='\t')
    
    #csv with variants in one family
    if args.family != "":
        fam = families[args.family]
        fam_variants = df.copy()
        for person in fam.people:
            filt = filter_zyg if person.phen == "Unaffected" else exclude_zyg
            fam_variants = filt(fam_variants, person.ID, "0/0")
        fam_variants.to_csv(fam.ID+".csv")

    # The following code calls all 4 models and then outputs a csv file
    # with the rows resulting for each one
    result = pd.DataFrame()
    
    for family in families.values():
        # result = pd.concat([result, MODEL_1_METHOD_NAME(df, family)])
        result = pd.concat([result, cmpd_het_model(df, family)])
        result = pd.concat([result, de_novo_model(df, family)])
        result = pd.concat([result, ad_model(df, family)])
    
    # organize result first by inh model and then by sample
    result = result.sort_values(['inh model', 'sample'])
    result.to_csv(args.output)
