import argparse
import pandas as pd
from family import *
from models import *

if __name__ == '__main__':
    argp = argparse.ArgumentParser()
    argp.add_argument('-p', '--pedfile', default="Test_Ped.txt")
    argp.add_argument('-d', '--data', default="Test_cleaned.txt")

    args = argp.parse_args()

    pedDf = pd.read_csv(args.pedfile,sep='\t')
    print(pedDf)

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
#    for family in families.values():
#        print(family.father.ID)
#        print(family.mother.ID)
#        print(family.child.ID)
#        for sibling in family.siblings:
#            print(sibling.ID)

    df = pd.read_csv(args.data, sep='\t')
    
    # The following code calls all 4 models and then outputs a csv file
    # with the rows resulting for each one
    result = pd.DataFrame()
    
    for family in families.values():
        # result = pd.concat([result, MODEL_1_METHOD_NAME(df, family)])
	    # result = pd.concat([result, MODEL_2_METHOD_NAME(df, family)])
        result = pd.concat([result, de_novo_model(df, family)])
        result = pd.concat([result, ad_model(df, family)])
    
    # organize result first by inh model and then by sample
    result = result.sort_values(['inh model', 'sample'])
    result.to_csv('filtered.csv')
