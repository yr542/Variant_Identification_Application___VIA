import argparse
import pandas as pd
import family

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
            fam = family.Family(Fam_ID)
            families[Fam_ID]=fam
        else:
            fam = families.get(Fam_ID)

        ID = pedDf["individual_ID"][i]
        status = pedDf["Status"][i]
        sex = pedDf["Sex"][i]
        phen = pedDf["Phenotype"][i]
        if status == "Father":
            fam.father=ID
            fam.father_phen=phen
        elif status == "Mother":
            fam.mother=ID
            fam.mother_phen=phen
        elif status == "Child":
            fam.child=ID
        elif status == "Sibling":
            sibling = family.Sibling(ID, sex, phen)
            fam.siblings.append(sibling)
#    for family in families.values():
#        print(family.father)
#        print(family.mother)
#        print(family.child)
#        for sibling in family.siblings:
#            print(sibling.ID)

    df = pd.read_csv(args.data, sep='\t')
#   for family in families.values():
#         ...
