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
            fam.hasFather = True
        elif status == "Mother":
            fam.mother = newperson
            fam.hasMother = True
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
    
    for fam in families.values():
        subfamilies = [fam]
        for person in fam.people:
            if person.affected and person != fam.child:
                subfamily = Family(fam.ID)
                subfamily.child = person
                opposite = ""
                if person == fam.father or person == fam.mother:
                    opposite = fam.mother if person == fam.father else fam.father
                for p in fam.people:
                    if p != opposite:
                        subfamily.people.append(p)
                subfamilies.append(subfamily)

        famresult = pd.DataFrame()
        for subfam in subfamilies:
            famresult = pd.concat([famresult, ad_model(df, subfam)])
            famresult = pd.concat([famresult, ar_model(df, subfam)])
            famresult = pd.concat([famresult, xl_model(df, subfam)])
            famresult = pd.concat([famresult, xldn_model(df, subfam)])
            famresult = pd.concat([famresult, de_novo_model(df, subfam)])
            famresult = pd.concat([famresult, cmpd_het_model(df, subfam)])
        famresult["loc"] = [chrom+str(start)+str(end) for chrom, start,end in zip(famresult['Chr'], famresult['Start'], famresult["End"])]
        uniquelocs = famresult["loc"].unique()
        combined = pd.DataFrame()
        for loc in uniquelocs:
            rows = famresult[famresult["loc"]==loc]
            output = rows.head(1).copy()
            samplestring = ""
            for sample in rows["sample"]:
                samplestring += sample+","
            samplestring = samplestring[:-1] #remove last comma
            modelstring = ""
            for model in rows["inh model"]:
                modelstring += model+","
            modelstring = modelstring[:-1] #remove last comma
            output["sample"] = [samplestring]
            output["inh model"] = [modelstring]
            del(output["loc"])
            combined = pd.concat([combined,output])
        result = pd.concat([result, combined])
    
    
    # organize result first by inh model and then by sample
    result = result.sort_values(['sample', 'inh model'])
    result.to_csv(args.output)
    print(result)
