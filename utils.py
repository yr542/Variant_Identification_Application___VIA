import pandas as pd
from family import *

def get_families(pedfile):
    pedDf = pd.read_csv(pedfile, sep='\t')

    families = {}
    for i in range(0, len(pedDf)):
        Fam_ID = pedDf["Family_ID"][i]
        if Fam_ID not in families.keys():
            fam = Family(Fam_ID)
            families[Fam_ID] = fam
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

    return families

import os
def load_phen(families, phenfile, mapfile):
    if not os.path.isfile(mapfile):
        answer = input("No phenotype-to-gene mapping found. Download one? [y/N]: ").lower()

        if answer=="y":
            import urllib.request
            url = "https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/phenotype_to_genes.txt"
            print("Downloading now. It might take a while...")
            urllib.request.urlretrieve(url,mapfile)
            print("Finished downloading.")
        else:
            print("exiting now")
            exit()
    phen_to_genes = pd.read_csv(mapfile, sep = '\t',
            header = None, comment = '#')

    phen_to_genes.columns = ["HPO-id", "HPO label", "gene-id",
            "gene-symbol", "additional-info",
            "source", "disease-ID"]
    #loads HPO numbers
    phenDf = pd.read_csv(phenfile, sep='\t')
    for j in families:
        for i in range(0, len(phenDf)):
            if phenDf["Family_ID"][i] == j:
                hpo = phenDf["HPO"][i]
                fam = families.get(j)
                fam.HPO = hpo.split(',')
                fam.genes = {}
                for HPO in fam.HPO:
                    genes = phen_to_genes[phen_to_genes["HPO-id"]==HPO]['gene-symbol'].tolist()
                    gene_nums = [fam.genes[gene]+1 if gene in fam.genes else 1 for gene in genes]
                    fam.genes.update(dict(zip(genes, gene_nums)))

def generate_subfamilies(fam):
    subfamilies = [fam]
    for person in fam.people:
        if person.affected and person != fam.child:
            subfamily = Family(fam.ID)
            subfamily.child = person
            opposite = ""
            if person == fam.father or person == fam.mother:
                opposite = fam.mother if person == fam.father else fam.father
            if person in fam.siblings:
                if fam.hasFather:
                    subfamily.father = fam.father
                    subfamily.hasFather = True
                if fam.hasMother:
                    subfamily.mother = fam.mother
                    subfamily.hasMother = True
                newsibs = fam.siblings.copy()
                newsibs.remove(person)
                subfamily.siblings = newsibs + [fam.child]
            for p in fam.people:
                if p != opposite:
                    subfamily.people.append(p)
            subfamilies.append(subfamily)
    return subfamilies

# combine multiple instances of the same variant into one row.
def combine_duplicates(df):
    # use a generated location string to determine if two variants
    # are the same or not.
    df["loc"] = [chrom + str(start) + str(end) for chrom, start, end in
                        zip(df['Chr'], df['Start'], df["End"])]
    uniquelocs = df["loc"].unique()
    combined = pd.DataFrame()
    for loc in uniquelocs:
        rows = df[df["loc"] == loc]
        output = rows.head(1).copy()
        samplestring = ""
        for sample in rows["sample"]:
            samplestring += sample + ","
        samplestring = samplestring[:-1]  # remove last comma
        modelstring = ""
        for model in rows["inh model"]:
            modelstring += model + ","
        modelstring = modelstring[:-1]  # remove last comma
        output["sample"] = [samplestring]
        output["inh model"] = [modelstring]
        del (output["loc"])
        combined = pd.concat([combined, output])
    return combined
