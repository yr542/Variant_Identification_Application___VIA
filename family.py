# Family class. Stores the strings and phenotype for every member of the family.

class Family:
    def __init__(self, ID):
        self.ID = ID
        self.siblings = []
        self.mother=""
        self.mother_phen=""
        self.father=""
        self.father_phen=""
        self.child=""

# Sibling class. Stores the string, sex, and phenotype for a sibling.
class Sibling:
    def __init__(self, ID, sex, phen):
        self.ID = ID
        self.sex = sex
        self.phen = phen
