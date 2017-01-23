#parse PSM results to get peptide and protein identifications


#this is written in Python 3 

HelpOutput = """
python IdentifyingPeptidesAndProteins.py -d /path/to/files -o file.txt

Required Parameters:
-d   <dir>  This is a path to the directory that holds all peptide identifications
-o   <file> This is a file that tells which PSM files belong to a specific organism

Example: python .\IdentifyingPeptidesAndProteins.py -d ..\MSGF.Searches -o FileToOrganisms.txt
"""


import os
import sys
import getopt
import string
import gzip



class Protein:
    def __init__(self, accession):
        self.ProteinAccession = accession
        self.KeggOrthologAccession = ""
        self.Peptides = {} #key = sequence, value = spectrum count

    def AddPeptide(self, Peptide):
        if not Peptide in self.Peptides:
            self.Peptides[Peptide] = 0
        self.Peptides[Peptide] += 1

    def IsModified(self, PTMChar):
        #this function simply looks for whether any peptide in this is modified
        for Peptide in self.Peptides.keys():
            if PTMChar in Peptide:
                return 1
        return 0 # we got here without returning, so that means no modifications

    def IsModifiedResidueTerminal(self, PTMChar):
        #this function looks for the PTM residue and then checks to see whether
        #that residue is the terminal residue.
        #use case: acetylated lysines should not be cleaved by trypsin.
        #--------- so we look at the acetyl lysines and then assess
        #--------- whether they are terminal or not.
        CountTerminal = 0
        TotalModifiedPeptides = 0
        for Peptide in self.Peptides.keys():
            if PTMChar in Peptide:
                TotalModifiedPeptides += 1
                #now look for the index(es) of the PTMChar and return letter at index-1
                AllPTMIndicies = [i for i, ltr in enumerate(Peptide) if ltr == PTMChar] #not sure how this black magic works
                #now we just see whether the index is the terminal index or not
                LastIndex = len(Peptide) -1
                if LastIndex in AllPTMIndicies:
                    CountTerminal += 1
        return (CountTerminal, TotalModifiedPeptides)
        

    def GetModifiedResidue(self, PTMChar):
        ToReturn = [] # list of modified residues
        for Peptide in self.Peptides.keys():
            if PTMChar in Peptide:
                #now look for the index(es) of the PTMChar and return letter at index-1
                AllPTMIndicies = [i for i, ltr in enumerate(Peptide) if ltr == PTMChar] #not sure how this black magic works
                for i in AllPTMIndicies:
                    ToReturn.append(Peptide[i-1])
        return ToReturn


class Organism:
    def __init__(self, Name):
        self.OrganismName = Name
        self.Proteins = {} #key = Accession, value= object


    def AddPeptide(self, Peptide, ProteinAcc):
        if not ProteinAcc in self.Proteins:
            self.Proteins[ProteinAcc] = Protein(ProteinAcc)
        self.Proteins[ProteinAcc].AddPeptide(Peptide)
    def GetProteinCount(self):
        #simply return the number of proteins
        return len(self.Proteins)

    def GetModifiedProteinList(self, PTMChar):
        #this function retuns a list of ProteinAccessions for all the 
        #proteins that are modified
        #PTMChar is the character that is used to denote PTMs in the text files that we are parsing
        ToReturn = []
        for Accession in self.Proteins.keys():
            ProteinObject = self.Proteins[Accession]
            if ProteinObject.IsModified(PTMChar):
                ToReturn.append(Accession)

        return ToReturn

    def GetNotModifiedProteins(self, PTMChar):
        ToReturn = []
        for Accession in self.Proteins.keys():
            ProteinObject = self.Proteins[Accession]
            if not ProteinObject.IsModified(PTMChar):
                ToReturn.append(Accession)

        return ToReturn

    def GetDecoyProteins(self):
        ### the decon string is always "XXX" prepended to the protein names
        #but I put it as an optinonal parameter. 
        DecoyCount = 0
        for Accession in self.Proteins.keys():
            #print (Accession)
            if Accession[:3] == "XXX":
                DecoyCount += 1
        return DecoyCount



class ParserClass:

    def __init__(self):
        "nothing to put in really"
        self.DirectoryOfPSMs = "" # this dir holds all the PSM results (txt files)
        self.FileToOrganismsPath  = ""
        self.QvalueCutoff = 0.001
        self.OrganismObjectsDictionary = {} #key = Name, value= object
        self.FileToOrganismDictionary = {} #key = PSM file stub, value = organism


    def GetOrganismObjects(self):
        return self.OrganismObjectsDictionary

    def SetQvalue(self, qvalue):
        #for the ipython notebook
        self.QvalueCutoff = qvalue

    def Main(self):
        
        #0. Figure out which PSM files belong with which organisms
        self.ParseFileToOrganism(self.FileToOrganismsPath)
        #1. parse all the PSM files just putting peptides into proteins
        ItemsInDir = os.listdir(self.DirectoryOfPSMs)
        for Item in ItemsInDir:
            #I'm expecting .txt.gz files
            if not Item[-7:] == ".txt.gz":
                continue
            Path = os.path.join(self.DirectoryOfPSMs, Item)
            FileStub = Item.replace(".txt.gz", "") #removes the extension
            #have to do some more cleanup to get back a good file name
            #print(FileStub)
            FileStub = FileStub.replace("_msgfdb_fht", "")
            FileStub = FileStub.replace("_msgfplus_fht", "")
            if not FileStub in self.FileToOrganismDictionary:
                print ("Can't associate this file %s with an organism."%Item)
                print ("In directory %s, but not listed in association file %s"%(self.DirectoryOfPSMs, self.FileToOrganismsPath))
                continue
            OrganismName = self.FileToOrganismDictionary[FileStub]
            if not os.path.isfile(Path):
                continue
            self.ParsePSMFile(Path, OrganismName)
        


    def ParseFileToOrganism(self, Path):
        #simple populating the FileToOrganismDictionary variable so that we can use this knowledge later on
        Handle = open(Path, 'r')
        Header = Handle.readline() #pop it off because I don't need it
        FileCounter = 0
        OrganismCounter = 0
        for Line in Handle:
            #[0] = organism name
            #[1] = comma separated list of file stubs
            Bits = Line.strip().split("\t")
            #make the organism object
            OrganismName = Bits[0]
            if not OrganismName in self.OrganismObjectsDictionary:
                self.OrganismObjectsDictionary[OrganismName] = Organism(OrganismName)
                OrganismCounter += 1
            #parse out the files
            Files = Bits[1].split(",")
            for File in Files:
                #probably have to strip out some white space here
                File = File.strip()
                self.FileToOrganismDictionary[File] = OrganismName
                FileCounter += 1
                #print (File, OrganismName)
        print("The associations file listed %s psm results files from %s organisms"%(FileCounter, OrganismCounter))

    def ParsePSMFile(self, Path, OrganismName):
        #print ("parsing %s"%Path)

        ##important subtle fact. these files are gziped. so we use the gzip reader
        Org = self.OrganismObjectsDictionary[OrganismName]
        Handle = gzip.open(Path, 'r')
        Header = Handle.readline()
        TargetPSMCount = 0
        DecoyPSMCount = 0
        DecoyPeptides = {} #key = peptide, value = spectra
        DecoyProteins = {} #key = protein, value = peptide
        for Line in Handle:
            decoded_Line = Line.decode() # the gzip reader puts things into 'byte' objects and not 'str'
            Bits = decoded_Line.strip().split("\t")
            # [4] is charge
            #[9] is the peptide string
            #[10] is protein
            #[14] is SPecEvalue
            #[17] is q value
            peptide = Bits[9][2:-2] #this last bit is to strip out the prefix/suffix
            qvalue = float(Bits[17])
            protein = Bits[10]
            #first a hard filter on the qvalue. We simply don't bother to look at crap
            if qvalue > self.QvalueCutoff:
                continue
            #now we don't want to include all the contaminants 
            if protein[:3] in ["Con"]:
                continue
            if protein[:3] in ["XXX"]:
                DecoyPSMCount += 1
            TargetPSMCount += 1
            ### Currently adding the entire protein fasta accession set
            ## sp|P54547|G6PD_BACSU
            ## ref|YP_001233830.1
            ## tr|J0JN52|J0JN52_ALCFA
            #print ("adding %s to %s"%(peptide, protein))
            Org.AddPeptide(peptide, protein)
        #print ("Target PSMs: %s, Decoy PSMs: %s"%(TargetPSMCount, DecoyPSMCount))


       
    def AddFiles(self, PSM_dir, FileToOrganisms):
        #this is called by the iPython Notebook to mimic the ParseCommandLine() function
        #but from the notebook
        self.FileToOrganismsPath = FileToOrganisms
        self.DirectoryOfPSMs = PSM_dir


    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:o:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                self.DirectoryOfPSMs = Value
            if Option == "-o":
                self.FileToOrganismsPath = Value 
        
if __name__ == "__main__":
    Robot = ParserClass()
    Robot.ParseCommandLine(sys.argv[1:])
    Robot.Main()
