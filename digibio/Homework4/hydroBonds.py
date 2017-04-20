#! /usr/bin/env python

import hsdistV2 as HSD
import numpy as np
import string
import os
import sys
import cStringIO
import math
import copy


class Bond:
    def __init__(self):
        self.one = False
        self.two = False
        self.three = False
        self.alpha = False
        self.pi = False
        self.six = False
    def addOne(self):
        self.one = True
    def addTwo(self):
        self.two = True
    def addThree(self):
        self.three = True
    def addAlpha(self):
        self.alpha = True
    def addPi(self):
        self.pi = True
    def addSix(self):
        self.six = True
    def addBond(self, dist):
        if dist == 1:
            self.addOne()
        elif dist == 2:
            self.addTwo()
        elif dist == 3:
            self.addThree()
        elif dist == 4:
            self.addAlpha()
        elif dist == 5:
            self.addPi()
        elif dist == 6:
            self.addSix()
        else:
            print("Bond out of range")
    def toString(self):
        return "+1: "+str(self.one)+" +2: "+str(self.two)+" +3: "+str(self.three)+" +4: "+str(self.alpha)+" +5: "+str(self.pi)+" +6: "+str(self.six)  
    def toTable(self):
        return "\t+1: "+boolToSymbol(self.one)+"\t+2: "+boolToSymbol(self.two)+"\t+3: "+boolToSymbol(self.three)+"\t+4: "+boolToSymbol(self.alpha)+"\t+5: "+boolToSymbol(self.pi)+"\t+6: "+boolToSymbol(self.six)

class SheetBond:
    def __init__(self,ori_):
        self.distance = list()
        self.ori = ori_
    def addDist(self, distance):
        self.distance.insert(0, distance)
    def toString(self):
        return "Orientation: "+str(self.ori)+" tDistances: "+str(self.distance)  
    def toTable(self):
        return "\tOrientation: "+str(self.ori)+"\tDistances: +"+str(self.distance)

# given a boolean,
# returns an X if True,
# or an O if False
def boolToSymbol(b):
    if b:
        return "X"
    else:
        return "O"

# given two 3D points A and B,
# returns the distance between the two points
def dist(a,b):
    distance = math.sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2])) 
    return distance

# given three 3D points A, B, and C,
# returns the angle between them in degrees
def angle(a, b, c): #find angle ABC
    a = np.array([a[0], a[1], a[2]])
    b = np.array([b[0], b[1], b[2]])
    c = np.array([c[0], c[1], c[2]])
    b_a = a - b
    b_c = c - b

    angle = np.arccos(np.dot(b_a,b_c) / (np.linalg.norm(b_a) * np.linalg.norm(b_c)))
    return np.degrees(angle)

# given 4 vectors, 
# returns a boolean describing if they pass the 4 criteria for hydrogen bonding
def isBond(d,h,a,b):
    test_one = (dist(d,a) < 3.5)
    test_two = (dist(h,a) < 2.5)
    test_three = (angle(d, h, a) > 90)
    test_four = (angle(d, a, b) > 90)
    test_five = (angle(h, a, b) > 90)
    return test_one & test_two & test_three & test_four & test_five


'''
# given pdb records and helices, 
# prints each hydrogen bond found and adds them to a dictionary
def countHelices(pdbRecs, pdbHelices):
    count = 0
    acceptorCount = 0
    bonds_formed = dict()

    num_res = int(pdbRecs[len(pdbRecs) -1]["resSeq"])
    for i in range(0, len(pdbHelices)):

        startLoc = int(pdbHelices[i]["resSeqBeg"]) 
        endLoc = int(pdbHelices[i]["resSeqEnd"]) 
        print ("START: %d        END: %d") % (startLoc, endLoc)
        
        for j in range (startLoc, endLoc + 1): #loop through the residues of the helix
            print("    TESTING residue %d for donor") % j
            donorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == j)] #get all of the atoms in the current residue j
            donorAtom = [rec for rec in donorRes if rec["atomName"] == "O"] #pull out O (donor)
            donor = [float(donorAtom[0]["x"]), float(donorAtom[0]["y"]), float(donorAtom[0]["z"])]
            antecedentAtom = [rec for rec in donorRes if rec["atomName"] == "C"] #pull out C (antecedent)
            antecedent = [float(antecedentAtom[0]["x"]), float(antecedentAtom[0]["y"]), float(antecedentAtom[0]["z"])] 
            
            
            for k in range (0, num_res):#(startLoc, endLoc + 1):#loop through the rest of the residues to see if their acceptor/hydrogen form hydrogen bonds
                acceptorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == k)]
                acceptorAtom = [rec for rec in acceptorRes if rec["atomName"] == "N"] #pull out N (acceptor)
                
                if((len(acceptorAtom) > 0) and (int(donorAtom[0]["resSeq"]) != (int(acceptorAtom[0]["resSeq"])))):
                    acceptor = [float(acceptorAtom[0]["x"]), float(acceptorAtom[0]["y"]), float(acceptorAtom[0]["z"])] 
                    hydrogenAtom = [rec for rec in acceptorRes if rec["atomName"] == "H"] #pull out H (hydrogen)
                    
                    if (len(hydrogenAtom) > 0):
                        hydrogen = [float(hydrogenAtom[0]["x"]), float(hydrogenAtom[0]["y"]), float(hydrogenAtom[0]["z"])] 
                        
                        if(isBond(donor, hydrogen, acceptor, antecedent)):
                            dist = int(acceptorAtom[0]["resSeq"]) - int(donorAtom[0]["resSeq"])
                            key = int(donorAtom[0]["resSeq"]) #edited
                            if not bonds_formed.get(key):
                                acceptorCount += 1
                                new_bond = Bond()
                                new_bond.addBond(dist)
                                bonds_formed.update({key:copy.deepcopy(new_bond)})
                            else:
                                edited_bond = bonds_formed[key]
                                edited_bond.addBond(dist)
                                bonds_formed.update({key:copy.deepcopy(edited_bond)})
                            print("   *** Hydrogen bond formed with +%d") % dist
                            print("    *Donor res num: %d") % (int(donorAtom[0]["resSeq"]))
                            print("    *Acceptor res num: %d") % key
                            count += 1
    print("TOTAL NUMBER OF HELICES: %d") % len(pdbHelices)
    print("TOTAL NUMBER OF ACCEPTORS: %d") % acceptorCount
    print("TOTAL NUMBER OF HYDROGEN BONDS: %d") % count
    
    for key in bonds_formed.keys():
        if key:
            print(str(key)+bonds_formed[key].toTable())
        else:
            print("Error happened")
'''

# given pdb records and helices, 
# prints each hydrogen bond found and adds them to a dictionary
def countHelices(filename, pdbRecs, pdbHelices):
    with open(filename+'.txt', 'w+') as f:
    

        count = 0
        acceptorCount = 0
        bonds_formed = dict()
        acceptorThreeCount = 0
        acceptorFourCount = 0
        acceptorFiveCount = 0
        acceptorBothCount = 0

        num_res = int(pdbRecs[len(pdbRecs) -1]["resSeq"])
        for i in range(0, len(pdbHelices)):

            startLoc = int(pdbHelices[i]["resSeqBeg"]) 
            endLoc = int(pdbHelices[i]["resSeqEnd"]) 
            print ("START: %d        END: %d") % (startLoc, endLoc)
        
            for j in range (startLoc, endLoc + 1): #loop through the residues of the helix
                print("    TESTING residue %d for donor") % j
                acceptorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == j)] #get all of the atoms in the current residue j
                acceptorAtom = [rec for rec in acceptorRes if rec["atomName"] == "O"] #pull out O (donor)
                if acceptorAtom:
                    acceptor = [float(acceptorAtom[0]["x"]), float(acceptorAtom[0]["y"]), float(acceptorAtom[0]["z"])]
                else:
                    continue
                antecedentAtom = [rec for rec in acceptorRes if rec["atomName"] == "C"] #pull out C (antecedent)
                antecedent = [float(antecedentAtom[0]["x"]), float(antecedentAtom[0]["y"]), float(antecedentAtom[0]["z"])] 
            
            
                for k in range (num_res):#(startLoc, endLoc + 1):#loop through the rest of the residues to see if their acceptor/hydrogen form hydrogen bonds
                    donorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == k)]
                    donorAtom = [rec for rec in donorRes if rec["atomName"] == "N"] #pull out N (acceptor)
                
                    if((len(donorAtom) > 0) and (int(donorAtom[0]["resSeq"]) != (int(acceptorAtom[0]["resSeq"])))):
                        donor = [float(donorAtom[0]["x"]), float(donorAtom[0]["y"]), float(donorAtom[0]["z"])] 

                        hydrogenAtom = [rec for rec in donorRes if rec["atomName"] == "H"] #pull out H (hydrogen)
                    
                        if (len(hydrogenAtom) > 0):
                            hydrogen = [float(hydrogenAtom[0]["x"]), float(hydrogenAtom[0]["y"]), float(hydrogenAtom[0]["z"])] 
                            

                            if(isBond(donor, hydrogen, acceptor, antecedent)):
                                dist = math.fabs(int(donorAtom[0]["resSeq"]) - int(acceptorAtom[0]["resSeq"]))
                                key = int(acceptorAtom[0]["resSeq"])
                                if not bonds_formed.get(key):
                                    acceptorCount += 1
                                    new_bond = Bond()
                                    new_bond.addBond(dist)
                                    if dist == 3:
                                        acceptorThreeCount += 1
                                    elif dist == 4:
                                        acceptorFourCount += 1
                                    elif dist == 5:
                                        acceptorFiveCount += 1
                                    bonds_formed.update({key:copy.deepcopy(new_bond)})
                                else:
                                    edited_bond = bonds_formed[key]
                                    edited_bond.addBond(dist)
                                    if dist == 3:
                                        acceptorThreeCount += 1
                                    elif dist == 4:
                                        acceptorFourCount += 1
                                    elif dist == 5:
                                        acceptorFiveCount += 1
                                    bonds_formed.update({key:copy.deepcopy(edited_bond)})
                                print("   *** Hydrogen bond formed with +%d") % dist
                                print("    *Donor res num: %d") % (int(donorAtom[0]["resSeq"]))
                                print("    *Acceptor res num: %d") % key
                                count += 1
        f.write("Number of helices searched through: "+str(len(pdbHelices))+"\n")
        f.write("Number of distinct acceptors: "+str(acceptorCount)+"\n\n")
        f.write("Total number of hydrogen bonds found: "+str(count)+"\n")
        f.write("Number of 3_10 bonds: "+str(acceptorThreeCount)+"\n")
        f.write("Number of alpha bonds: "+str(acceptorFourCount)+"\n")
        f.write("Number of pi bonds: "+str(acceptorFiveCount)+"\n")
        
        for key in bonds_formed.keys():
            if (bonds_formed[key].three and bonds_formed[key].alpha) or (bonds_formed[key].alpha and bonds_formed[key].pi) or (bonds_formed[key].three and bonds_formed[key].pi):
                acceptorBothCount += 1
        f.write("Number of acceptors with double bonds formed: "+str(acceptorBothCount)+"\n\n")

        for key in bonds_formed.keys():
            if key:
                f.write(str(key)+bonds_formed[key].toTable()+"\n")
            else:
                print("Error happened")



# given pdb records and sheets, 
# prints each hydrogen bond found and adds them to a dictionary
def countSheets(filename, pdbRecs, pdbSheets):
    with open(filename+'B.txt', 'w+') as f:

        count = 0
        acceptorCount = 0
        bonds_formed = dict()
        acceptorBothCountPara = 0
        acceptorBothCountAnti = 0

        num_res = int(pdbRecs[len(pdbRecs) -1]["resSeq"])
        for i in range(0, len(pdbSheets)):

            startLoc = int(pdbSheets[i]["resSeqBeg"]) 
            endLoc = int(pdbSheets[i]["resSeqEnd"]) 
            orientation = int(pdbSheets[i]["SheetType"])
            print ("START: %d        END: %d") % (startLoc, endLoc)
        
            for j in range (startLoc, endLoc + 1): #loop through the residues of the helix
                print("    TESTING residue %d for donor") % j
                acceptorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == j)] #get all of the atoms in the current residue j
                acceptorAtom = [rec for rec in acceptorRes if rec["atomName"] == "O"] #pull out O (donor)
                if acceptorAtom:
                    acceptor = [float(acceptorAtom[0]["x"]), float(acceptorAtom[0]["y"]), float(acceptorAtom[0]["z"])]
                else:
                    continue
                antecedentAtom = [rec for rec in acceptorRes if rec["atomName"] == "C"] #pull out C (antecedent)
                antecedent = [float(antecedentAtom[0]["x"]), float(antecedentAtom[0]["y"]), float(antecedentAtom[0]["z"])] 
            
            
                for k in range (0, num_res):#(startLoc, endLoc + 1):#loop through the rest of the residues to see if their acceptor/hydrogen form hydrogen bonds
                    donorRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == k)]
                    donorAtom = [rec for rec in donorRes if rec["atomName"] == "N"] #pull out N (acceptor)
                
                    if((len(donorAtom) > 0) and (int(donorAtom[0]["resSeq"]) != (int(acceptorAtom[0]["resSeq"])))):
                        donor = [float(donorAtom[0]["x"]), float(donorAtom[0]["y"]), float(donorAtom[0]["z"])] 

                        hydrogenAtom = [rec for rec in donorRes if rec["atomName"] == "H"] #pull out H (hydrogen)
                    
                        if (len(hydrogenAtom) > 0):
                            hydrogen = [float(hydrogenAtom[0]["x"]), float(hydrogenAtom[0]["y"]), float(hydrogenAtom[0]["z"])] 
                            
                            if(isBond(donor, hydrogen, acceptor, antecedent)):
                                dist = math.fabs(int(donorAtom[0]["resSeq"]) - int(acceptorAtom[0]["resSeq"]))
                                key = int(acceptorAtom[0]["resSeq"]) #edited
                                if not bonds_formed.get(key):
                                    acceptorCount += 1
                                    new_bond = SheetBond(orientation)
                                    new_bond.addDist(dist)
                                    bonds_formed[key] = new_bond
                                else:
                                    edited_bond = bonds_formed[key]
                                    edited_bond.addDist(dist)
                                    bonds_formed[key] = edited_bond                             
                                print("   *** Hydrogen bond formed with +%d") % dist
                                print("    *Donor res num: %d") % (int(donorAtom[0]["resSeq"]))
                                print("    *Acceptor res num: %d") % key
                                count += 1
        f.write("Number of sheets searched through: "+str(len(pdbSheets))+"\n")
        f.write("Number of distinct acceptors: "+str(acceptorCount)+"\n\n")
        f.write("Total number of hydrogen bonds found: "+str(count)+"\n")

        for key in bonds_formed.keys():
            if len(bonds_formed[key].distance) > 1:
                if bonds_formed[key].ori == 1:
                    acceptorBothCountPara += 1
                elif bonds_formed[key].ori == -1:
                    acceptorBothCountAnti += 1
        f.write("Number of acceptors with double bonds formed over parallel sheets: "+str(acceptorBothCountPara)+"\n\n")
        f.write("Number of acceptors with double bonds formed over antiparallel sheets: "+str(acceptorBothCountAnti)+"\n\n")
    
        for key in bonds_formed.keys():
            if key:
                f.write(str(key)+bonds_formed[key].toTable()+"\n")
            else:
                print("Error happened")


def countAmines(filename,pdbRecs,pdbHelices):
    
    f = open(filename+"A.txt", 'w+')
    countsThree = dict()
    countsAlpha = dict()
    countsPi = dict()

    countThree = 0
    countAlpha = 0
    for i in range(len(pdbHelices)):
        helixType = int(pdbHelices[i]["SheetType"])
        startLocation = int(pdbHelices[i]["resSeqBeg"])
        endLocation = int(pdbHelices[i]["resSeqEnd"])
   
        print(str(helixType))
 
        for j in range(startLocation, endLocation + 1):
             currRes = [rec for rec in pdbRecs if (int(rec["resSeq"]) == j)]
             currResName = None
             if len(currRes) > 0:
                 currResName = currRes[0]["resName"]
             else:
                 print("no residue here")
     
             if helixType == 1:
                 if not countsAlpha.get(currResName):
                     countsAlpha.update({currResName:1})
                 else:
                     countsAlpha[currResName] = countsAlpha[currResName] + 1
                 countAlpha += 1
             elif helixType == 5:
                 if not countsThree.get(currResName):
                     countsThree.update({currResName:1})               
                 else:
                     countsThree[currResName] = countsThree[currResName] + 1
                 countThree += 1
    f.write("Alpha helices counts: "+str(countAlpha)+"\n")
    for key in countsAlpha.keys():
        f.write(key+" "+str(countsAlpha[key])+" "+str(float(countsAlpha[key])/float(countAlpha))+"\n")

    f.write("3_10 helices counts: "+str(countThree)+"\n")
    for key in countsThree.keys():
        f.write(key+" "+str(countsThree[key])+" "+str(float(countsThree[key])/float(countThree))+"\n")
    f.close()
    

files = ("1BGK_A","1BHU_A","2A2B_A","2A3D_A","2AB3_A","2AJJ_A","2AJW_A","2AK0_A","2AMN_A", "2B7E_A","2BBL_A","2BF9_A","2BIC_A","2BL6_A", "2CBP_A", "2CP8_A", "2CP9_A","2DKZ_A", "2JM2_A", "3NUL_A", "8TFV_A", "7AHL_B")

def main():
    # saves the pdb file argument
    #for filename in files:


    filename = sys.argv[1]
    dotpdb=".pdb"
    full_filename=filename+".pdb"

    	# opens the pdb file and scans for records, helices, and sheets
    pdbFile = open(full_filename, "r")
    pdbRecs = HSD.retrievePDBRecords(pdbFile)
    pdbFile.seek(0)
    pdbHelices = HSD.retrievePDBHelices(pdbFile)
    pdbFile.seek(0)
    pdbSheets = HSD.retrievePDBSheets(pdbFile)
    pdbFile.close()

    print("***********HELICES*****************")
    #countHelices(filename, pdbRecs, pdbHelices)
    print("************SHEETS*****************")
    countSheets(filename, pdbRecs, pdbSheets)
    #countAmines(filename, pdbRecs, pdbHelices)\

    
if __name__ == '__main__':
    main()

########################################################################################################################################

########################################################################################################################################

