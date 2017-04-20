#! /usr/bin/env python


################################################################################
# Copyright (c) 2009 Christopher M. Fraser and L. Ridgway Scott                #
# All Rights Reserved.                                                         #
#                                                                              #
# This software is the confidential and proprietary information of Christopher #
# M. Fraser and L. Ridgway Scott ("Confidential Information"). You shall not   #
# disclose such Confidential Information and shall use it only in accordance   #
# with the terms of the license agreement you entered into with Christopher M. #
# Fraser and L. Ridgway Scott. The use, reproduction, distribution, or sale of #
# this software either in whole or in part without the written consent of      #
# Christopher M. Fraser and L. Ridgway Scott is strictly prohibited.           #
#                                                                              #
# CHRISTOPHER M. FRASER AND L. RIDGWAY SCOTT MAKE NO REPRESENTATIONS OR        #
# WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, #
# INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. NEITHER CHRISTOPHER   #
# M. FRASER NOR L. RIDGWAY SCOTT SHALL BE LIABLE FOR ANY DAMAGES SUFFERED BY   #
# LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE       #
# OR ITS DERIVATIVES.                                                          #
################################################################################


################################################################################
# DESCRIPTION:                                                                 #
# This program takes as input one PDB file and prints the minimum distance     #
# between neighboring CA's to look for possible cis configurations             #
#                  python cadistmin.py 1NGL                                    #
################################################################################


################################################################################
# IMPORTS #
###########

import string
import os
import sys
import cStringIO
import math
import copy

################################################################################
# METHODS #
###########

#BEGIN: extractPDBAtomRecord()
def extractPDBAtomRecord(pdbRecord):
    if len(pdbRecord) >= 8:
        if  (pdbRecord[0:3] == "TER") or (pdbRecord[0:4] == "ATOM") or (pdbRecord[0:6] == "HETATM"):
            return {"recordType": pdbRecord[0:6].strip(),   \
                    "atomSerial": pdbRecord[6:11].strip(),  \
                    "atomName":   pdbRecord[12:16].strip(), \
                    "altLoc":     pdbRecord[16:17].strip(),    \
                    "resName":    pdbRecord[17:20].strip(), \
                    "chainID":    pdbRecord[21:22].strip(),    \
                    "resSeq":     pdbRecord[22:26].strip(), \
                    "iCode":      pdbRecord[26:27].strip(),    \
                    "x":          pdbRecord[30:38].strip(), \
                    "y":          pdbRecord[38:46].strip(), \
                    "z":          pdbRecord[46:54].strip(), \
                    "occupancy":  pdbRecord[54:60].strip(), \
                    "tempFactor": pdbRecord[60:66].strip(), \
                    "element":    pdbRecord[76:78].strip(), \
                    "charge":     pdbRecord[78:80].strip()  \
                   }
#END:  extractPDBAtomRecord()

#BEGIN: retrievePDBRecords
def retrievePDBRecords(pdbFile):
    pdbAtomRecords = []
    model = "0"
    for pdbRecord in pdbFile:
        if len(pdbRecord) >= 14:
            if pdbRecord[0:5] == "MODEL":
                model = pdbRecord[10:14].strip()
            else:
                if (pdbRecord[0:4] == "ATOM") or (pdbRecord[0:6] == "HETATM"):
                    pdbAtomRecord = extractPDBAtomRecord(pdbRecord)
                    pdbAtomRecord["model"] = str(model)
                    pdbAtomRecords.append(pdbAtomRecord)
    return pdbAtomRecords

#BEGIN: extractPDBhelixRecord()
def extractPDBhelixRecord(pdbRecord):
    if len(pdbRecord) >= 8:
        if pdbRecord[0:5] == "HELIX":
            return {"recordType": pdbRecord[0:5].strip(),   \
                    "helixNumbr": pdbRecord[6:10].strip(),  \
                    "helixName":  pdbRecord[11:14].strip(), \
                    "altLocBeg":  pdbRecord[14].strip(),    \
                    "resNameBeg": pdbRecord[15:18].strip(), \
                    "chainIDBeg": pdbRecord[19].strip(),    \
                    "resSeqBeg":  pdbRecord[20:25].strip(), \
                    "iCodeBeg":   pdbRecord[25].strip(), \
                    "altLocEnd":  pdbRecord[26].strip(),    \
                    "resNameEnd": pdbRecord[27:30].strip(), \
                    "chainIDEnd": pdbRecord[31].strip(),    \
                    "resSeqEnd":  pdbRecord[32:37].strip(), \
                    "iCodeEnd":   pdbRecord[37].strip(), \
                    "SheetType":  pdbRecord[38:40].strip(),    \
                   }
#END:  extractPDBhelixRecord()

#BEGIN: extractPDBsheetRecord()
def extractPDBsheetRecord(pdbRecord):
    if len(pdbRecord) >= 8:
        if pdbRecord[0:5] == "SHEET":
            return {"recordType": pdbRecord[0:5].strip(),   \
                    "sheetNumbr": pdbRecord[6:10].strip(),  \
                    "sheetName":  pdbRecord[11:16].strip(), \
                    "altLocBeg":  pdbRecord[16].strip(), \
                    "resNameBeg": pdbRecord[17:20].strip(),\
                    "chainIDBeg": pdbRecord[21].strip(), \
                    "resSeqBeg":  pdbRecord[22:26].strip(),\
                    "iCodeBeg":   pdbRecord[26].strip(), \
                    "altLocEnd":  pdbRecord[27].strip(), \
                    "resNameEnd": pdbRecord[28:31].strip(),\
                    "chainIDEnd": pdbRecord[32].strip() ,\
                    "resSeqEnd":  pdbRecord[33:37].strip(),\
                    "iCodeEnd":   pdbRecord[37].strip(), \
                    "SheetType":  pdbRecord[38:40].strip(),    \
                   }
#END:  extractPDBsheetRecord()

#BEGIN: retrievePDBHelices
def retrievePDBHelices(pdbFile):
    pdbHelixRecords = []
    for pdbRecord in pdbFile:
        if len(pdbRecord) >= 81:
            if pdbRecord[0:5] == "HELIX":
                pdbHelixRecord = extractPDBhelixRecord(pdbRecord)
                pdbHelixRecords.append(pdbHelixRecord)
    return pdbHelixRecords
#END: retrievePDBHelices

#BEGIN: retrievePDBSheets
def retrievePDBSheets(pdbFile):
    pdbSheetRecords = []
    for pdbRecord in pdbFile:
        if len(pdbRecord) >= 8:
            if pdbRecord[0:5] == "SHEET":
                pdbSheetRecord = extractPDBsheetRecord(pdbRecord)
                pdbSheetRecords.append(pdbSheetRecord)
    return pdbSheetRecords
#END: retrievePDBSheets

#BEGIN: iCodeOrder
def iCodeOrder(icname):
    ico=12
    if (icname==""):  ico=0
    if (icname=="A"): ico=1
    if (icname=="B"): ico=2
    if (icname=="C"): ico=3
    if (icname=="D"): ico=4
    if (icname=="E"): ico=5
    if (icname=="F"): ico=6
    if (icname=="G"): ico=7
    if (icname=="H"): ico=8
    if (icname=="I"): ico=9
    if (icname=="J"): ico=10
    if (icname=="X"): ico=11
    return ico
#END: iCodeOrder

#BEGIN: HelixCall
def HelixCall(pdbrec,pdbSSdb):
    call="L"
    for i in range(0, len(pdbSSdb)):
        if (pdbrec["chainID"]==pdbSSdb[i]["chainIDBeg"]) and\
           (eval(pdbrec["resSeq"])+0.1*iCodeOrder(pdbrec["iCode"]) >= \
            eval(pdbSSdb[i]["resSeqBeg"])+0.1*iCodeOrder(pdbSSdb[i]["iCodeBeg"])) and\
           (eval(pdbrec["resSeq"])+0.1*iCodeOrder(pdbrec["iCode"]) <= \
            eval(pdbSSdb[i]["resSeqEnd"])+0.1*iCodeOrder(pdbSSdb[i]["iCodeEnd"])):
            call="H"
#           if iCodeOrder(pdbSSdb[i]["iCodeBeg"])>0:
#               print "iCode",pdbrec["resSeq"],pdbrec["iCode"],call,iCodeOrder(pdbrec["iCode"])
    return call
#END: HelixCall

#BEGIN: SheetCall
def SheetCall(pdbrec,pdbSSdb):
    call="L"
    for i in range(0, len(pdbSSdb)):
        if (pdbrec["chainID"]==pdbSSdb[i]["chainIDBeg"]) and\
           (eval(pdbrec["resSeq"])+0.1*iCodeOrder(pdbrec["iCode"]) >= \
            eval(pdbSSdb[i]["resSeqBeg"])+0.1*iCodeOrder(pdbSSdb[i]["iCodeBeg"])) and\
           (eval(pdbrec["resSeq"])+0.1*iCodeOrder(pdbrec["iCode"]) <= \
            eval(pdbSSdb[i]["resSeqEnd"])+0.1*iCodeOrder(pdbSSdb[i]["iCodeEnd"])):
#print pdbrec["model"],pdbrec["chainID"],pdbrec["resSeq"],pdbSSdb[i]["resSeqBeg"],pdbSSdb[i]["resSeqEnd"]
            call="S"
#           if iCodeOrder(pdbSSdb[i]["iCodeBeg"])>0:
#           print "iCode",pdbrec["resSeq"],pdbrec["iCode"],call,iCodeOrder(pdbrec["iCode"]),\
#               iCodeOrder(pdbSSdb[i]["iCodeBeg"]), iCodeOrder(pdbSSdb[i]["iCodeEnd"]),\
#               "resSeqBeg",pdbSSdb[i]["resSeqBeg"], "iCodeBeg",iCodeOrder(pdbSSdb[i]["iCodeBeg"]),\
#               "resSeqEnd",pdbSSdb[i]["resSeqEnd"], "iCodeEnd",iCodeOrder(pdbSSdb[i]["iCodeEnd"])
    return call
#END: SheetCall

#BEGIN: SecStrucCall
def SecStrucCall(pdbrec,pdbHelices,pdbSheets):
    call="L"
    ch=HelixCall(pdbrec,pdbHelices)
    if (ch=="H"): call="H"
    cs=SheetCall(pdbrec,pdbSheets)
    if (cs=="S"): call="S"
    if (ch=="H") and (cs=="S"): 
        call="F"
#       print "iCodecall",pdbrec["resName"],pdbrec["chainID"], pdbrec["resSeq"],pdbrec["iCode"]
    return call
#END: SecStrucCall



########################################################################################################################################

########################################################################################################################################

