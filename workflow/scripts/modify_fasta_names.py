#!/usr/bin/python
import sys
import re
from fasta import readfasta

fasta = open(sys.argv[1],"r")

################################################################################
# Read the fasta file into a dictionary
################################################################################
fastaDict1 = readfasta(fasta)

# sp_remove = ["P26503_106", "P26503_118", "P28566_101", "P28566_107", "P28566_108", "P28566_114"]

# names = {"P26503_101":"VOL_Bys06",
# "P26503_102":"VOLAUR_K1013",
# "P26503_103":"PANMOR_K1379",
# "P26503_104":"PAN_Kgh192",
# "P26503_105":"EUDELE_NIVA14",
# "P26503_107":"CHLSCH_Pil01",
# "P26503_108":"CHLPET_SAG10_73",
# "P26503_109":"CHLZEB_SAG10_83",
# "P26503_110":"CHLASY_SAG11-41",
# "P26503_111":"VIT_SAG12_93",
# "P26503_112":"CHLCRI_SAG13_72",
# "P26503_113":"CHLPAR_SAG14_88",
# "P26503_114":"PAUPSE_SAG17_94",
# "P26503_117":"CHLFOT_SAG21_83",
# "P26503_119":"CHLPIL_SAG39_72",
# "P26503_120":"LOBROS_SAG45-1",
# "P26503_122":"TETCYL_SAG63_94",
# "P26503_123":"VITAUL_SAG69_72",
# "P26503_124":"CHLGLO_SAG7_73",
# "P26503_125":"CHLLAT_SAG70_81",
# "P26503_126":"VOLTER_SAG88-3",
# "P26503_127":"CHLCAL_SAG9_72",
# "P28566_102":"TETSOC_Nbe02",
# "P28566_103":"CHLPAR_Nbe06",
# "P28566_104":"MICGLO_Lie01",
# "P28566_105":"CHLTYP_NIVA21",
# "P28566_106":"CHLPAR_Osg09",
# "P28566_110":"VITORD_Kgh02",
# "P28566_112":"CHLPAR_Osg16",
# "P28566_113":"HAELAC_NIV9",
# "P28566_115":"CHLISA_Bae13", # This one is merged with P29912_115
# "P28566_116":"CHLsp_Kgh28",
# "P28566_117":"VOLGLO_SAG199_80",
# "P28566_119":"CHLTYP_SAG61_72",
# "P28566_120":"VOLVBAR_NIES730",
# "P28566_121":"VOLGIG_NIES867",
# "P29912_105":"CHLTYP_NIVA21",
# "S104_merge":"MICGLO_Lie01",
# "S115_merge":"CHLISA_Bae13"}

for key in fastaDict1.keys():
    if "--" in key:
        seq_name = key.split("--")[0]#.replace(">", "")
        print(">"+seq_name)
        print(fastaDict1[key])
    else:
        seq_name = re.sub(r'^[^_]+_', '>', key)
        print(seq_name)
        print(fastaDict1[key])

