
from Bio import SeqIO

VIC_HA_OR = ""
for record in SeqIO.parse("../data/reference/VIC_PlasmidControl_OR.fa", "fasta"):
    if str(record.id) == "HA":
        VIC_HA_OR = record.seq

YAM_HA_OR = ""
for record in SeqIO.parse("../data/reference/YAM_PlasmidControl_OR.fa", "fasta"):
    if str(record.id) == "HA":
        YAM_HA_OR = record.seq

VIC_HA = ""
for record in SeqIO.parse("../data/reference/VIC_PlasmidControl.fa", "fasta"):
    if str(record.id) == "HA":
        VIC_HA = record.seq

YAM_HA = ""
for record in SeqIO.parse("../data/reference/YAM_PlasmidControl.fa", "fasta"):
    if str(record.id) == "HA":
        YAM_HA = record.seq

VIC_HA_OR_aa = VIC_HA_OR.translate()
YAM_HA_OR_aa = YAM_HA_OR.translate()


# These are the antigenic sites based on alignment
VIC_120 = "HIRLSTHNVINAENAPGGPYKI"
VIC_150 = "GSCPNVTNGN"
VIC_160 = "KNDKNKTAT"
VIC_190 = "NETQMAKLY"

YAM_120 = "NIRLSTQNVIDAEKAPGGPYRL"
YAM_150 = "GSCPNATSKS"
YAM_160 = "KDNNKNAT"
YAM_190 = "NKTQMKNLY"

# First VIC
VIC_120_start = VIC_HA_OR_aa.find(VIC_120)
VIC_150_start = VIC_HA_OR_aa.find(VIC_150)
VIC_160_start = VIC_HA_OR_aa.find(VIC_160)
VIC_190_start = VIC_HA_OR_aa.find(VIC_190)

VIC_120_end = VIC_120_start + len(VIC_120)
VIC_150_end = VIC_150_start + len(VIC_150)
VIC_160_end = VIC_160_start + len(VIC_160)
VIC_190_end = VIC_190_start + len(VIC_190)

VIC_120_nt = VIC_HA_OR[VIC_120_start*3:VIC_120_end*3]
VIC_150_nt = VIC_HA_OR[VIC_150_start*3:VIC_150_end*3]
VIC_160_nt = VIC_HA_OR[VIC_160_start*3:VIC_160_end*3]
VIC_190_nt = VIC_HA_OR[VIC_190_start*3:VIC_190_end*3]

VIC_120_start_abs = VIC_HA.find(VIC_120_nt)
VIC_150_start_abs = VIC_HA.find(VIC_150_nt)
VIC_160_start_abs = VIC_HA.find(VIC_160_nt)
VIC_190_start_abs = VIC_HA.find(VIC_190_nt)

VIC_120_end_abs = VIC_120_start_abs + len(VIC_120_nt)
VIC_150_end_abs = VIC_150_start_abs + len(VIC_150_nt)
VIC_160_end_abs = VIC_160_start_abs + len(VIC_160_nt)
VIC_190_end_abs = VIC_190_start_abs + len(VIC_190_nt)

# Now YAM
YAM_120_start = YAM_HA_OR_aa.find(YAM_120)
YAM_150_start = YAM_HA_OR_aa.find(YAM_150)
YAM_160_start = YAM_HA_OR_aa.find(YAM_160)
YAM_190_start = YAM_HA_OR_aa.find(YAM_190)

YAM_120_end = YAM_120_start + len(YAM_120)
YAM_150_end = YAM_150_start + len(YAM_150)
YAM_160_end = YAM_160_start + len(YAM_160)
YAM_190_end = YAM_190_start + len(YAM_190)

YAM_120_nt = YAM_HA_OR[YAM_120_start*3:YAM_120_end*3]
YAM_150_nt = YAM_HA_OR[YAM_150_start*3:YAM_150_end*3]
YAM_160_nt = YAM_HA_OR[YAM_160_start*3:YAM_160_end*3]
YAM_190_nt = YAM_HA_OR[YAM_190_start*3:YAM_190_end*3]

YAM_120_start_abs = YAM_HA.find(YAM_120_nt)
YAM_150_start_abs = YAM_HA.find(YAM_150_nt)
YAM_160_start_abs = YAM_HA.find(YAM_160_nt)
YAM_190_start_abs = YAM_HA.find(YAM_190_nt)

YAM_120_end_abs = YAM_120_start_abs + len(YAM_120_nt)
YAM_150_end_abs = YAM_150_start_abs + len(YAM_150_nt)
YAM_160_end_abs = YAM_160_start_abs + len(YAM_160_nt)
YAM_190_end_abs = YAM_190_start_abs + len(YAM_190_nt)

header = "Lineage,Region,Start,End\n"
line1 = "B/VIC,120_loop," + str(VIC_120_start_abs) + "," + str(VIC_120_end_abs) + "\n"
line2 = "B/VIC,150_loop," + str(VIC_150_start_abs) + "," + str(VIC_150_end_abs) + "\n"
line3 = "B/VIC,160_loop," + str(VIC_160_start_abs) + "," + str(VIC_160_end_abs) + "\n"
line4 = "B/VIC,190_helix," + str(VIC_190_start_abs) + "," + str(VIC_190_end_abs) + "\n"
line5 = "B/YAM,120_loop," + str(YAM_120_start_abs) + "," + str(YAM_120_end_abs) + "\n"
line6 = "B/YAM,150_loop," + str(YAM_150_start_abs) + "," + str(YAM_150_end_abs) + "\n"
line7 = "B/YAM,160_loop," + str(YAM_160_start_abs) + "," + str(YAM_160_end_abs) + "\n"
line8 = "B/YAM,190_helix," + str(YAM_190_start_abs) + "," + str(YAM_190_end_abs) + "\n"

with open("../data/processed/antigenic_positions.csv", 'a') as pos_file:
    pos_file.write(header)
    pos_file.write(line1)
    pos_file.write(line2)
    pos_file.write(line3)
    pos_file.write(line4)
    pos_file.write(line5)
    pos_file.write(line6)
    pos_file.write(line7)
    pos_file.write(line8)
