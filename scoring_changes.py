#Score Base Changes
#If Transitions: T --> C  or G --> A

# If Transversions:
# T --> A or G
# C --> A or G
# G --> T or C
# A --> T or C

# High score means tolerable

#Checks for transition. If not transition, means it is transversion
def check_transition (previous_base, new_base):
  if previous_base == "T" or previous_base == "C":
    if new_base == "T" or new_base == "C":
      return True
  else:
    if new_base == "G" or new_base == "A":
      return True
  return False

#Checks the tolerability of aminoacid change according to the chemical properties 

def chemical (previousaa, newaa):
  chemical = {"nonpolar": ["G","A","V","C","P","L","I","M","W","F"], "polar" : ["S","T","Y","N","Q"], "neg" : ["D","E"], "pos" : ["K","R","H"]}
  if previousaa in chemical["nonpolar"]:
    if newaa in chemical["nonpolar"]:
      return True
  elif previousaa in chemical["polar"]:
    if newaa in chemical["polar"]:
      return True
  elif previousaa in chemical["neg"]:
    if newaa in chemical["neg"]:
      return True
  elif previousaa in chemical["pos"]:
    if newaa in chemical["pos"]:
      return True
  return False

#Contains the properties of VUS

VUS = {"Variant #0000325290" : {"previousaa" : "V", "newaa" : "A" , "previousbase" : "T" , "newbase" : "C"},
       "Variant #0000255886" : {"previousaa" : "D", "newaa" : "G" , "previousbase" : "A" , "newbase" : "G"},
       "Variant #0000262212" : {"previousaa" : "D", "newaa" : "E" , "previousbase" : "T" , "newbase" : "G"},
       "Variant #0000560653" : {"previousaa" : "A", "newaa" : "P" , "previousbase" : "G" , "newbase" : "C"},
       "Variant #0000560656" : {"previousaa" : "E", "newaa" : "G" , "previousbase" : "A" , "newbase" : "G"},
       "Variant #0000807853" : {"previousaa" : "I", "newaa" : "V" , "previousbase" : "A" , "newbase" : "G"},
       "Variant #0000256069" : {"previousaa" : "I", "newaa" : "V" , "previousbase" : "A" , "newbase" : "G"},
       "Variant #0000325293" : {"previousaa" : "R", "newaa" : "C" , "previousbase" : "C" , "newbase" : "T"},
       "Variant #0000560657" : {"previousaa" : "S", "newaa" : "L" , "previousbase" : "C" , "newbase" : "T"}}

#Contains the scores of VUS
VUS_points = {"Variant #0000325290" : 0,
              "Variant #0000255886" : 0,
              "Variant #0000262212" : 0,
              "Variant #0000560653" : 0,
              "Variant #0000560656" : 0,
              "Variant #0000807853" : 0,
              "Variant #0000256069" : 0,
              "Variant #0000325293" : 0,
              "Variant #0000560657" : 0}

#Iterates over each VUS to assign a score

for variant,properties in VUS.items():
  if chemical(properties['previousaa'],properties["newaa"]):
    VUS_points[variant] += 1
    if check_transition(properties['previousbase'],properties["newbase"]):
      VUS_points[variant] += 1
    else:
      VUS_points[variant] -= 1
  else:
    VUS_points[variant] -= 1
    if check_transition(properties['previousbase'],properties["newbase"]):
      VUS_points[variant] += 1
    else:
      VUS_points[variant] -= 1


for variant, score in VUS_points.items():
    print(variant, " : ", score)