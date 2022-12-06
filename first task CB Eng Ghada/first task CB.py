from pyopenms import *
edb = ElementDB()
import math
import pyopenms
oxygen = edb.getElement("O")
#oxxygen infos
print(oxygen.getName())
print(oxygen.getSymbol())
print(oxygen.getMonoWeight())
print(oxygen.getAverageWeight())
print("One mole of oxygen weight", 2*oxygen.getAverageWeight(), "grams")

#sulfur
sulfur=edb.getElement("S")
print(sulfur.getName())
print(sulfur.getSymbol())
print(sulfur.getMonoWeight())
print(sulfur.getAverageWeight())
isotopes = sulfur.getIsotopeDistribution()
print("One mole of 16O2 weight", 2*sulfur.getAverageWeight(), "grams")



#oxygen isotubes and its masses
edb=ElementDB()
oxygen=edb.getElement("O")
isotopes=oxygen.getIsotopeDistribution()
oxygen_isoDist = {"mass": [], "abundance": []}

for iso in isotopes.getContainer():
    print("Oxygen isotope", iso.getMZ(),
          "has abundance", iso.getIntensity()*100, "%")
    oxygen_isoDist["mass"].append(iso.getMZ())
    oxygen_isoDist["abundance"].append((iso.getIntensity() * 100))

#sulufr istopes and its masses


Sulfur_isoDist = {"mass": [], "abundance": []}
Sulfur = edb.getElement("S")
isotopes = Sulfur.getIsotopeDistribution()
for iso in isotopes.getContainer():
    print("Sulfur isotope", iso.getMZ(),
          "has abundance", iso.getIntensity()*100, "%")
    Sulfur_isoDist["mass"].append(iso.getMZ())
    Sulfur_isoDist["abundance"].append((iso.getIntensity() * 100))


# graphs

import matplotlib.pyplot as plt
import math

def adjustText(x1, y1, x2, y2):
    if(y1 > y2):
        plt.annotate("%0.3f" % (y1), xy=(x1, y2), xytext=(x1+0.5, y1+9),
                     textcoords="data",
                     arrowprops=dict(arrowstyle="->", color="r", lw=0.5),
                     horizontalalignment="right", verticalalignment="top")

    else:
        plt.annotate("%0.3f" % (y1), xy=(x1, y1), xytext=(x1+0.5, y1+9),
                     textcoords="data",
                     arrowprops=dict(arrowstyle="->", color="r", lw=0.5),
                     horizontalalignment="right", verticalalignment="top")

def plotDistribution(distribution):
        n = len(distribution["mass"])
        for i in range(0, n):
            plt.vlines(x=distribution["mass"][i], ymin=0,
                       ymax=distribution["abundance"][i])
            if int(distribution["mass"][i - 1]) == int(distribution["mass"][i]) \
                    and i != 0:
                adjustText(distribution["mass"][i - 1], distribution["abundace"][i - 1],
                           distribution["mass"][i], distribution["abundance"][i])

            else:
                plt.text(x=distribution["mass"][i],
                         y=(distribution["abundance"][i] + 2),
                         s='%0.3f' % (distribution["abundance"][i]), va='center',
                         ha='center')
                plt.ylim([0, 110])
        plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                         math.ceil(distribution["mass"][-1]) + 2))

    plt.figure(figsize=(10, 7))

    plt.subplot(1, 2, 1)
    plt.title("Isotopic distribution of oxygen")
    plotDistribution(oxygen_isoDist)
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")

    plt.subplot(1, 2, 2)
    plt.title("Isotopic distribution of sulfur")
    plotDistribution(Sulfur_isoDist)
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")

    plt.show()



#Carpon istopes and the diff between them

edb = ElementDB()
isotopes = edb.getElement( "C" ).getIsotopeDistribution().getContainer()
carbon_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()
isotopes = edb.getElement( "N" ).getIsotopeDistribution().getContainer()
nitrogen_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()
print( "Mass difference between 12C and 13C:", carbon_isotope_difference )
print( "Mass difference between 14N and N15:", nitrogen_isotope_difference )
print( "Relative deviation:", 100 * (carbon_isotope_difference -
                                     nitrogen_isotope_difference) / carbon_isotope_difference, "%" )



Methanol = EmpiricalFormula( "CH3OH" )

Water = EmpiricalFormula( "H2O" )

Ethanol = EmpiricalFormula( "CH2" ) + Methanol

print( "Ethanol chemical formula:", Ethanol.toString() )
print( "Ethanol composition:", Ethanol.getElementalComposition() )
print( "Ethanol has", Ethanol.getElementalComposition()

#aminoacids names
lys = ResidueDB().getResidue( "Lysine" )

print( lys.getName() )
print( lys.getThreeLetterCode() )
print( lys.getOneLetterCode() )
print( lys.getAverageWeight() )
print( lys.getMonoWeight() )
print( lys.getPka() )
print( lys.getFormula().toString() )


ox = ModificationsDB().getModification("Oxidation")
print(ox.getUniModAccession())
print(ox.getUniModRecordId())
print(ox.getDiffMonoMass())
print(ox.getId())
print(ox.getFullId())
print(ox.getFullName())
print(ox.getDiffFormula())

isotopes = ox.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for iso in isotopes.getContainer():
    print (iso.getMZ(), ":", iso.getIntensity())

uridine = RibonucleotideDB().getRibonucleotide( b"U" )
print(uridine.getName())
print(uridine.getCode())
print(uridine.getAvgMass())
print(uridine.getMonoMass())
print(uridine.getFormula().toString())
print(uridine.isModified())
methyladenosine = RibonucleotideDB().getRibonucleotide(b"m1A")
print(methyladenosine.getName())
print(methyladenosine.isModified())
















