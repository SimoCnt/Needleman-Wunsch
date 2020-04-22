"""
Es.1:
    Seq1: TFDERILGVQTYWAECLA
    Seq2: QTFWECIKGDNATY

Es.2:
    Seq1: LTGARDWEDIPLWTDWDIEQESDFKTRAFGTANCHK
    Seq2: TGIPLWTDWDLEQESDNSCNTDHYTREWGTMNAHKA
"""


import numpy as np
import Bio.SubsMat.MatrixInfo as bio
import matplotlib.pyplot as plt
import seaborn as sns
import sys


class my_dictionary(dict):

    def __init__(self):
        self = dict()

    def add(self, key, value):
        self[key] = value



################################
# Dizionario dei valori Blosum #
################################
# Contiene i valori di tutte le possibili coppie di amminoacidi per costruire la
# matrice di sostituzione
X = bio.blosum62
#X



######################################################
# Creazione matrice Blosum (matrice di sostituzione) #
######################################################
print("\nSequenza 1: ")
seq1 = input()
print("\nSequenza 2: ")
seq2 = input()
while(True):
    print("\nGap opening penalty: ")
    delta = float(input())
    print("\nGap extension penalty: ")
    gamma = float(input())
    if(delta >= gamma):
        break
    else:
        print("\nAttenzione! Il valore Gap extension penalty deve essere minore di Gap opening penalty")

bl = np.zeros((len(seq1), len(seq2)))

for i in range(bl.shape[0]):
    for j in range(bl.shape[1]):
        if(seq1[i], seq2[j]) in X:
            bl[i][j] = X[(seq1[i], seq2[j])]
        else:
            bl[i][j] = X[(seq2[j], seq1[i])]

bl = np.transpose(bl)
#bl



##################################
# Normalizzazione matrice blosum #
##################################
m = bl.min()
for i in range(len(seq2)):
    for j in range(len(seq1)):
        bl[i][j] = bl[i][j] - m

#bl



############################################################
# Creazione matrice dei punteggi (matrice di allineamento) #
############################################################
pos = my_dictionary() # Creo un dizionario che mi conserva la provenienza di ciascun punteggio di ogni cella
mScore = np.zeros((len(seq1), len(seq2)))
mScore = np.transpose(mScore)

mScore[0][0] = bl[0][0]
pos.add((0, 0), ('end', mScore[0][0]))


# Celle prima riga
for j in range(1, mScore.shape[1]):
    lenGap = 0

    # Controllo se la cella [0][j-1] presenta gap e calcolo la lunghezza di quest'ultimo
    if(pos[(0, j-1)][0] == 'sx'):
        z = j-1
        while(pos[(0, z)][0] == 'sx'):
            lenGap = lenGap + 1
            z = z - 1

    # Identifico il valore massimo tra le celle confrontate
    if((bl[0][j]) > (mScore[0][j-1] - delta - (gamma * (lenGap - 1)))):
        mScore[0][j] = bl[0][j]
        pos.add((0, j), ('end', mScore[0][j]))
    else:
        mScore[0][j] = mScore[0][j-1] - delta - (gamma * (lenGap - 1))
        pos.add((0, j), ('sx', mScore[0][j]))


# Celle prima colonna
for i in range(1, mScore.shape[0]):
    lenGap = 0

    # Controllo se la cella [i-1][0] presenta gap e calcolo la lunghezza di quest'ultimo
    if(pos[(i-1, 0)][0] == 'up'):
        z = i-1
        while(pos[(z, 0)][0] == 'up'):
            lenGap = lenGap + 1
            z = z - 1

    # Identifico il valore massimo tra le celle confrontate
    if((bl[i][0]) > (mScore[i-1][0] - delta - (gamma * (lenGap - 1)))):
        mScore[i][0] = bl[i][0]
        pos.add((i, 0), ('end', mScore[i][0]))
    else:
        mScore[i][0] = mScore[i-1][0] - delta - (gamma * (lenGap - 1))
        pos.add((i, 0), ('up', mScore[i][0]))


# Celle centrali
for i in range(1, mScore.shape[0]):
    for j in range(1, mScore.shape[1]):
        lenGapSx = 0
        lenGapUp = 0

        if(pos[(i, j-1)][0] == 'sx'):
            z = j-1
            while(pos[(i, z)][0] == 'sx'):
                lenGapSx = lenGapSx + 1
                z = z - 1

        if(pos[(i-1, j)][0] == 'up'):
            z = i-1
            while(pos[(z, j)][0] == 'up'):
                lenGapUp = lenGapUp + 1
                z = z - 1

        diag = mScore[i-1][j-1] + bl[i][j]
        sx = mScore[i][j-1] - delta - (gamma * (lenGapSx - 1))
        up = mScore[i-1][j] - delta - (gamma * (lenGapUp - 1))

        mv = max(diag, sx, up) 
        mScore[i][j] = mv

        if(mv == diag):
            pos.add((i,j), ('diag', mv))
        elif(mv == sx):
            pos.add((i,j), ('sx', mv))
        else:
            pos.add((i,j), ('up', mv))

#pos

"""
plt.figure(figsize=(9,9))
sns.heatmap(mScore, annot=True, cmap="Blues_r", linewidths=.5, square=True, xticklabels=seq1, yticklabels=seq2)
"""


# Inizializzo i due indici alla posizione del punteggio massimo della matrice mScore
max = 0
for i in range(mScore.shape[0]):
    for j in range(mScore.shape[1]):
        if(mScore[i][j] >= max):
            max = mScore[i][j]
            i_max = i
            j_max = j


# Punteggio dell'allineamento
punteggio = mScore[i_max][j_max]



#############################
# Ricerca percorso migliore #
#############################
percorso = []

i = i_max
j = j_max

ax = []
ay = []


if(pos[(i, j)][0] == 'diag'):
    ax.append(i+1)
    ay.append(j+1)
elif(pos[(i, j)][0] == 'sx'):
    ax.append(i)
    ay.append(j+1)
else:
    ax.append(i+1)
    ay.append(j)


while True:
    if(pos[(i, j)][0] == 'diag'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        i = i - 1
        j = j - 1
    elif(pos[(i, j)][0] == 'sx'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        j = j - 1
    elif(pos[(i, j)][0] == 'up'):
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        i = i - 1
    else: # Termina quando pos[(i, j)][0] == 'end'
        percorso.append([[seq1[j], seq2[i]], pos[(i, j)][0]])
        ax.append(i)
        ay.append(j)
        break

percorso.reverse()
#percorso



################################
# Allineamento tra le stringhe #
################################
newSeq1 = []
newSeq2 = []
confr = []
index1 = 0
index2 = 0

while(seq1[index1] != percorso[0][0][0]):
    newSeq2.append(' ')
    index1 = index1 + 1

while(seq2[index2] != percorso[0][0][1]):
    newSeq1.append(' ')
    index2 = index2 + 1

if(index1 > 0):
    for i in range(index1):
        newSeq1.append(seq1[i])

if(index2 > 0):
    for i in range(index2):
        newSeq2.append(seq2[i])

for i in range(len(percorso)):
    if(percorso[i][1] == 'diag'):
        newSeq1.append(percorso[i][0][0])
        newSeq2.append(percorso[i][0][1])
    elif(percorso[i][1] == 'up'):
        newSeq1.append('-')
        newSeq2.append(percorso[i][0][1])
    elif(percorso[i][1] == 'sx'):
        newSeq1.append(percorso[i][0][0])
        newSeq2.append('-')
    else:
        newSeq1.append(percorso[i][0][0])
        newSeq2.append(percorso[i][0][1])

i = i_max + 1
j = j_max + 1

while(i < len(seq2)):
    newSeq2.append(seq2[i])
    i = i + 1

while(j < len(seq1)):
    newSeq1.append(seq1[j])
    j = j + 1

if(len(newSeq1) < len(newSeq2)):
    count = len(newSeq1)
else:
    count = len(newSeq2)

for i in range(count):
    if(newSeq1[i] == newSeq2[i] and newSeq1[i] != '-' and newSeq2 != '-'):
        confr.append('|')
    else:
        confr.append(' ')



#####################################
# Stampa dell'allineamento migliore #
#####################################
print('\n')
print(*newSeq1, sep=" ")
print(*confr, sep=" ")
print(*newSeq2, sep=" ")
print('\n')



#####################
# Plot del Percorso #
#####################
#punteggio = round(punteggio)

plt.figure(figsize=(12,12))
plt.title("Percorso migliore (punteggio: {})".format(punteggio))
sns.heatmap(mScore, annot=True, cmap="Blues_r", linewidths=.5, square=True, xticklabels=seq1, yticklabels=seq2)
plt.plot(ay, ax, 'r-')
plt.show()
