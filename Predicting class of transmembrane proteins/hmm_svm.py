## packages/ modules
#
import csv
import numpy as np
from sklearn import hmm
from sklearn import svm
import itertools
from matplotlib import pyplot as plt
from matplotlib import patches as pat

####### EXTRA FUNCTIONS - NO NEED FOR MODIFICATION #######

def extractAllRecords(fileName):
    """
    Extract all CSV records and return as single LIST-of-RECORDS

    Note that the first element in this list contains the
    description of the columns as defined in the CSV file
    """
    file = open (fileName, 'r')
    reader = csv.reader(file)
    results = []
    for row in reader:
        results.append(row)
    file.close()
    return (results)


def summarize (values):
    """Return the total number of occurrences of unique elements in dictionary format."""
    d = {}
    unique = []
    counts = []
    for val in values:
        if val not in unique:
            unique.append(val)

    for e in unique:
        counts.append(values.count(e))

    for i in range(len(unique)):
        key = unique[i]
        value = counts[i]
        d[key] = value

    return d

####### END EXTRA FUNCTIONS #######


## amino acid distribution preparation
#

# load in sequences.txt created in separate script
inputFile = extractAllRecords('sequences.txt')

# initialize sequences array
sequences = []

# load sequences array
for line in inputFile:
   sequences.append(line)

# parse for beta and alpha sequences, split into arrays
beta = sequences[5:34]
alpha = sequences[35:166]

# combine list of lists into list of strings
beta = list(itertools.chain(*beta))
alpha = list(itertools.chain(*alpha))

# define test array of non-beta/alpha sequences (3rd state)
# add in whatever you want here - I gave one sequence as starting point
untrainedSeq = ['RTDKYKRLKAEVEKQSKKLEKKKETITESAGRQQKKKIERQEEKLKNNNRDLSMVRMK',
                'DVGIFGSKTQPLNPSKGVEVAKEVLRGMCIARYCLDHIKDSTDFCNNTGDLRMRSSFEGG',
                'WALPFFLMGDPEAQAEVVQYNSHNRDLVILVYKKINLALIKYRMLDFQNLAEAEYVKASL',
                'KTWVNRVGIETEAPRTRKHQTDLTAIPDLTPERIMQLGGHFTGPVGSQSFSAHSTIVAPA']
totalUntrained = len(untrainedSeq)   # used later for defining test set

# initialize initial probability arrays of parameters
A = []
R = []
N = []
D = []
C = []
E = []
Q = []
G = []
H = []
I = []
L = []
K = []
M = []
F = []
P = []
S = []
T = []
W = []
Y = []
V = []

# start counter of amino acid residues at zero
scoreA = 0
scoreR = 0
scoreN = 0
scoreD = 0
scoreC = 0
scoreE = 0
scoreQ = 0
scoreG = 0
scoreH = 0
scoreI = 0
scoreL = 0
scoreK = 0
scoreM = 0
scoreF = 0
scoreP = 0
scoreS = 0
scoreT = 0
scoreW = 0
scoreY = 0
scoreV = 0

# combine arrays and change list of lists to list of strings
tmp = [beta,alpha,untrainedSeq]
tmp = list(itertools.chain(*tmp))

# for each sequence in tmp, count occurrences of amino acid residues,
# and store percentage of each residue in arrays
for e in tmp:
   scoreA = 0     # start each counter at 0 when looping through
   scoreR = 0     # a new sequence in tmp
   scoreN = 0
   scoreD = 0
   scoreC = 0
   scoreE = 0
   scoreQ = 0
   scoreG = 0
   scoreH = 0
   scoreI = 0
   scoreL = 0
   scoreK = 0
   scoreM = 0
   scoreF = 0
   scoreP = 0
   scoreS = 0
   scoreT = 0
   scoreW = 0
   scoreY = 0
   scoreV = 0
   for i in range(len(e)):    # at each residue in sequence add to
      if e[i] == 'A':         # that residue's counter
         scoreA += 1
      elif e[i] == 'R':
         scoreR += 1
      elif e[i] == 'N':
         scoreN += 1
      elif e[i] == 'D':
         scoreD += 1
      elif e[i] == 'C':
         scoreC += 1
      elif e[i] == 'E':
         scoreE += 1
      elif e[i] == 'Q':
         scoreQ += 1
      elif e[i] == 'G':
         scoreG += 1
      elif e[i] == 'H':
         scoreH += 1
      elif e[i] == 'I':
         scoreI += 1
      elif e[i] == 'L':
         scoreL += 1
      elif e[i] == 'K':
         scoreK += 1
      elif e[i] == 'M':
         scoreM += 1
      elif e[i] == 'F':
         scoreF += 1
      elif e[i] == 'P':
         scoreP += 1
      elif e[i] == 'S':
         scoreS += 1
      elif e[i] == 'T':
         scoreT += 1
      elif e[i] == 'W':
         scoreW += 1
      elif e[i] == 'Y':
         scoreY += 1
      else:
         scoreV += 1
         
   A.append(scoreA/float(len(e)))      # append the total number of 
   R.append(scoreR/float(len(e)))      # occurrences of each residue
   N.append(scoreN/float(len(e)))      # over the sequence length to
   D.append(scoreD/float(len(e)))      # previously defined initial
   C.append(scoreC/float(len(e)))      # probability arrays
   E.append(scoreE/float(len(e)))
   Q.append(scoreQ/float(len(e)))
   G.append(scoreG/float(len(e)))
   H.append(scoreH/float(len(e)))
   I.append(scoreI/float(len(e)))
   L.append(scoreL/float(len(e)))
   K.append(scoreK/float(len(e)))
   M.append(scoreM/float(len(e)))
   F.append(scoreF/float(len(e)))
   P.append(scoreP/float(len(e)))
   S.append(scoreS/float(len(e)))
   T.append(scoreT/float(len(e)))
   W.append(scoreW/float(len(e)))
   Y.append(scoreY/float(len(e)))
   V.append(scoreV/float(len(e)))


## training and testing set preparation
#

# initialize training arrays (only training with known beta and alpha sequences)
betaTrain = []
alphaTrain = []

# intialize testing arrays (testing 3 states)
betaTest = []
alphaTest = []
untrainedTest = []


# define the beta train set
for i in range(25):
   betaTrain.append([A[i],R[i],N[i],D[i],C[i],E[i],Q[i],G[i],H[i],I[i],
                  L[i],K[i],M[i],F[i],P[i],S[i],T[i],W[i],Y[i],V[i]])

print betaTrain

# define the alpha train set
for i in range(30,55):
   alphaTrain.append([A[i],R[i],N[i],D[i],C[i],E[i],Q[i],G[i],H[i],I[i],
                   L[i],K[i],M[i],F[i],P[i],S[i],T[i],W[i],Y[i],V[i]])

# define the beta test set
for i in range(25,30):
   betaTest.append([A[i],R[i],N[i],D[i],C[i],E[i],Q[i],G[i],H[i],I[i],
                     L[i],K[i],M[i],F[i],P[i],S[i],T[i],W[i],Y[i],V[i]])

# define the alpha test set
for i in range(55,60):
   alphaTest.append([A[i],R[i],N[i],D[i],C[i],E[i],Q[i],G[i],H[i],I[i],
                   L[i],K[i],M[i],F[i],P[i],S[i],T[i],W[i],Y[i],V[i]])

# define the untrained test set
for i in range((len(A)-totalUntrained),len(A)):
   untrainedTest.append([A[i],R[i],N[i],D[i],C[i],E[i],Q[i],G[i],H[i],I[i],
                         L[i],K[i],M[i],F[i],P[i],S[i],T[i],W[i],Y[i],V[i]])
         

## set up 3-state HMM
#

# define model
model = hmm.GaussianHMM(3, "tied")     # 3 states, use tied because it
                                       # can learn models from less data

# creat numpy matrices of appropriate "shape" for fitting and predicting
X = np.vstack([betaTrain,alphaTrain])  # for fitting
C = np.vstack([betaTest,alphaTest,untrainedTest])  # for predicting

# fit model with X
model.fit([X])

# predict X to see success of training
training = model.predict(X)

# predict C to see success of testing
untrainedTesting = model.predict(C)

## obtain accuracy of model
#

# put known states into arrays
betas = []
for i in range(len(betaTest)):
   betas.append(untrainedTesting[i])

alphas = []
for i in range(len(betaTest),(len(alphaTest)+len(betaTest))):
   alphas.append(untrainedTesting[i])

randoms = []
for i in range((len(alphaTest)+len(betaTest)),len(untrainedTesting)):
   randoms.append(untrainedTesting[i])

# determine number of times correct/incorrect states predicted
betaCount = summarize(betas)
alphaCount = summarize(alphas)
randomCount = summarize(randoms)

# determine number of times correct states predicted
betaMax =  max(betaCount, key=betaCount.get)
alphaMax = max(alphaCount, key=alphaCount.get)
randomMax = max(randomCount, key=randomCount.get)

# obtain total accuracy
accuracy = (betaCount[betaMax] + alphaCount[alphaMax] + randomCount[randomMax]) / float(len(untrainedTesting))

## print results
#

# print training results
print("Training: ", training)

# print testing results
print("Testing: ", untrainedTesting)

# print accuracy
print ("Accuracy: ", accuracy)

# print transition probability
print ("Transition probability matrix:")
print model.transmat_
print ""

# print means and vars for each state
print "means and vars of each hidden state"
for i in xrange(3):
    print "%dth hidden state:" % i
    print "mean = ", model.means_[i]
    print "var = ", np.diag(model.covars_[i])
    print ""

# plot testing results
plt.figure()
plt.plot(untrainedTesting)
plt.show()

### END ###

### SVM - predicting betas or alphas

Ytr = [0] * 25 + [1] * 25
Yte = [0] * 5 + [1] * 5

Csvm = np.vstack([betaTest,alphaTest])

clf = svm.SVC()
clf.fit(X,Ytr)
print(clf.predict(Csvm))
