
import os,sys
import numpy as np
from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.descriptor.cutoffs import Cosine,Polynomial
from pyamff.mlModels.pytorchNN import NeuralNetwork
from pyamff.pyamff import PyAMFF
from amp.utilities import make_filename, hash_images, check_images, Logger
#from amp.model.pytorchNN_ import NeuralNetwork
import torch


#Define NN parameters
weights = {}
bias = {}
params = []

#Define log file
log = Logger(make_filename('test', '-log.txt'))

#Read in images
images = hash_images(images='train.traj', log=log)
check_images(images, forces=True)

parallel = {'envcommand': None, 'cores':1}

#Define fingerprints
Gs = {"H": [{"type":"G2", "element":"H", "eta":0.05, "Rs":0.},
            {"type":"G2", "element":"Pd", "eta":160., "Rs":0.},
           ],
      "Pd": [{"type":"G2", "element":"H", "eta":0.05, "Rs":0.},
             {"type":"G2", "element":"Pd", "eta":400., "Rs":6.0},
             ]}
descriptor=Gaussian(Gs=Gs,
                    cutoff=Cosine(6.0)
                    )
#Calculate fingerpints
descriptor.calculate_fingerprints(
   images=images,
   parallel=parallel,
   log=log,
   calculate_derivatives=True)

var = -1.
l1_weight = [[-0.06245638, -0.02177071], [ 0.01106360, -0.02641876]]
l1_bias  = [-0.15927899, 0.08015668]
l2_weight = [[-0.01062437, -0.00506942],[0.00031429, -0.12331749]]
l2_bias  = [0.18475361, 0.17465138]
l3_weight = [[-0.10928684, 0.02496535]]
l3_bias  = [0.12033310]

l1_weight1 = [[-0.12125415, -0.09574964], [-0.12175528,-0.01690174]]
l1_bias1   = [-0.00525357, 0.16651967]
l2_weight1 = [[0.14360507, -0.03717202],[0.18517374,-0.01645058]]
l2_bias1   = [0.03095685, -0.14551263]
l3_weight1 = [[-0.12728526, 0.03061981]]
l3_bias1   = [-0.00860174]

params = [l1_weight, l1_bias, l2_weight, l2_bias, l3_weight, l3_bias,
          l1_weight1, l1_bias1, l2_weight1, l2_bias1, l3_weight1, l3_bias1]

#params = None
#Define the NN model
model=NeuralNetwork(
          hiddenlayers=(2, 2),
          nFPs = {'H': 2, 'Pd':2},
          energy_coefficient=1.0,
          force_coefficient=0.1,
          learningrate = 0.00001,
          params =params
          )

#Define the pyamff calculator
calc = PyAMFF(model = model,optimizer='SD', 
              debug=None,
              weight_decay=0.1)

#run the calculation
calc.fit(trainingimages=images, descriptor=descriptor, maxEpochs=10, log=log)
#calc.fitSingleModel(trainingimages=images, descriptor=descriptor, maxEpochs=5000, log=log)
#print(var, loss)
