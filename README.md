# Modelling-adaptation-Ultraslow-oscillations
Masters thesis for NTNU


## **1. Simulate LIF model:**

go to folder named second model (multiple neuron)

use ***network_model.m***, that calls the functions LIF.m and Gen_Current.m among others.


## **2. Simulate single neuron HH model:**

go to folder named Hodgkin-Huxley_model

use ***singleHH.m***, that calls the function excitatory_HH.m

Swap with inhibitory_HH.m to simulate an inhibitory neuron. 


## **3. Simulate network w. sequences:**

In folder named HH_network, use ***Sompolinsky_Network.m*** This will need to call most of the functions in the folder.
Simulating sequences over 550 ms takes approx. 15-20 minutes with a network of 6 of each neurons. 
This can be adjusted by changing the step size dt.


