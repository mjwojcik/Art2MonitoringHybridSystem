ART-2 Neural Network & Machine Monitoring Hybrid System
=======================================================

### Matlab tools for Online Clustering


This contribution is a Matlab library implementing the ART-2 Neural Network and the Monitoring Hybrid System which is a complex system of artificial intelligence (ART-2, EMMoG, VBMoG) used to analyze and cluster multidimensional data in real-time. Originally, it was applied to the intelligent monitoring of wind turbines.



ART-2 Neural Network
--------------------

#### What it is

ART-2 network is an unsupervised neural network, based on the adaptive resonance theory (ART). It was introduced by Carpenter and Grossberg in [1].

**Art2.m** file is a Matlab class representing the model of ART-2 network. The Art2 class members represent network weights and parameters.

The Art2.m depends on an external **log4m.m** logger class which provides the ability to log network operations and states on a console or into a specified file.

#### Does it work?

The contribution contains the **Art2TestCases.m** suite which contains several unit tests based on the examples presented in section 5.3.3 of [2] as well as functional test cases. They confirm the high quality of the implementation and can show an example system running.  

#### How it works

The implementation is based on the detailed description of the ART-2 network provided by Fausett in section 5.3 of [2].

#### How to use it

* Add Art2MonitoringHybridSystem folder to Matlab search path
* Initialize the log4m logger  
 *  logger = log4m.getLogger(LOG_FILE_NAME);
* Initialize the Art2 object using structure of parameters
 * art2params = Art2.getDefaultParams(DATA_DIMENSION);
 * art2params.vigilance = 0.97; % setup network parameters, e.g. vigilance
 * art2 = Art2(art2params);
* Loop over all data points and execute the following method for each point X
	 [is_new_cluster, cluster_idx] = art2.processPoint(X);


#### External references

[1] Carpenter, G. A. & Grossberg S. (1987). *ART2: self-organization of stable category recognition codes for analog input pattern*. Applied Optics 26, 4919-4930

[2] Fausett L. (1994). *Fundamentals of Neural Networks: Architectures, Algorithms, and Applications*. Prentice-Hall, Inc., Upper Saddle River, NJ, USA.



Monitoring Hybrid System
------------------------

#### What it is

The Monitoring Hybrid System is a complex system of artificial intelligence (ART-2, EMMoG, VBMoG) used to analyze multidimensional data in real-time.
The system was designed to perform an online clustering task as a monitoring of machine states. It can be considered as an early warning tool which is able to notify about potential machine drawbacks. This approach is a data driven algorithm, which decides on a similarity of current data to the data already known by the system. In other words, the data from machine is compared to one of known states and in case when a new state is discovered, a human expert is alarmed.

**MonitoringHybridSystem.m** file is the Matlab class which is a hybrid of the following components:
* ART-2 neural network – to make general online classification
* stereographic projection and scaling unit – for data preprocessing
* Mixture of Gaussian (MoG) elements – to remember already discovered data areas
* Variational Bayesian Mixture of Gaussian (VBMoG) element – to determine how to set vigilance parameter of the ART-2 network.

**I/O: **
The system input is a stream of multidimensional data points.
The system returns numerical identifiers of recognized internal clusters (states) and it generates notifications about new clusters.

The implementation depends on the **log4m.m** logger class and Michael Chen's libraries of EMMoG and VBMoG algorithms (**GaussianMixtureLab** package).

#### Does it work?

The contribution contains the **MonitoringHybridSystemTestCases.m** suite which contains several unit tests as well as functional test cases. They confirm the high quality of the implementation and can show an example system running.

#### How it works

The core algorithm initializes ART-2 network using VBMoG unit and then manage the running network. If the number of known states is unchanged for specified period of time then it means that ART-2 network is stable and it is recorded as new MoG element. ART-2 network is used only for input points which are not matched by any saved MoG element. The algorithm can also manipulate the ART-2 vigilance parameter.
The system inform the user which ART-2 data pattern and which MoG element were used for each input data point. System produces a notification about creating new ART-2 data pattern or new MoG element, what means that a new event was observer in a processing data stream.  

For more information see Reference section.

#### How to use it

* Add Art2MonitoringHybridSystem folder to Matlab search path
* Initialize the log4m logger
 * logger = log4m.getLogger(LOG_FILE_NAME);
* Initialize the MonitoringHybridSystem object using structures of parameters
 * hybridparams = MonitoringHybridSystem.getDefaultParams(DATA_DIMENSION);
 * hybridparams.DeltaTstart = 1000; % start monitoring after 1000 data points
 * art2params = Art2.getDefaultParams(DATA_DIMENSION);
 * art2params.theta = 0.001; % setup ART-2 parameters, e.g. theta
 * hsystem = MonitoringHybridSystem(hybridparams, art2params);
* Loop over all data points and execute the following method for each point X  
  [is_new_area, is_new_cluster, area_idx, pattern_idx] = hsystem.processPoint(X);

The file runSimulation.m contains an example procedure how to run the hybrid system.

Reference
---------

**If you use this tools in some scientific research, please cite the following paper:**

**Bielecki, A. & Wójcik M. (2017).**  
***Hybrid system of ART and RBF neural networks for online clustering.***  
**Applied Soft Computing, vol. 58, 1-10, ISSN 1568-4946**

The research was supported by the National Centre for Research and Development in Poland, grant number WND-DEM-1-153/01.

The system implementation is a part of the Mateusz Wójcik’s PhD thesis - “*Hybrid neural system for intelligent monitoring of wind turbine*”.
