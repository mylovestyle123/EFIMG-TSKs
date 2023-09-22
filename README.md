# Code-of-EFIMG-TSKs
The codes of EFIMG-TSKs
The codes of EFIMG-TSKs include the main EFIMG.m file and five other functional functions
(1)The CCPP as one of our experimental dataset can be available at http://archive.ics.uci.edu/dataset/294.
(2)The features in both ccpp_train and ccpp_test are normalized in this case.
(3)As a case of Algorithm 1 in our paper, the hyperparameters are directly given without grid searching for simplicity.
(4)The codes of EFIMG-TSKs on ccpp do not include the random feature selecting in the module of preprocess and the calculation of the p. For comparisons and more in-depth research, please contact us. We will provide more detailed codes after passing your assessment. 
(5)In order to make the codes runable, we adopt all the features and set the weight of each rule in each base FIMG-TSK as 1/Nc, where Nc is the number of rules in each base FIMG-TSK for simplicity.
(6)Without the random feature selecting, the optimized p and the optimal hyperparameters, the average testing RMSEs of this EFIMG-TSKs is only slightly higher than the reports in our paper, keeping it around 4.15.
