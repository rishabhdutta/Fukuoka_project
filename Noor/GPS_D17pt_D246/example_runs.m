% example_run

clear all; close all; clc

NonlinFO_1


configCluster
noor1 = parcluster('noor1');
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('fukuoka')
ClusterInfo.setQueueName('rh6_q24hr')
ClusterInfo.setWallTime('10:00')

job = batch(noor1,'run_nice_bayesian_med', 'pool',7, 'AttachedFiles', 'extras/');