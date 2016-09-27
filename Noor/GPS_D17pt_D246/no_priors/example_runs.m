% example_run

NonlinFO_1

configCluster
noor1 = parcluster('noor1');
%ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('fukuoka')
ClusterInfo.setQueueName('defaultq')
ClusterInfo.setWallTime('600')

job = batch(noor1,'run_nice_bayesian_med', 'pool',7, 'AttachedFiles', 'extras/');

noor1 = parcluster('noor1'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('run_nogeo')
ClusterInfo.setQueueName('rh6_q2hr')
ClusterInfo.setWallTime('2:00')

job = batch(noor1, 'run_tmcmc2', 'pool',255, 'AttachedFiles', 'additional_scripts/');

noor2 = parcluster('noor2'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('run_noprior')
ClusterInfo.setQueueName('defaultq')
ClusterInfo.setWallTime('1200')

job = batch(noor2, 'run_nice_bayesian_med', 'pool',15, 'AttachedFiles', 'extras/');

