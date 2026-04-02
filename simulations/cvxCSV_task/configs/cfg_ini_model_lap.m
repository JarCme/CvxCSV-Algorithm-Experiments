N = 5000;
d = 5;
T = 10;

%ini_pert = logspace(-2,2,5)/10;

ini_pert = [0.001, 0.01, 0.1, 0.5, 1, 5]

soi_a = 0.5;
soi_y = 0;
tau = 0;

non_soi_a = 0.5;
non_soi_y = 0;

target_SIR = nan;

trials = 5000;

a1aT_type.name = 'model';
a1aT_type.angle = -1;

 params_fastdiva.maxit = 10000;
 params_fastdiva.initype = 'a';
 params_fastdiva.nonln = 'sign';
 params_fastdiva.scaling = 'n';
 params_fastdiva.approach = 'u';
 params_fastdiva.numsig = 1;
 params_fastdiva.T = T;
 params_fastdiva.precond = [];

params_full.maxit = 10000;
params_full.mu = 0.001;
params_full.T = T;
params_full.nonln = 'cggd';

params_bogice.maxit = 10000;
params_bogice.mu = 0.001;
params_bogice.T = T;
params_bogice.nonln = 'cggd';

methods_params = dictionary();

methods_params("full") = params_full;
methods_params("bogice") = params_bogice;
methods_params("fastdiva") = params_fastdiva;

method_list = keys(methods_params).';
