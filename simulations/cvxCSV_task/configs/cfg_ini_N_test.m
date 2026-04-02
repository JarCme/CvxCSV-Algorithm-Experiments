N = [100, 250, 500, 1000, 2500, 5000, 7500];
d = 15;
T = 10;

ini_pert = 0.5;

soi_a = 0.5;
soi_y = 0;
tau = 0;

non_soi_a = 0.5;
non_soi_y = 0;

trials = 5000;

target_SIR = 0;

a1aT_type.name = 'angle_Q_static';
a1aT_type.angle = 10;

params_fastdiva1.maxit = 10000;
params_fastdiva1.initype = 'a';
params_fastdiva1.nonln = 'sign';
params_fastdiva1.scaling = 'n';
params_fastdiva1.approach = 'u';
params_fastdiva1.numsig = 1;
params_fastdiva1.T = T;
params_fastdiva1.precond = [];

params_fastdiva2.maxit = 10000;
params_fastdiva2.initype = 'a';
params_fastdiva2.nonln = 'sign';
params_fastdiva2.scaling = 'n';
params_fastdiva2.approach = 'u';
params_fastdiva2.numsig = 1;
params_fastdiva2.T = 1;
params_fastdiva2.precond = [];

params_full.maxit = 10000;
params_full.mu = 0.001;
params_full.T = T;
params_full.nonln = 'sign';

params_bogice.maxit = 10000;
params_bogice.mu = 0.001;
params_bogice.T = T;
params_bogice.nonln = 'sign';

methods_params = dictionary();

methods_params("full") = params_full;

methods_params("bogice") = params_bogice;

methods_params("fastdiva_T10") = params_fastdiva1;
methods_params("fastdiva_T1") = params_fastdiva2;

method_list = keys(methods_params).';