N = 5000;
d = 5;
T = 10;

ini_pert = 0.3;

soi_a = 1;
soi_y = 0:0.1:0.9;
tau = 0;

trials = 500;

a1aT_type.name = 'model_eye';
a1aT_type.angle = -1;

% params_fastdiva.maxit = 10000;
% params_fastdiva.initype = 'a';
% params_fastdiva.nonln = 'sign';
% params_fastdiva.scaling = 'n';
% params_fastdiva.approach = 'u';
% params_fastdiva.numsig = 1;
% params_fastdiva.T = T;
% params_fastdiva.precond = [];

params_full.maxit = 50000;
params_full.mu = 0.001;
params_full.T = T;
params_full.nonln = 'cggd';

params_bogice.maxit = 50000;
params_bogice.mu = 0.001;
params_bogice.T = T;
params_bogice.nonln = 'cggd';

methods_params = dictionary();

methods_params("full") = params_full;
methods_params("bogice") = params_bogice;

method_list = keys(methods_params).';
