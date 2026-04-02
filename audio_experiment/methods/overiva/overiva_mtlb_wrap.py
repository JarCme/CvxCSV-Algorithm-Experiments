import importlib
import numpy as np
import overiva
# importlib.reload(overiva)

# X = np.asarray(X)

# print(X.shape)
# print(X.dtype)
# print(n_src)                                 
# print(n_iter)
# print(proj_back)
# print(W0.shape)
# print(model)
# print(init_eig)
# print(return_filters)

if W0 is not None:
    # print(W0.shape)
    if len(W0.shape)==2:
        W0 = np.expand_dims(W0,2);

if return_filters:
    [S, W] = overiva.overiva(np.array(X),n_src=n_src,n_iter=n_iter,proj_back=proj_back,model=model,init_eig=init_eig,return_filters=return_filters,W0=W0)
    W = np.array(W);
else:
    S      = overiva.overiva(np.array(X),n_src=n_src,n_iter=n_iter,proj_back=proj_back,model=model,init_eig=init_eig,return_filters=return_filters,W0=W0)
S = np.array(S);
