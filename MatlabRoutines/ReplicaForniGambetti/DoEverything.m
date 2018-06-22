function [irf,v,stdv,meanv,chi,sh,cimp,bootcimp,imp,bootimp] = ...
    DoEverything(X, idvariables, r, k, L,Transf, nrepli,conf,h,sho,opt,interestvariables)
if nargin <= 10
    opt = 0;
end
if nargin <= 9
sho = 3;
end
if nargin <= 8
    h = 48;
end
if nargin <= 7
    conf = 0.8;
end
if nargin <= 6
    nrepli = 100;
end
[imp, chi,sh] = DfmCholImp(X, idvariables, r,k, h);
cimp = CumImp(imp, Transf);
bootimp = DfmCholBlockBoot(X, idvariables, r, k,h, L, nrepli);
bootcimp = CumImp(bootimp, Transf);
[v, stdv, meanv] = ExplainedVariances(bootcimp,cimp, sho);
irf = confband(bootcimp, cimp, conf);
if opt == 1
DoGraphs(irf,interestvariables,sho);
end