function [irf, cimp, bootcimp] = FAVAR(X,Z,Transf,idvariables,k,h,L,nrepli,conf,sho,opt,interestvariables)

[imp, chi,sh] = FAVARCholImp(X, Z,idvariables,k, h);
cimp = CumImp(imp, Transf);
bootimp = FAVARCholBlockBoot(X,Z, idvariables,  k,h, L, nrepli);
bootcimp = CumImp(bootimp, Transf);
%[v, stdv, meanv] = ExplainedVariances(bootcimp,cimp, sho);
irf = confband(bootcimp, cimp, conf);
if opt == 1
DoGraphs(irf,interestvariables,sho);
end