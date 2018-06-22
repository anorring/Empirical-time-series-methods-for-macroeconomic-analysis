function [imp, chi,sh] = DfmCholImp(data, variables, r,k, h)
q = length(variables);
[BB, chi,rsh] = DfmRawImp(data,q,r,k,h);
[imp sh] = DfmCholIdent(BB,variables,rsh);
