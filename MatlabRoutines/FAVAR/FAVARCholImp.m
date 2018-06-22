function [imp, chi,sh] = FAVARCholImp(X, Z, variables,k, h)
[BB, chi,rsh] = FAVARRaw(X, Z, k,h);
[imp sh] = FAVARCholIdent(BB,variables,rsh);
