function [pc,A,L]=PC(data,npc)

% Function that computes the principal components:

Sigma = cov(data);
opts.disp = 0;
[ A, L ] = eigs(Sigma,npc,'LM',opts);
pc = data*A;

end
