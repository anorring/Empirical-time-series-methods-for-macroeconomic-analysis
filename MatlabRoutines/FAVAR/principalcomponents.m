%
% pc   = principalcomponents(x , npc) computes the first npc ordinary principal
%    components of x, a  matrix having series on the
%    columns. 
%
function [pc,R,D,chi] = principalcomponents(x , npc)
S = cov(x);
opts.disp = 0;
[ R, D ] = eigs(S,npc,'LM',opts);
pc = x*R;
chi = x*R*R';