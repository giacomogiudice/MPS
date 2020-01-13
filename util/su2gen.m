function [sx,sy,sz,sp,sm] =  su2gen(d)
% SU2GEN Return a d-dimensional representation of the generators of SU(2).
%
% [sx,sy,sz,sp,sm] =  su2gen(d)
% Returns the spin matrices in a (d x d)-dimensional representation of
% SU(2).

% Total spin
S = (d-1)/2;
% Elements of sp, defined by sp|S,m> = sqrt(S*(S+1) - m*(m+1))|S,m+1>
ladderdiag = @(s,m) sqrt(s.*(s+1)-m.*(m+1));
sp = diag(ladderdiag(S,S -(1:d-1)),1);
sm = sp';
% Generate operators along x,y,z
sz = diag(S:-1:-S);
sx = 0.5*(sp + sm);
sy = 1i*0.5*(sm - sp);
end