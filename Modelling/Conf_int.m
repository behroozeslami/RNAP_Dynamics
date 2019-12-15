function [B,M,L,U,STD] = Conf_int(x,P,alpha)

[Max, I] = max(P);

B = x(I);

Psum =cumsum(P);

M = dot(x,P);

L = min(x(Psum > alpha));

U = min(x(Psum > (1-alpha)));

STD = sqrt(dot(x.^2,P)-M^2);