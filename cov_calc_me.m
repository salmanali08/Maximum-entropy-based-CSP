function [ S_i_me ] = cov_calc_me( S_i,Sp )

% COV_CALC_ME Summary of this function goes here
% This code is the implementation of an algorithm proposed in 
%%
% C. E. Thomaz, D. F. Gillies and R. Q. Feitosa, "A new covariance estimate for Bayesian classifiers in biometric recognition," 
% in IEEE Transactions on Circuits and Systems for Video Technology, vol. 14, no. 2, pp. 214-223, Feb. 2004.
% doi: 10.1109/TCSVT.2003.821984
%%

[phi_i_me,lmds]=eig((S_i+Sp));

[lmds,Idxs] = sort(diag(lmds),'descend'); 
phi_i_me= phi_i_me(:,Idxs);

lmds_i=diag(phi_i_me'*S_i*phi_i_me);
lmds_i_p=diag(phi_i_me'*Sp*phi_i_me);
lmds_i=[lmds_i,lmds_i_p];
lmds_i_me=max(lmds_i,[],2);
Z_i_me=diag(lmds_i_me);
S_i_me=phi_i_me*Z_i_me*phi_i_me';

end

