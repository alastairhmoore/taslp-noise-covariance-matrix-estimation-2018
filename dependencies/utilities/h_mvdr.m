function h = h_mvdr(Rvv,steering_vec_ifreq)
% H_MVDR returns the frequency-domain mvdr beamformer.
%
% INPUTS:
% Rvv: M*M noise covariance matrix, where M is the number of channels.
% steering_vec_ifreq: M*1 steering vector of one frequency.
%
% OUTPUT:
% h:  M*1 frequency-domain mvdr beamformer
%
% REFERENCES: 
% [1] Benesty, Jacob, Jingdong Chen, and Yiteng Huang.
% Microphone Array Signal Processing. Springer Science & Business Media,
% 2008.
% [2] Pan, Chao, Jingdong Chen, and J. Benesty. Performance Study
% of the MVDR Beamformer as a Function of the Source Incidence Angle.
% IEEE/ACM Transactions on Audio, Speech, and Language Processing 22, no. 1
% (2014): 67-79.
% 
% Wei Xue, Imperial College, Feb 14, 2017
% rvvinv = invmat(Rvv);
inv_product = Rvv\steering_vec_ifreq;
h = inv_product/(steering_vec_ifreq'*inv_product);