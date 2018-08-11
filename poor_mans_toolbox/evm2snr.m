function [ snr ] = evm2snr( evm )
%function [ snr ] = evm2snr( evm )
%
%   Takes EVM in RMS and converts it to decibel SNR.
%
% (c) 2015 ryan@guerra.rocks

snr = 10*log10(1./x.^2);

end
