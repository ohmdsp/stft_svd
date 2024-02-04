function mag=cabs2(x)

% Computes sum of complex vector x element squared magnitudes
%
%     mag=cabs2(x)

mag = sum(real(x).^2 + imag(x).^2);
%
