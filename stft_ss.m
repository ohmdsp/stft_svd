function tfr = stft_ss(num_freqs,window,T,p,factor,sigs,x)

%load helicopter.mat
% num_freqs = 256;
% window = 0;
% T = 1/1024;
% p = 10;
% factor = 2;
% sigs = 5;
% x = helicopter;

% modified 17 June 1995 to revert to LP if 0 order occurs

% Short-time Fourier transform using original data + extrapolated
% data obtained using the signal subspace approach by truncating
% the number of eigenvectors of the SVD from the covariance method
% data matrix that are used, then determining the forward and
% backward linear prediction parameters from the truncated SVD
%   factor = data extension factor (2 means twice data length, 3
%            means three times data length)
%   sigs = maximum number of signal subspace eigenvectors to save

show = 0;
order_lp = 0;

[m,n] = size(x);
if m == 1
   y = x.';
   N = n;
elseif n == 1
   y = x;
   N = m;
else
   disp('Data passed to STFT_LP is a matrix!')
end

[Uf,Sf,Vf] = svd(toeplitz(y(p:N-1),y(p:-1:1)),0);
[Ub,Sb,Vb] = svd(hankel(y(2:N-p+1),y(N-p+1:N)),0);
for k = sigs:-1:0
   order = k;
   if order == 0, break, end
   a_f = -Vf(:,1:k)*diag(1./diag(Sf(1:k,1:k)))*Uf(:,1:k)'*y(p+1:N);
   if show == 0
      root_for = abs(roots([1 a_f.']));
   else
      root_for = abs(roots([1 a_f.']));
   end
   if all( root_for < 1.02 )
      a_b = -Vb(:,1:k)*diag(1./diag(Sb(1:k,1:k)))*Ub(:,1:k)'*y(1:N-p);
      if all( abs(roots([1 a_b'])) < 1.02 ), break, end
   end
end

if order == 0
   for k = sigs:-1:0
      order_lp = k;
      if order_lp == 0, break, end
      [rho_f,a_f,rho_b,a_b] = covar(k,y);
      if show == 0
         root_for = abs(roots([1 a_f.']));
      else
         root_for = abs(roots([1 a_f.']));
      end
      if all( root_for < 1.02 )
         if all( abs(roots([1 a_b'])) < 1.02 ), break, end
      end
   end
end

extension = (factor-1)*N/2;
if order_lp ~= 0
   pp = order_lp;
else
   pp = p;
end
if (order ~= 0) | (order_lp ~= 0)
   for k = 1:extension
      y = [y; -a_f.'*y(N+k-1:-1:N+k-pp)];
   end
   for k= 1:extension
      y = [-a_b.'*y(1:pp); y];
   end
end

tfr = perogram_marple(num_freqs,window,0,length(y),T,y);

disp(['order = ',int2str(order),'   order_lp = ',int2str(order_lp)])

if show == 1
   if (order == 0) & (order_lp == 0)
      y = [zeros(extension,1); y; zeros(extension,1)];
   end
   subplot(2,1,1)
   plot((-extension+1:N+extension),real(y),'y')
   hold on
   plot((1:N),real(x),'r')
   hold off
   subplot(2,1,2)
   plot((-extension+1:N+extension),imag(y),'y')
   hold on
   plot((1:N),imag(x),'r')
   hold off
   pause
end
%
