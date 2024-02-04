function [rho_f,a_f,rho_b,a_b]=covar(p,x)

% Covariance least squares autoregressive parameter estimation algorithm using
% QR-decomposition type of solution.
%
%     [rho_f,a_f,rho_b,a_b]=covar(p,x)                               
%
% p     -- order of linear prediction/autoregressive filter
% x     -- vector of data samples
% rho_f -- least squares estimate of forward linear prediction variance
% a_f   -- vector of forward linear prediction/autoregressive parameters
% rho_b -- least squares estimate of backward linear prediction variance
% a_b   -- vector of backward linear prediction/autoregressive parameters

%*****************  Initialization [ eq (8.D.42) ]  *****************

n = length(x);
if  2*p+1 > n, error('Order too high; will make solution singular.'), end
a_f = [];
a_b = [];
r1 = cabs2(x(2:n-1));
r2 = cabs2(x(1));
r3 = cabs2(x(n));
if p <= 0
    rho_f = (r1 + r2 + r3)/n;
    rho_b = rho_f;
    return
end
rho_f = r1 + r3;
rho_b = r1 + r2;
r1 = 1/(rho_b + r3);
c = conj(x(n))*r1;
d = conj(x(1))*r1;
ef = x;                                                 % eq (8.D.1)
eb = x;                                                 % eq (8.D.2)
ec = c*x;
ed = d*x;

%*************************  Main Recursion  *************************

for k=1:p

    if (rho_f <= 0) | (rho_b <= 0)
       error('A prediction squared error was less than or equal to zero.')
    end
    gam = 1 - real(ec(n-k+1));                          % eq (8.D.29)
    del = 1 - real(ed(1));                              % eq (8.D.28)
    if (gam <= 0) | (gam > 1) | (del <= 0) | (del > 1)
       error('GAM or DEL gain factor not in range 0 to 1.')
    end

    % computation for k-th order reflection coefficients
    [eff,ef_k] = splitoff(ef,'top');
    [ebb,eb_n] = splitoff(eb,'bottom');
    delta = ebb'*eff;                                   % eq (8.D.32)
    k_f = -delta/rho_b;                                 % eq (8.D.19)
    k_b = -conj(delta)/rho_f;                           % eq (8.D.22)

    % order updates for squared prediction errors  rho_f and rho_b
    rho_f = rho_f*(1 - real(k_f*k_b));                  % eq (8.D.20)
    rho_b = rho_b*(1 - real(k_f*k_b));                  % eq (8.D.23)

    % order updates for linear prediction parameter arrays  a_f and a_b
    temp = a_f;
    a_f = [a_f; 0] + k_f*[flipud(a_b);  1];             % eq (8.D.18,19)
    a_b = [a_b; 0] + k_b*[flipud(temp); 1];             % eq (8.D.21,22)

    % check if maximum order has been reached
    if k == p
        rho_f = rho_f/(n-p);
        rho_b = rho_b/(n-p);
        return
    end

    % order updates for prediction error arrays  ef and eb
    eb = ebb + k_b*eff;                                 % eq (8.D.24)
    ef = eff + k_f*ebb;                                 % eq (8.D.24)

    % coefficients for next set of updates
    c1 = ec(1);
    c2 = c1/del;
    c3 = conj(c1)/gam;

    % time updates for gain arrays  c' and d"
    temp = c;
    c = c + c2*d;                                       % eq (8.D.38)
    d = d + c3*temp;                                    % eq (8.D.39)

    % time updates for ec' and ed"
    temp = ec;
    ec = ec + c2*ed;                                    % eq (8.D.40)
    ed = ed + c3*temp;                                  % eq (8.D.41)
    [ecc,ec_k] = splitoff(ec,'top');
    [edd,ed_n] = splitoff(ed,'bottom');

    if (rho_f <= 0) | (rho_b <= 0)
       error('A prediction squared error was less than or equal to zero.')
    end
    gam = 1 - real(ecc(n-k));                           % eq (8.D.29)
    del = 1 - real(edd(1));                             % eq (8.D.28)
    if (gam <= 0) | (gam > 1) | (del <= 0) | (del > 1)
        error('GAM or DEL gain factor not in range 0 to 1.')
    end

    % coefficients for next set of updates
    c1 = ef(1);
    c2 = eb(n-k);
    c3 = conj(c2)/rho_b;
    c4 = conj(c1)/rho_f;
    c5 = c1/del;
    c6 = c2/gam;

    % order updates for c and d; time updates for a_f' and a_b"
    temp = flipud(a_b);
    a_b = a_b + c6*flipud(c);                           % eq (8.D.35)
    c = [c; 0] + c3*[temp; 1];                          % eq (8.D.25,27)
    temp=a_f;
    a_f = a_f + c5*d;                                   % eq (8.D.33)
    d = [0; d] + c4*[1; temp];                          % eq (8.D.25,26)

    % time updates for rho_f' and rho_b"
    rho_f = rho_f - real(c5*conj(c1));                  % eq (8.D.34)
    rho_b = rho_b - real(c6*conj(c2));                  % eq (8.D.36)

    % order updates for ec and ed; time updates for ef' and eb"
    ec = ecc + c3*eb;                                   % eq (8.D.30)
    eb = eb + c6*ecc;                                   % eq (8.D.37)
    ed = edd + c4*ef;                                   % eq (8.D.31)
    ef = ef + c5*edd;                                   % eq (8.D.37)
end
% 
