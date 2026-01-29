function X_MF = Linear_solution(varargin)
% compute the quasi-periodic solution of linear dynamical system
% 1ST ODEs: B\dot{x} = AX + E;
% 2ND ODEs: M\ddot{x} + D\dot{x} + KX = E;

if nargin == 4
    B = varargin{1};
    A = varargin{2};
    E_h = varargin{3}.harmonic_order;
    E_a = varargin{3}.amplitude;
    OM =  varargin{4};

    assert( numel(OM) == size(E_h,2) &&  numel(OM) == size(E_h,1) , ...
     '%s: the harmonic order of excitation is wrong'); 

    n = size(B,1);
    d = numel(OM);
    e_num = size(E_h,1);

    X_MF = zeros(e_num,d+n);

    for ii = 1:e_num
        X_MF(ii,1:d) = E_h(ii,:);
        X_MF(ii,d+1:end) = transpose( (1i*OM(ii)*B - A) \ transpose(E_a(ii,:)) );
    end


elseif nargin == 5
    M = varargin{1};
    D = varargin{2};
    K = varargin{3};
    E_h = varargin{4}.harmonic_order;
    E_a = varargin{4}.amplitude;
    OM =  varargin{5};

    assert( numel(OM) == size(E_h,2) &&  numel(OM) == size(E_h,1) , ...
     '%s: the harmonic order of excitation is wrong'); 

    n = size(M,1);
    d = numel(OM);
    e_num = size(E_h,1);

    X_MF = zeros(e_num,d+n);

    for ii = 1:e_num
        X_MF(ii,1:d) = E_h(ii,:);
        X_MF(ii,d+1:end) = transpose( ( - OM(ii)^2 *M +  1i*OM(ii)*D +  K) \ transpose(E_a(ii,:)) );
    end

end

end
