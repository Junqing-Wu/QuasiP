% 2ND ODEs: M*\ddot{x} + D*\dot{x} + K*x + n(x,\dot{x}) = e  
% 1ST ODEs: B*\dot{X} = A*X + N(X) + E  ----> X = [] 
% B = [D M; M 0]; A = [-K 0; 0 M]; N = [-n(X); 0]; E = [e; 0];

%  \ddot{x} - mu (1-x^2)*\dot{x} + Om0^2*x +alp*x^3 = \sum_{i=1}^d epsd*cos(OMd*t)
%  M =1; + D = -mu; K = Om0^2; n(x,\dot{x}) = mu*x^2*\dot{x} + alp*x^3; e = \sum_{i=1}^d epsd*cos(OMd*t);
%  B = [-mu 1; 1 0]; A = [-Om0^2 0; 0 1]
%  N = [-(mu*X1^2*X2 + alp*X1^3) ; 0]
%  E = [e = \sum_{i=1}^d epsd*cos(OMd*t); 0]

clearvars;
close all;
clc;
%%

RHO = [1/5^0.5;1/11^0.5];
OM0 = [10;RHO*10];
pname = {'mu'};
p0 = 0.2;
x0 = zeros(2,1);
tau0 = zeros(3,1);
%%

B = VD_B(p0);
F0 = VD_F(tau0,x0,p0,'','');
X_MF0 = Linear_solution(B,F0.A,F0.E_MF,OM0);

prob = coco_prob();
ode_opts = struct('d',3, 'e',3 );
prob = coco_set(prob, 'cont','PtMX',1000);
prob = coco_set(prob, 'QuasiP', 'order', 1, 'DOF', 2, 'ode_opts', ode_opts); 


PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
    'type',    'mVCF'           , ... 
    'd',       ode_opts.d,...
    'e',       ode_opts.e,...    
    'method',  { {'HB','HB','HB'}   }, ...   
    'isCO1',   0,...                     
    'k',       { {1:1:1,1:1:1,1:1:1} }, ...        
    'U',       { {1*2+1,1*2+1,1*2+1} }, ...    
    'S',       { {2^5,2^5,2^5}     }, ...
    'S_wave',  { {2^8,2^8,2^8}     });
prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 



% PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
%     'type',    'mVCF'           , ... 
%     'd',       ode_opts.d,...
%     'e',       ode_opts.e,...    
%     'method',  { {'HB','FD'}   }, ...   
%     'isCO1',   0,...                     
%     'k',       { {1:1:5,-2:1:1} }, ...        
%     'U',       { {5*2+1,2^5} }, ...    
%     'S',       { {2^5,2^5}     }, ...
%     'S_wave',  { {2^8,2^8}     });
% prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 


% PHI_INITIAL = struct( 'parametrization_tori', 'MP', ...
%     'type',    'mVCF'           , ... 
%     'd',       ode_opts.d,...
%     'e',       ode_opts.e,...    
%     'method',  { {'HB',' CO'}   }, ...   
%     'isCO1',   1,...                     
%     'k',       { {1:1:5,[4,1]} }, ...        
%     'U',       { {5*2+1,2^5} }, ...    
%     'S',       { {2^5,2^5}     }, ...
%     'S_wave',  { {2^8,2^8}     });
% prob = coco_set(prob, 'QuasiP', 'PHI_INITIAL', PHI_INITIAL); 


% par_opts.parallel = true;
% par_opts.ncores = 6;
% prob = coco_set(prob, 'QuasiP', 'par_opts', par_opts);

funcs = {@VD_B, @VD_B_P, @VD_F, @VD_F_X, @VD_F_P};

prob = ode_isol2QuasiP(prob, '', funcs{:}, X_MF0, OM0, pname, p0);

[data, uidx] = coco_get_func_data(prob, 'QuasiP', 'data', 'uidx');

% fc2_funcs = {@fc_2};
% prob = coco_add_func(prob, 'OM2_1', fc2_funcs{:}, [], 'zero','uidx',uidx( [data.OM_idx(1),data.OM_idx(2)] ) );
% 
% fc3_funcs = {@fc_3};
% prob = coco_add_func(prob, 'OM3_1', fc3_funcs{:}, [], 'zero','uidx',uidx( [data.OM_idx(1),data.OM_idx(3)] ) );

prob = coco_add_func(prob, 'OM', @fc, [], 'zero', 'uidx', uidx(data.OM_idx) );


prob = coco_add_pars(prob, 'OM1', uidx(data.OM_idx(1)), 'Om1');

coco(prob, 'qp', [], 1, { 'Om1' 'mu'}, [0.2 10]);



function [data, y, J] = fc(prob, data, u)
RHO_2 = 1/5^0.5;
RHO_3 = 1/11^0.5;

Om1 = u(1);
Om2 = u(2);
Om3 = u(3);

y = [RHO_2* Om1 - Om2; RHO_3* Om1 - Om3];
J = [RHO_2 , -1,  0
     RHO_3 ,  0, -1];
end


function [data, y, J] = fc_2(prob, data, u)
RHO_2 = 1/5^0.5;
Om1 = u(1);
Omi = u(2);
y = RHO_2* Om1 - Omi;
J = [RHO_2 , -1];
end

function [data, y,J] = fc_3(prob, data, u)
RHO_3 = 1/11^0.5;
Om1 = u(1);
Omi = u(2);
y = RHO_3* Om1 - Omi;
J = [RHO_3 , -1];
end



