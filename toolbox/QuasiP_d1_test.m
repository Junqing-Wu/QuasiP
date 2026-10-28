function prob = QuasiP_d1_test(prob, tbid, data)
% QUASIP_D1_TEST  Add Floquet test functions for d=1 solutions.
tfid = coco_get_id(tbid, 'test');
tfps = coco_get_id(tfid, {'FP' 'PD' 'TR'});
uidx = coco_get_func_data(prob, tbid, 'uidx');
data.tfid = tfid;
prob = coco_add_chart_data(prob, tfid, [], []);
prob = coco_add_func(prob, tfid, @QuasiP_d1_TF, data, 'regular', tfps, 'uidx', uidx, 'fdim', 3, 'passChart');
prob = coco_add_event(prob, 'FP', tfps{1}, 0);
prob = coco_add_event(prob, 'PD', tfps{2}, 0);
prob = coco_add_event(prob, @QuasiP_d1_evhan_TR, data, tfps{3}, 0);
end
function [data, chart, y] = QuasiP_d1_TF(prob, data, chart, u)
% QUASIP_D1_TF  Monitor functions for FP, PD and TR bifurcations.
cdata = coco_get_chart_data(chart, data.tfid);
if ~isempty(cdata) && isfield(cdata, 'la') 
    la  = cdata.la;
else
    Cii  = u(data.C_idx);
    Omii = u(data.OM_idx);
    pii  = u(data.p_idx);
    data_eval = data;
    if strcmp(data.PHI.param, 'MT') && ...
            strcmp(data.PHI.method, 'SHT')
        [data_eval, ~, ~] = MT_SHT_QuasiP(prob, data_eval, u);
    end
    Amp = PostP_amp(data_eval, Cii, Omii, 1:data.ode_opts.n, false);
    Sta = PostP_sta(data_eval, Omii, pii, Amp);
    la  = Sta.Ly(:,1);
    if data.ode_opts.e == 0 && ~isempty(la)
        [~, idx] = min(abs(la-1));
        la(idx) = [];
    end
    cdata.la  = la;
end
y = zeros(3,1);
y(1) = real(prod(la-1));
y(2) = real(prod(la+1));
if numel(la) > 1
    idx1 = [];
    idx2 = [];
    for ii = 1:numel(la)-1
        idx1 = [idx1; repmat(ii, numel(la)-ii, 1)]; 
        idx2 = [idx2; (ii+1:numel(la))']; 
    end
    y(3) = real(prod(la(idx1).*la(idx2)-1));
else
    y(3) = 1;
end
cdata.test      = y;
cdata.nunstable = sum(abs(la) > 1);
chart = coco_set_chart_data(chart, data.tfid, cdata);
end
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
function [data, cseg, msg] = QuasiP_d1_evhan_TR(prob, data, cseg, cmd, msg) 
% QUASIP_D1_EVHAN_TR  Keep unit-circle crossings of a multiplier pair.
switch cmd
    case 'init'
        if isfield(msg, 'finish') || strcmp(msg.action, 'warn')
            msg.action = 'finish';
        elseif strcmp(msg.action, 'locate')
            msg.action = 'warn';
        else
            cdata0 = coco_get_chart_data(cseg.ptlist{1}, data.tfid);
            cdata1 = coco_get_chart_data(cseg.ptlist{end}, data.tfid);
            s0 = sum(sign(abs(cdata0.la)-1));
            s1 = sum(sign(abs(cdata1.la)-1));
            if abs(s0-s1) == 4
                msg.point_type = 'TR';
                msg.action     = 'locate';
            else
                msg.action = 'finish';
            end
            msg.idx = 1;
        end
    case 'check'
        msg.action = 'add';
        msg.finish = true;
end
end