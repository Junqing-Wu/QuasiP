function [data, res] = PostP_Amp_Sta_bddat(prob, data, command, varargin)
%POSTP_AMP_STA_BDDAT  Export Processor results to COCO branch data.
%
% This callback is registered as a COCO bddat slot. During initialization it
% declares amplitude and stability column names. During chart evaluation it
% extracts the current QuasiP variables, calls PostP_amp and PostP_sta, and
% returns the values stored on the continuation branch.
%
% On input:
%   prob      : COCO continuation problem structure.
%   data      : QuasiP slot data.
%   command   : COCO slot command such as ''init'' or ''data''.
%   varargin  : Command-specific chart or slot arguments.
%
% On output:
%   data      : Updated slot data.
%   res       : Cell array of branch-data values.

res = {};
my_tbid = data.my_tbid;
prefix_tbid = data.prefix;
switch command
    case 'init'
        if ~data.amp_opts.isamp && ~data.sta_opts.issta
            temp = {'C', 'Om', 'p'};
        elseif data.amp_opts.isamp && ~data.sta_opts.issta
            if data.amp_opts.isuamp
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'HT_u' 'T_x' 'T_u' 'MF_x' 'MF_f' 't'};
            else
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'T_x' 'MF_x' 'MF_f' 't'};
            end
        elseif data.sta_opts.issta
            if data.amp_opts.isuamp
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'HT_u' 'T_x' 'T_u' 'MF_x' 'MF_f' 't' 'Psi' 'Ly' 'Stab'};
            else
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'T_x' 'MF_x' 'MF_f' 't' 'Psi' 'Ly' 'Stab'};
            end
        end
        for j = 1:length(temp)
            res{end+1} = [prefix_tbid, temp{j}];
        end
    case 'data'
        res = {};
        num_QuasiP = numel(prob.slots.save_full);
        chart = varargin{1};
        [data_ii, uidx_ii] = coco_get_func_data(prob, my_tbid, 'data', 'uidx');
        Cii   = chart.x(uidx_ii(data_ii.C_idx));
        Omii  = chart.x(uidx_ii(data_ii.OM_idx));
        pii   = chart.x(uidx_ii(data_ii.p_idx));
        PHI   = data_ii.PHI;
        n     = PHI.n;
        if isequal(prob.cont.atlas, @atlas_1d.create)
            iscont =  true;
            if data.isp
                cont_para  = chart.x(num_QuasiP*data_ii.p_idx(end)+1);
            else
                cont_para  = chart.x(num_QuasiP*data_ii.OM_idx(end)+1);
            end
        else
            iscont = false;
        end
        figure(11111)
        if ~data.amp_opts.isamp && ~data.sta_opts.issta
            res = {reshape(Cii,[],n), Omii', pii'};
            if iscont
                plot(cont_para,norm(Cii),'ko');
                xlabel('cont para');
                ylabel('C');
                hold on
            end
        elseif data.amp_opts.isamp && ~data.sta_opts.issta
            outdof = data.amp_opts.outdof;
            Amp    = PostP_amp(data_ii,Cii,Omii,outdof,data_ii.amp_opts.isuamp);
            if ~data_ii.amp_opts.isht
                Amp.HT_x = [];Amp.HT_u = [];
            end
            if data_ii.amp_opts.isuamp
                res = {reshape(Cii,[],n),Omii',pii',Amp.Maxht_x,Amp.maxx,Amp.minx,...
                    Amp.HT_x, Amp.HT_u,Amp.T_x, Amp.T_u,...
                    Amp.MF_x, Amp.MF_f, Amp.t};
            else
                res = {reshape(Cii,[],n),Omii',pii',Amp.Maxht_x,Amp.maxx,Amp.minx,...
                    Amp.HT_x,Amp.T_x,Amp.MF_x, Amp.MF_f, Amp.t};
            end
            if iscont
                plot(cont_para,Amp.Maxht_x(1),'ko');
                xlabel('cont para');
                ylabel('||x_1||')
                hold on
            end
        elseif data.sta_opts.issta
            outdof = 1:n;
            Amp    = PostP_amp(data_ii,Cii,Omii,outdof,data_ii.amp_opts.isuamp);
            Sta    = PostP_sta(data_ii,Omii,pii,Amp);
            % Keep the transition history at TR for ode_TR2QuasiP restarts.
            % Floquet multipliers and stability are stored at every point.
            if data_ii.ode_opts.d==1 && isfield(data_ii.sta_opts,'store_all_psi') ...
                    && ~data_ii.sta_opts.store_all_psi ...
                    && ~contains(chart.pt_type,'TR')
                Sta.Psi = [];
            end
            if ~data_ii.amp_opts.isht
               Amp.HT_x = [];Amp.HT_u = [];
            else
               Amp.HT_x = Amp.HT_x(:,data_ii.amp_opts.outdof);
               if data_ii.amp_opts.isuamp
                  Amp.HT_u = Amp.HT_u(:,data_ii.amp_opts.outdof);
               end
            end
            if data_ii.amp_opts.isuamp
                res = {reshape(Cii,[],n),Omii',pii',Amp.Maxht_x,Amp.maxx,Amp.minx,...
                    Amp.HT_x, Amp.HT_u,Amp.T_x, Amp.T_u,...
                    Amp.MF_x, Amp.MF_f, Amp.t,Sta.Psi, Sta.Ly, Sta.Stab};
            else
                res = {reshape(Cii,[],n),Omii',pii',Amp.Maxht_x,Amp.maxx,Amp.minx,...
                    Amp.HT_x,Amp.T_x,Amp.MF_x, Amp.MF_f, Amp.t,Sta.Psi, Sta.Ly, Sta.Stab};
            end
            if iscont
                if Sta.Stab(1)
                    plot(cont_para,Amp.Maxht_x(1),'ko');
                else
                    plot(cont_para,Amp.Maxht_x(1),'r^');
                end
                xlabel('cont para');
                ylabel('||x_1||')
                hold on
            end
        end
end
end
