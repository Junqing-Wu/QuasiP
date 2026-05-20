function [data, res] = PostP_forQP(prob, data, command, varargin)
res = {};
my_tbid = data.my_tbid;
prefix_tbid = data.prefix;
switch command
    case 'init'
        if ~data.amp_opts.isamp && ~data.stable_opts.issta
            temp = {'C', 'Om', 'p'};
        elseif data.amp_opts.isamp && ~data.stable_opts.issta
            if data.amp_opts.isuamp
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'HT_u' 'T_x' 'T_u' 'MF_x' 'MF_f' 't'};
            else
                temp   = {'C' 'Om' 'p' 'Max_x' 'max.x' 'min.x' 'HT_x' 'T_x' 'MF_x' 'MF_f' 't'};
            end
        elseif data.stable_opts.issta
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
        chart = varargin{1};
        [data_ii, uidx_ii] = coco_get_func_data(prob, my_tbid, 'data', 'uidx');
        Cii   = chart.x(uidx_ii(data_ii.C_idx));
        Omii  = chart.x(uidx_ii(data_ii.OM_idx));
        pii   = chart.x(uidx_ii(data_ii.p_idx));
        PHI   = data_ii.PHI;
        n     = PHI.n;
        if ~data.amp_opts.isamp && ~data.stable_opts.issta
            res = {reshape(Cii,[],n), Omii', pii'};
        elseif data.amp_opts.isamp && ~data.stable_opts.issta
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
        elseif data.stable_opts.issta
            outdof = 1:n;
            Amp    = PostP_amp(data_ii,Cii,Omii,outdof,data_ii.amp_opts.isuamp);
            Sta    = PostP_sta(data_ii,Omii,pii,Amp);
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
        end
end
end