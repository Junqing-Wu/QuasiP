function [data, y, J] = QuasiP(prob, data, u)
% Generator dispatch relation:
%   y = R(C,Om,p),   J = dR/d(C,Om,p).
% The selected PHI.param and PHI.method choose MP-mVCF or MT-SHT.
%QUASIP  Dispatch residual evaluation to the selected QuasiP formulation.
%
% Read data.PHI.param and data.PHI.method and call the corresponding private
% residual. The returned y is the zero-problem residual and J is its Jacobian.
%
% On input:
%   prob : COCO problem structure.
%   data : QuasiP model and PHI data.
%   u    : Unknown vector [C(:);Om(:);p(:)].
%
% On output:
%   data : Updated formulation data.
%   y    : QuasiP residual.
%   J    : Residual Jacobian when supported by the selected formulation.
switch data.PHI.param
    case{'MP'}
        switch data.PHI.method
            case{'mVCF'}
                if nargout<3 && isfield(data.PHI,'sparse_jacobian') && data.PHI.sparse_jacobian
                    isvar = data.PHI.isvar;
                    data.PHI.isvar = false;
                    [data,y] = MP_mVCF_QuasiP(prob,data,u);
                    data.PHI.isvar = isvar;
                else
                    [data,y,J] = MP_mVCF_QuasiP(prob,data,u);
                end
            case{'sCCF'}
        end
    case{'MT'}
        switch data.PHI.method
            case{'DM'}

            case{'SHT'}
                if nargout<3 && ~data.sta_opts.issta && ...
                        isfield(data.PHI,'residual_only') && data.PHI.residual_only
                    isvar = data.PHI.isvar;
                    data.PHI.isvar = false;
                    [data,y] = MT_SHT_QuasiP(prob,data,u);
                    data.PHI.isvar = isvar;
                else
                    [data,y,J] = MT_SHT_QuasiP(prob,data,u);
                end
        end
end
end

