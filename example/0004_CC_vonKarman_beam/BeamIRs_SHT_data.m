function [bd,runs] = BeamIRs_SHT_data()
% Read the continued SHT branch, including a saved restart if needed.
run = 'beamIRs_qp_sht';
bd = coco_bd_read(run);
runs = repmat({run},size(bd,1)-1,1);
resume = 'beamIRs_qp_sht_resume';
if isempty(coco_bd_labs(bd,'SHT3')) && isfile(fullfile('data',resume,'bd.mat'))
    tail = coco_bd_read(resume);
    om = coco_bd_col(bd,'Om');
    omTail = coco_bd_col(tail,'Om');
    assert(isequal(bd(1,:),tail(1,:)) && norm(om(end,:)-omTail(1,:),inf)<1e-7, ...
        'The SHT restart does not match the end of the original run.');
    bd = [bd;tail(3:end,:)];
    runs = [runs;repmat({resume},size(tail,1)-2,1)];
end
end
