% Copyright MIT-LICENSE 2022 Antonio POLITI<COPYRIGHT HOLDER>
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


%% run_fit_CI
% Run a fit of the model 'single_npc_logn_hill' to the data.
% Compute the confidence interval using the Chi2 countour
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

function run_fit_CI()
%% File paths
if ispc
    main_path ='P:\Code\';
else
    main_path = '/Users/apoliti/owncloud/Shotaro_Otsuka_collaboration/';
end
usestd = 1;
[indata, poi_names, outfiles, max_pois] = configure(usestd, main_path);

%% initialize model
par_ini = [0 0,  ...  % path1_reg1
    2 0, ... path1_reg2 // 2 minutes delay from one region to the other.
    0 0, ... path2_reg1
    2 0];  % path2_reg2
par_kin =  [10*rand, 50*rand, ...    % path1
    10*rand, 50*rand];      % path2
frac = [0.9,  0.3]; ...  %  path1_reg1, path1_reg2 % data is 0.90-0.94 and  0.385 - 0.617 noncore and core
    
% optimal parameters
par_offset = [0, 0];
single_npc_model = 'single_npc_logn_hill';

clear('ND')
timemod = [0:0.01:140]';
timeout = [0:0.5:140]';
ND = nuclei_dynamic(single_npc_model, par_ini,  par_kin, frac, ...
    'par_offset' , par_offset, 'usestd', usestd, ...
    'frac_lb' , [0.7 0.1], 'frac_hb', [1 0.5],...
    'par_ini_hb', 10*ones(1, 8));
% region 1 is noncore, region2 is core
reg_names = {'noncore', 'core'};

% files that contain optimal parameters. This gives the starting values and
% range for the CI computation
filename_optimalpar_frac = fullfile(main_path, ['NPCMaturation/results/singlepore/20210417_bgsub_nuclei_dynamic_usestd' num2str(usestd)], '20210417_single_npc_logn_hill_2min__frac_par.txt');
filename_optimalpar =  fullfile(main_path, ['NPCMaturation/results/singlepore/20210417_bgsub_nuclei_dynamic_usestd' num2str(usestd)], '20210417_single_npc_logn_hill_2min__par.txt');

%% Fit procedure %%%

%% fit fraction for each protein and kinetics
% table_opti_par = readtable(filename_optimalpar_frac);
% Chi2_diff= 3.84;
% for ipoi = 1:max_pois
%     data = indata.(poi_names{ipoi});
%     assert(strcmp( poi_names{ipoi}, table_opti_par(ipoi,1).POI{:}));
%     timeexp = data.core(:,1);
%     data = [data.noncore(:,2:3), data.core(:,2:3)];
%     opti_par = table2array(table_opti_par(ipoi,2:end));
%     opti_norm = opti_par(end);
%     vect_ci =  [9:12 13 14]
%     for ici = 1:length(vect_ci)
%         par_fit_id = [9:12 13 14]; % [9:12 13 14] fit kinetics and fraction for each protein
%         id_delete = find(par_fit_id == vect_ci(ici));
%         par_fit_id(id_delete) = [];
%         ref_par_value = opti_par(vect_ci(ici));
%         par_vec = opti_par;
%         
%         % Compute upper boundary
%         for iborder = 1:2
%             if iborder == 1
%                 CI_par = [ref_par_value*0.1, ref_par_value, ref_par_value];
%             else
%                 CI_par = [ref_par_value, ref_par_value, ref_par_value*5];
%             end
%             CI_par(2) = (CI_par(1)+CI_par(3))/2;
%             
%             for i = 1:3
%                 par_vec(vect_ci(ici)) = CI_par(i);
%                 ND.set_par_fromvec(par_vec);
%                 [parfit(:,i), norm(i)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
%             end
%             if abs(max(norm) - opti_norm) < Chi2_diff
%                 fit_result.(poi_names{ipoi}).par(:,iborder + (ici-1)*2) = [opti_par(1:end-1)'; vect_ci(ici)];
%                 fit_result.(poi_names{ipoi}).norm(:,iborder + (ici-1)*2) = opti_norm;
%                 continue;
%             end
%             while (abs(CI_par(2) - CI_par(1)) > 1e-3)
%                 if iborder == 1
%                 if abs(norm(2) - opti_norm) > Chi2_diff 
%                     CI_par(1) = CI_par(2);
%                 else
%                     CI_par(3) = CI_par(2);
%                 end
%                 else
%                 if abs(norm(2) - opti_norm) > Chi2_diff 
%                     CI_par(3) = CI_par(2);
%                 else
%                     CI_par(1) = CI_par(2);
%                                   end
%                 end
%                 CI_par(2) = (CI_par(1)+CI_par(3))/2;
%                 for i = 1:3
%                     par_vec(vect_ci(ici)) = CI_par(i);
%                     ND.set_par_fromvec(par_vec);
%                     [parfit(:,i), norm(i)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
%                 end
%             end
%             fit_result.(poi_names{ipoi}).par(:, iborder  + (ici-1)*2) = [parfit(:,2); vect_ci(ici)];
%             fit_result.(poi_names{ipoi}).norm(:, iborder  + (ici-1)*2) = norm(2);
%         end
%     end
% end
% 
% 
% writeout_par(outfiles{1}, ND, par_fit_id, fit_result);

%% keep fraction for all proteins constant
table_opti_par = readtable(filename_optimalpar);
Chi2_diff= 3.84;
for ipoi = 1:max_pois
    data = indata.(poi_names{ipoi});
    assert(strcmp( poi_names{ipoi}, table_opti_par(ipoi,1).POI{:}));
    timeexp = data.core(:,1);
    data = [data.noncore(:,2:3), data.core(:,2:3)];
    opti_par = table2array(table_opti_par(ipoi,2:end));
    opti_norm = opti_par(end);
    vect_ci =  [9:12]
    for ici = 1:length(vect_ci)
        par_fit_id = [9:12]; 
        id_delete = find(par_fit_id == vect_ci(ici));
        par_fit_id(id_delete) = [];
        ref_par_value = opti_par(vect_ci(ici));
        par_vec = opti_par;
        
        % Compute upper boundary
        for iborder = 1:2
            if iborder == 1
                CI_par = [ref_par_value*0.1, ref_par_value, ref_par_value];
            else
                CI_par = [ref_par_value, ref_par_value, ref_par_value*5];
            end
            CI_par(2) = (CI_par(1)+CI_par(3))/2;
            
            for i = 1:3
                par_vec(vect_ci(ici)) = CI_par(i);
                ND.set_par_fromvec(par_vec);
                [parfit(:,i), norm(i)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
            end
            if abs(max(norm) - opti_norm) < Chi2_diff
                fit_result.(poi_names{ipoi}).par(:,iborder + (ici-1)*2) = [opti_par(1:end-1)'; vect_ci(ici)];
                fit_result.(poi_names{ipoi}).norm(:,iborder + (ici-1)*2) = opti_norm;
                continue;
            end
            while (abs(CI_par(2) - CI_par(1)) > 1e-3)
                if iborder == 1
                if abs(norm(2) - opti_norm) > Chi2_diff 
                    CI_par(1) = CI_par(2);
                else
                    CI_par(3) = CI_par(2);
                end
                else
                if abs(norm(2) - opti_norm) > Chi2_diff 
                    CI_par(3) = CI_par(2);
                else
                    CI_par(1) = CI_par(2);
                                  end
                end
                CI_par(2) = (CI_par(1)+CI_par(3))/2;
                for i = 1:3
                    par_vec(vect_ci(ici)) = CI_par(i);
                    ND.set_par_fromvec(par_vec);
                    [parfit(:,i), norm(i)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
                end
            end
            fit_result.(poi_names{ipoi}).par(:, iborder  + (ici-1)*2) = [parfit(:,2); vect_ci(ici)];
            fit_result.(poi_names{ipoi}).norm(:, iborder  + (ici-1)*2) = norm(2);
        end
    end
end


writeout_par(outfiles{2}, ND, par_fit_id, fit_result);



end


function [indata, poi_names, outfiles, max_pois] = configure(usestd, main_path)
% input data


bgsub = '_bgsub'; %bgsub = ''; Uses data where for the core region the mean values of the first 3 time-points have been subtracted
indata_file = fullfile(main_path , ...
    ['NPCMaturation/expdata/livecellCalibrated/HeLa4D-core-noncore-normalized-summary-190225_norm'  bgsub '.txt']);

% output data

out_path = fullfile(main_path , ['NPCMaturation/results/singlepore/20210417' bgsub '_nuclei_dynamic_usestd' num2str(usestd)]);
if ~exist(out_path, 'dir')
    mkdir(out_path);
end
file_id = '20210417_single_npc_logn_hill_CI_';

%% First set of outfiles is when the fraction is varied for each protein, second set the fraction is identical for all proteins
outfiles = {struct('setup', [out_path '/' file_id '_frac_setup.txt'], ...
    'par', [out_path '/' file_id '_frac_par.txt'], ...
    'totkin', [out_path '/' file_id '_frac_totkin.txt'], ...
    'singlepath', [out_path '/' file_id '_frac_singlepath.txt']), ...
    struct('setup', [out_path '/' file_id '_setup.txt'], ...
    'par', [out_path '/' file_id '_par.txt'], ...
    'totkin', [out_path '/' file_id '_totkin.txt'], ...
    'singlepath', [out_path '/' file_id '_singlepath.txt'])}
[indata, poi_names] = nuclei_dynamic.loaddata(indata_file);
n_pois = length(poi_names);
max_pois = n_pois;
end


function write_singlepath(ND, single_path, time, fileout, reg_names, poi_name)
% output order
% time_min path_int path_id frac_path region POI;
% single path has dimension time, path, region
fid = fopen(fileout, 'a');
for i_reg = 1:size(single_path, 3)
    for i_path = 1:size(single_path, 2)
        for i_t = 1:size(single_path, 1)
            if i_path < 2 %holds only for 2 paths
                frac =  ND.frac(i_reg);
            else
                frac = 1 - ND.frac(i_reg); % hold only for
            end
            fprintf(fid, '%.3f\t%.3f\t%d\t%.3f\t%s\t%s\n', time(i_t), single_path(i_t, i_path, i_reg), i_path, frac,  reg_names{i_reg}, poi_name);
        end
    end
end

end


function write_totkin(tot_kin, time, fileout, reg_names, poi_name)

fid = fopen(fileout, 'a');
for i_reg = 1:size(tot_kin, 2)
    for i_t = 1:size(tot_kin, 1)
        fprintf(fid, '%.3f\t%.3f\t%s\t%s\n', time(i_t), tot_kin(i_t, i_reg), reg_names{i_reg}, poi_name);
    end
end
fclose(fid);
end

function writeout_par(outfile, ND, par_fit_id, fit_result)
poi_names = fieldnames(fit_result);

fid = fopen(outfile.setup, 'w');

fprintf(fid, 'usestd %d\nsingle_npc_model %s\nfrac_lb %s\nfrac_hb %s\npar_fit_id %s' , ...
    ND.usestd, ND.model_name, num2str(ND.frac_lb), num2str(ND.frac_hb), ...
    num2str(par_fit_id));%, ND.us
fclose(fid);

fid = fopen(outfile.par, 'w');
fprintf(fid, 'POI\tini1_p1_r1\tini2_p1_r1\tini1_p1_r2\tini2_p1_r2\tini1_p2_r1\tini2_p2_r1\tini1_p2_r2\tini2_p2_r2\t');
fprintf(fid, 'kin1_p1_r12\tkin2_p1_r12\tkin1_p2_r12\tkin2_p2_r12\t');
fprintf(fid, 'frac_p1_r1\tfrac_p1_r2\t');
fprintf(fid, 'offset_r1\toffset_r2\tpar_varied\tnorm\n');

for ipoi = 1:length(poi_names)
    
    parvec = fit_result.(poi_names{ipoi}).par;
    
    for icol = 1:size(parvec,2)
        numstr = poi_names{ipoi};
        for irow = 1:size(parvec,1)
            numstr = sprintf('%s\t%.3f', numstr, parvec(irow, icol));
        end 
        numstr = sprintf('%s\t%.3f', numstr, fit_result.(poi_names{ipoi}).norm(icol));
        fprintf(fid, '%s\n', numstr);%, ND.us
    end
end
fclose(fid);
end



