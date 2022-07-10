% Copyright MIT-LICENSE 2022 Antonio POLITI<COPYRIGHT HOLDER>
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% run_fit
% Run a fit of the model 'single_npc_logn_hill' to the data.
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

function run_fit()

%% File paths
if ispc
    main_path ='P:\Code\';
else
    main_path = '/Users/apoliti/ownCloud/Shotaro_Otsuka_collaboration/';
end
% input data

usestd = 1;
bgsub = '_bgsub'; 
%bgsub = '';
indata_file = fullfile(main_path , ...
    ['NPCMaturation/expdata/livecellCalibrated/HeLa4D-core-noncore-normalized-summary-190225_norm'  bgsub '.txt']);
% output data

out_path = fullfile(main_path , ['NPCMaturation/results/singlepore/20210417' bgsub '_nuclei_dynamic_usestd' num2str(usestd)]);
if ~exist(out_path, 'dir')
    mkdir(out_path);
end
delay_core_noncore = 2; % in minutes


max_rnd = 10;  % maximal repetitions of the fit with random initial condistions

file_id = ['20210417_single_npc_logn_hill_' num2str(delay_core_noncore) 'min_'];
outfiles = {struct('setup', [out_path '/' file_id '_frac_setup.txt'], ...
    'par', [out_path '/' file_id '_frac_par.txt'], ...
    'totkin', [out_path '/' file_id '_frac_totkin.txt'], ...
    'singlepath', [out_path '/' file_id '_frac_singlepath.txt']), ...
    struct('setup', [out_path '/' file_id '_setup.txt'], ...
    'par', [out_path '/' file_id '_par.txt'], ...
    'totkin', [out_path '/' file_id '_totkin.txt'], ...
    'singlepath', [out_path '/' file_id '_singlepath.txt'])}

%% load data
[indata, poi_names] = nuclei_dynamic.loaddata(indata_file);
n_pois = length(poi_names);
max_pois = n_pois;

%% initialize model
par_ini = [0, 1,  ... % path1_reg1
          0, 1, ...     % path1_reg2
          0, 1, ...     % path2_reg1
          0, 1];            % path2_reg2
par_kin =  [10*rand, 50*rand, ...    % path1
    10*rand, 50*rand];      % path2
frac = [0.9,  0.3]; ...  %  path1_reg1, path1_reg2 % data is 0.90-0.94  noncore  and  0.385 - 0.617 and core from EM postmitotic
    
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

%% Fit procedure %%%

%% fit fraction for each protein and kinetics
par_fit_id = [9:12 13 14]; % [9:12 13 14] fit kinetics and fraction for each protein
fid = fopen(outfiles{1}.totkin, 'w');
fprintf(fid, 'time_min\ttot_int\tregion\tPOI\n');
fclose(fid);
fid = fopen(outfiles{1}.singlepath, 'w');
fprintf(fid, 'time_min\tpath_int\tpath_id\tfrac_path\tregion\tPOI\n');
fclose(fid);
for ipoi = 1:max_pois
    data = indata.(poi_names{ipoi});
    timeexp = data.core(:,1);
    data = [data.noncore(:,2:3), data.core(:,2:3)];
    clear('norm')
    for rnd = 1:max_rnd
        % initial values
        par_ini = [0 0,  ...  % path1_reg1 
                   delay_core_noncore 0, ... path1_reg2 // 2 minutes delay from one region to the other.
                   0 0, ... path2_reg1
                   delay_core_noncore 0];  % path2_reg2
        par_kin =  [10*rand, 100*rand, ...    % path1
                        10*rand, 100*rand];      % path2
        frac = [0.7 + 0.3*rand,  0.1 + 0.4*rand]; %  path1_reg1, path1_reg2 % data is 0.90-0.94 and  0.385 - 0.617 noncore and core
        par_vec = ND.get_parvec(par_ini, par_kin,  frac, par_offset);
        ND.set_par_fromvec(par_vec);
        
        [parfit(:,rnd), norm(rnd)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
        
        [val, pos] = min(norm);
        
    end
    fit_result.(poi_names{ipoi}).par = parfit(:,pos);
    fit_result.(poi_names{ipoi}).norm = norm(pos);
    
    ND.set_par_fromvec(fit_result.(poi_names{ipoi}).par);
    [residuals_tmp, norm_tmp, tot_kin_tmp, R2] = ND.dist(timemod, timeexp , data);
    fit_result.(poi_names{ipoi}).R2 = R2;
    
    ND.plot_kinetics(1, timemod, timeexp, data, poi_names{ipoi});
    %saveas(1, fullfile(out_path, 'figures', [file_id '_frac_' poi_names{ipoi} '.png']));
    [tot_kin, kin_ini_conv_mom] = ND.nuclei_kin(timeout);
    write_totkin(tot_kin, timeout, outfiles{1}.totkin, reg_names, poi_names{ipoi});
    write_singlepath(ND, kin_ini_conv_mom, timeout, outfiles{1}.singlepath, reg_names, poi_names{ipoi});
end
writeout_par(outfiles{1}, ND, par_fit_id, fit_result)

%% find mean fractions and repeat the fit
frac_all = [];
for ipoi = 1:max_pois
    frac_all = [frac_all, fit_result.(poi_names{ipoi}).par(13:14)];
end

frac_all =  mean(frac_all, 2)';
init_all = [];
for ipoi = 1:max_pois
    init_all = [init_all, fit_result.(poi_names{ipoi}).par(1:8)];
end

init_all =  mean(init_all, 2)';
%% fit only single kinetics with fixed fractions
par_fit_id = [9:12];
fid = fopen(outfiles{2}.totkin, 'w');
fprintf(fid, 'time_min\ttot_int\tregion\tPOI\n');
fclose(fid);
fid = fopen(outfiles{2}.singlepath, 'w');
fprintf(fid, 'time_min\tpath_int\tpath_id\tfrac_path\tregion\tPOI\n');
fclose(fid);
for ipoi = 1:max_pois
    data = indata.(poi_names{ipoi});
    timeexp = data.core(:,1);
    data = [data.noncore(:,2:3), data.core(:,2:3)];
    clear('norm')
    for rnd = 1:max_rnd
        par_ini = init_all
        par_kin =  [10*rand, 100*rand, ...    % path1
            10*rand, 100*rand];      % path2
        par_vec = ND.get_parvec(par_ini, par_kin,  frac, par_offset);
        ND.set_par_fromvec(par_vec);
        % first round is done with default parameters other round the mean
        % is used.
        frac =  frac_all;  %  path1_reg1, path1_reg2 % data is 0.90-0.94 and  0.385 - 0.617 noncore and core
        
        data = indata.(poi_names{ipoi});
        timeexp = data.core(:,1);
        data = [data.noncore(:,2:3), data.core(:,2:3)];
        [parfit(:,rnd), norm(rnd)] = ND.lsqnonlin_fit(par_fit_id, timemod, timeexp, data);
    end
    [val, pos] = min(norm);
    fit_result.(poi_names{ipoi}).par = parfit(:,pos);
    fit_result.(poi_names{ipoi}).norm = norm(pos);
    ND.set_par_fromvec(fit_result.(poi_names{ipoi}).par);
    [residuals, norm, tot_kin, R2] = ND.dist(timemod, timeexp , data);
    fit_result.(poi_names{ipoi}).R2 = R2;
    ND.plot_kinetics(1, timemod, timeexp, data, poi_names{ipoi});

    %saveas(1, fullfile(out_path, 'figures', [file_id '_' poi_names{ipoi} '.png']))
    [tot_kin, kin_ini_conv_mom] = ND.nuclei_kin(timeout);
    write_totkin(tot_kin, timeout, outfiles{2}.totkin, reg_names, poi_names{ipoi});
    write_singlepath(ND, kin_ini_conv_mom, timeout, outfiles{2}.singlepath, reg_names, poi_names{ipoi});

end
writeout_par(outfiles{2}, ND, par_fit_id, fit_result)
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
fprintf(fid, 'offset_r1\toffset_r2\tnorm\tR2\n');

for ipoi = 1:length(poi_names)
    parvec = fit_result.(poi_names{ipoi}).par;
    numstr = poi_names{ipoi} ;
    for i = 1:length(parvec)
        numstr = sprintf('%s\t%.3f', numstr, parvec(i));
    end
    numstr = sprintf('%s\t%.3f\t%.3f', numstr, fit_result.(poi_names{ipoi}).norm, fit_result.(poi_names{ipoi}).R2);
    
    fprintf(fid, '%s\n', numstr);%, ND.us
end
fclose(fid);
end



