% Copyright MIT-LICENSE 2022 Antonio POLITI
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Class nuclei_dynamic
% function to simulate the combined dynamics observed in the different
% regions of the nuclei. Reads a model for a single_pore kinetics and
% compute distance to data.
% In this model the initialization is independent for the path, but
% identical for each path. 
% So all in all there are 2 set of initialization parameters one per each
% pathway.
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef nuclei_dynamic_limitinit <  abs_nuclei_dynamic
    properties
        n_regions = 2;    % number of regions e.g. core, non-core
        n_assembly_paths = 2;  % number of different assembly paths, e.g. interphase and postmitotic assembly
    end
    methods
        function ND = nuclei_dynamic_limitinit(single_npc_model, par_ini, par_kin, frac, varargin)
             %% Constructor
            % single_npc_model: a string corresponding to model of npc assembly
            % par_ini: a vector of initial distribution parameters. The
            % number of parameters is par_ini or of model*ND.n_regions*ND.n_assembly_paths
            % par_kin: matrix of kinetic parameters
            % frac: fraction of postmito, inter
            % optional arguments
            %   par_offset: offset parameters
            %   par_ini_lb: low boundary initialization
            %   par_ini_hb: high boundary initialization
            %   par_kin_lb: low boundary kinetics
            %   par_kin_hb: high boundary kinetics
            
            %% consistency check
              %% consistency check
            ND = ND@abs_nuclei_dynamic();
            % Handle of class
            try
                clsHdl = str2func(single_npc_model);
            catch
                error(sprintf('No class named %s', single_npc_model));
            end
            %initialize model to access properties
            ND.MO = clsHdl();
            ND.n_par_ini_all = ND.MO.n_par_ini*ND.n_assembly_paths;
            ND.n_par_kin_all = ND.MO.n_par_kin*ND.n_assembly_paths;
            ND.n_frac_all = (ND.n_assembly_paths-1)*ND.n_regions;
            ND.assign_parameters(par_ini, par_kin, frac, varargin{:});

        end
        
 
        function [tot_kin, kin_ini_conv_mom] = nuclei_kin(ND, time)
            %% Combine kinetics of two assembly pathways to create nuclei kinetics
            
            % consistency check column wise time
            if size(time,2)> 1
                time = time';
            end
            dt = time(2)-time(1);
            % compute offset function
            offset = ND.offset(time);
            % kinetics convoluted with initiation probability
            kin_ini_conv_mom = zeros(length(time), ND.n_assembly_paths, ND.n_regions);
            tot_kin = zeros(length(time), ND.n_regions);
            %%
            for ip = 1:ND.n_assembly_paths
                for ir = 1:ND.n_regions 
                    idx_start = 1 + (ip-1)*ND.MO.n_par_ini;
                    idx_end = ND.MO.n_par_ini + (ip-1)*ND.MO.n_par_ini;
                    ND.MO.setpar_ini(ND.par_ini(idx_start:idx_end));
                    
                    idx_start = 1 + (ip-1)*ND.MO.n_par_kin;
                    idx_end = ND.MO.n_par_kin + (ip-1)*ND.MO.n_par_kin;
                    ND.MO.setpar_kin(ND.par_kin(idx_start:idx_end));
                    kin_ini_conv_mom(:, ip, ir) = ND.MO.kin_ini_conv_mom(time, 1)*dt;
                end
            end
            %%
            for ir = 1:ND.n_regions
                tot_kin(:, ir) = offset(:,ir);
                for ip = 1:ND.n_assembly_paths
                    idx_frac = ip + (ND.n_assembly_paths - 1)*(ir-1);
                    idx_start = (ND.n_assembly_paths - 1)*(ir-1) + 1;
                    idx_end = (ND.n_assembly_paths -1)  + (ND.n_assembly_paths - 1)*(ir-1);
                    if ip == ND.n_assembly_paths
                        frac = 1 - sum(ND.frac(idx_start:idx_end));
                    else
                        frac = ND.frac(idx_frac);
                    end
                    tot_kin(:, ir) = frac*kin_ini_conv_mom(:, ip, ir) + tot_kin(:, ir);
                end
            end
        end
    end 
        

    
    methods (Static)            
        %% method to load the data

        function test_get_set_parvec()
            %% test if linearization of parameter vectors is correct or not
            par_ini = [3, 1,  ...  % path1_reg1_reg2
                       10, 1];     % path2_reg1_reg2  
            par_kin =  [2, 5, ...    % path1
                        2, 30];      % path2
            frac = [0.9,  0.3];   %  path1_reg1, path1_reg2 % data is 0.90-0.94 and  0.385 - 0.617
            par_offset = [0, 0];
            ND = nuclei_dynamic_limitinit( 'single_npc_logn_hill', par_ini, par_kin, frac, 'par_offset', par_offset);
            frac = [1, 0.1];
            par_ini = [1, 2,  ...  % path1_reg1_reg2
                       4, 3];     % path2_reg1_reg2  
            
            par_vec = ND.get_parvec(par_ini, par_kin, frac, par_offset);
            ND.set_par_fromvec(par_vec);
            assert(all(par_ini == ND.par_ini));
            assert(all(par_kin == ND.par_kin));
            assert(all(frac == ND.frac));
            assert(all(par_offset == ND.par_offset));
        end
        
        function test_model_fit()
            par_ini = [5, 0.5,  ...  % path1_reg1_reg2
                       5, 0.5];     % path2_reg1_reg2  
            par_kin =  [2, 5, ...    % path1
                        2, 30];      % path2
            frac = [0.9,  0.3]; ...  %  path1_reg1, path1_reg2 % data is 0.90-0.94 and  0.385 - 0.617
            par_offset = [0, 0];

            single_npc_model = 'single_npc_logn_hill';
            ND = nuclei_dynamic_limitinit(single_npc_model, par_ini,  par_kin, frac, 'par_offset' , par_offset, 'usestd', 1);
            indata_file = 'P:\Code\NPCMaturation\expdata\livecellCalibrated\HeLa4D-core-noncore-normalized-summary-190225_norm.txt';
            if ~exist(indata_file, 'file')
                indata_file = uigetfile();
            end
            indata = ND.loaddata(indata_file);
            time = [0:0.01:140];
            timeexp = indata.TPR.core(:,1);
            data = [indata.TPR.noncore(:,2:3), indata.TPR.core(:,2:3)];
            [norm, distvec, tot_kin] = ND.dist(time, timeexp , data);
            
            if ND.usestd 
                norm_man = [(data(:,1) - tot_kin(:,1))./data(:,2); (data(:,3) - tot_kin(:,2))./data(:,4)];
            else
                norm_man = [(data(:,1) - tot_kin(:,1)); (data(:,3) - tot_kin(:,2))];
            end
            
            norm_man  = sum(norm_man.^2);
            par_vec = ND.get_parvec(par_ini, par_kin, frac, par_offset);
            ND.lsqnonlin_fit([1 3 5:8], time, timeexp, data);
            [norm, distvec, tot_kin] = ND.dist(time, timeexp , data);
            figure(2)
            clf
            hold
            plot(timeexp, data(:,[1 3]), 'o')
            plot(timeexp, tot_kin)
            
        end
    end
end
