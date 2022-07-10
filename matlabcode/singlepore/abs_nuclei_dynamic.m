% Copyright MIT-LICENSE 2022 Antonio POLITI
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


%% Abstract class abs_nuclei_dynamic
% function to simulate the combined dynamics observed in the different
% regions of the nuclei. Reads a model for a single_pore kinetics and
% compute distance to data.
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef abs_nuclei_dynamic <  handle
    properties
        MO;             % Model for single NPC assembly
        par_ini;        % A vector, [par_ini_path1_reg1, par_ini_path1_reg2, ..., par_ini_path2_reg1,...]
        par_kin;        % A vector [par_kin_path1, par_kin_path2, ...]
        frac;           % Fraction of different type of (n_assembly_paths-1)*n_regions
        % [frac_path1_reg1,frac_path2_reg1] or [frac_path1_reg1, frac_path1_reg2, ..., frac_path2_reg1, ...] ]
        par_offset;     % An offset to cope with variabilities in initiation part. This should not be fitted but computed directly from the data! size n_regions
        par_ini_hb;     % A vector high-boundary of initialization parameters.
        par_ini_lb;     % A vector low-boundary of initialization parameters.
        par_kin_lb;     % A vector low-boundary of kinetics parameters. 
        par_kin_hb;     % A vector low-boundary of kinetics parameters. 
        frac_lb;        % A vector fractions lowboundary
        frac_hb;        % A vector fractions highboundary 
        par_offset_lb;  % offset lowboundary
        par_offset_hb;  % offset higboundary
        usestd;     % scale data with STD. Using the complete data set would be more precise.
        model_name;           % name of model
    end
    
    properties (Access = protected)
        n_par_ini_all;  % total number of initialization parameters (path and regions)
        n_par_kin_all;  % total number of kinetic parameters (path and regions)
        n_frac_all;     % number of frac variables all is set in the constructor
    end
    
    properties(Abstract)
        n_regions;    % number of regions e.g. core, non-core
        n_assembly_paths;  % number of different assembly paths, e.g. interphase and postmitotic assembly
    end
    methods(Abstract)
        [tot_kin, kin_ini_conv_mom] = nuclei_kin(ND, time) % computes the core and non-core dynamic
    end
    
    methods
        function ND = abs_nuclei_dynamic(single_npc_model, par_ini, par_kin, frac, varargin)
            % Handle of class

  
            
        end
        
        function assign_parameters(ND, par_ini, par_kin, frac, varargin)
            p = inputParser;
            addOptional(p, 'usestd', 0, @(x) (x == 0 || x == 1));
            addRequired(p, 'par_ini', @(x) all(size(x) == [1, ND.n_par_ini_all]));
            addRequired(p, 'par_kin', @(x) all(size(x) == [1, ND.n_par_kin_all]));
            addRequired(p, 'frac', @(x) all(size(x) == [1, ND.n_frac_all]));
            addOptional(p, 'par_offset', zeros(1, ND.n_regions), @(x) all(size(x) == [1, ND.n_regions]));
            addOptional(p, 'par_ini_lb', repmat(ND.MO.get_par_ini_lb, 1, ND.n_par_ini_all/ND.MO.n_par_ini) ...
                , @(x) all(size(x) == [1, ND.n_par_ini_all]));
            addOptional(p, 'par_ini_hb', repmat(ND.MO.get_par_ini_hb, 1, ND.n_par_ini_all/ND.MO.n_par_ini) ...
                , @(x) all(size(x) == [1, ND.n_par_ini_all]));
            addOptional(p, 'par_kin_lb', repmat(ND.MO.get_par_kin_lb, 1, ND.n_par_kin_all/ND.MO.n_par_kin), ...
                @(x) all(size(x) == [1, ND.n_par_kin_all]));
            addOptional(p, 'par_kin_hb', repmat(ND.MO.get_par_kin_hb, 1, ND.n_par_kin_all/ND.MO.n_par_kin), ...
                @(x) all(size(x) == [1, ND.n_par_kin_all]));
            addOptional(p, 'frac_lb', zeros(1,  ND.n_frac_all), @(x) all(size(x) == [1,  ND.n_frac_all]));
            addOptional(p, 'frac_hb', ones(1, ND.n_frac_all), @(x) all(size(x) == [1, ND.n_frac_all]));
            addOptional(p, 'par_offset_lb', zeros(1, ND.n_regions), @(x) all(size(x) == [1, ND.n_regions]));
            addOptional(p, 'par_offset_hb', 0.5*ones(1, ND.n_regions), @(x) all(size(x) == [1, ND.n_regions]));
            parse(p, par_ini, par_kin, frac, varargin{:});
            p = p.Results;
            
            %% Initiate model values
            ND.par_ini = p.par_ini;
            ND.par_kin = p.par_kin;
            ND.par_offset = p.par_offset;
            ND.usestd = p.usestd;
            ND.frac = p.frac;
            % boundaries
            ND.frac_lb = p.frac_lb;
            ND.frac_hb = p.frac_hb;
            ND.par_ini_lb = p.par_ini_lb;
            ND.par_ini_hb =  p.par_ini_hb;
            ND.par_kin_lb =  p.par_kin_lb;
            ND.par_kin_hb = p.par_kin_hb;
            ND.par_offset_lb = p.par_offset_lb;
            ND.par_offset_hb = p.par_offset_hb;
            % set initial parameters and test
            ND.MO.setpar_ini(p.par_ini(1:ND.MO.n_par_ini));
            ND.MO.setpar_kin(p.par_kin(1:ND.MO.n_par_kin));

           

        end
        function mat = offset(ND, time)
            %% An offset to the data of the two regions. This can be a constant offset or a time dependent offset like kinetochore appearance
            for ir = 1:ND.n_regions
                mat(:, ir) = ND.par_offset(ir);
            end
        end
        
        
        function [residuals, norm, tot_kin, R2] = dist(ND, time, time_data, data)
            %% distance vector to data
            % Data is a vector containing region1, std1, region2, std2
            
            tot_kin = ND.nuclei_kin(time);
            % resample so that it matches the data points
            tot_kin = tot_kin(ismember(time, time_data), :);
            % normalize the data in case it has not reached an en
            % this is not the way to compute it. One region is normalized
            % to the value of the other. We should convert it back to what
            % it was before or use a ratio noncore/core to map this
            % difference
            residuals = [];
            SStot = 0;
            for ir = 1:ND.n_regions
                tot_kin(:,ir) = tot_kin(:,ir)/tot_kin(end, ir);
                datar_mean = data(:, 1 + (ir-1)*ND.n_regions);
                if ND.usestd
                    datar_sd = data(:, 2 + (ir-1)*ND.n_regions);
                else
                    datar_sd = ones(size(data, 1), 1);
                end
                residuals = [residuals; (tot_kin(:,ir) - datar_mean)./datar_sd];
                SStot = SStot + sum(((datar_mean - mean(datar_mean))./datar_sd).^2);
            end
            
            norm = sum(residuals.^2);
            R2 = 1- norm/SStot;
        end
        
        
        function [distvec, norm, tot_kin, R2] = dist_par(ND, par, par_idx, time, time_data, data)
            %% distance vector to data given parameter in par vector with index par_idx
            par_vec = ND.get_parvec();
            par_vec(par_idx) = par;
            ND.set_par_fromvec(par_vec);
            
            % Data is a vector containing region1, std1, region2, std2
            [distvec, norm, tot_kin, R2] = ND.dist(time, time_data, data);
            
        end
        
        function par_tot = get_parvec(ND, par_ini, par_kin,  frac, par_offset)
            %% Create a parameter vector from entries. Use default values if some entries are missing
            par_tot = [];
            if nargin < 2
                par_ini = ND.par_ini;
            end
            if nargin < 3
                par_kin = ND.par_kin;
            end
            if nargin < 4
                frac = ND.frac;
            end
            if nargin < 5
                par_offset = ND.par_offset;
            end
            
            %linearize the vector
            par_tot = [par_ini(:); par_kin(:); frac(:); par_offset(:)];
        end
        
        function par_tot = get_parvec_hb(ND)
            par_tot = ND.get_parvec(ND.par_ini_hb, ND.par_kin_hb, ND.frac_hb, ND.par_offset_hb);
        end
        
        function par_tot = get_parvec_lb(ND)
            par_tot = ND.get_parvec(ND.par_ini_lb, ND.par_kin_lb, ND.frac_lb, ND.par_offset_lb);
        end
        
        function par_tot = set_par_fromvec(ND, par_tot)
            %% assign value to parameters from vector
            a = size(ND.par_ini,1);
            b = size(ND.par_ini,2);
            for i = 1:b
                ND.par_ini(:, i) = par_tot(1 + (i-1)*a:a*i);
            end
            
            idx_offset = a*b;
            a = size(ND.par_kin,1);
            b = size(ND.par_kin,2);
            
            for i = 1:b
                ND.par_kin(:, i) = par_tot(idx_offset + 1 + (i-1)*a:idx_offset + a*i);
            end
            
            idx_offset = idx_offset + a*b;
            a = size(ND.frac,1);
            b = size(ND.frac,2);
            
            for i = 1:b
                ND.frac(:, i) = par_tot(idx_offset + 1 + (i-1)*a: idx_offset + a*i);
            end
            
            idx_offset = idx_offset + a*b;
            a = size(ND.par_offset,1);
            b = size(ND.par_offset,2);
            
            for i = 1:b
                ND.par_offset(:, i) = par_tot(idx_offset + 1 + (i-1)*a: idx_offset + a*i);
            end
        end
        
        
        
        
        function [par, norm] = lsqnonlin_fit(ND, par_idx, time, time_data, data)
            %% Run a lsq_non linear fit given the fit to par_idx (these are the idx of the linearized vector)
            par = ND.get_parvec();
            par_hb = ND.get_parvec_hb;
            par_lb = ND.get_parvec_lb;
            
            parfit = par(par_idx);
            parfit_lb = par_lb(par_idx);
            parfit_hb = par_hb(par_idx);
            options = optimset('Display', 'iter', 'MaxIter', 100);
            [parfit, norm] = lsqnonlin(@ND.dist_par, parfit, parfit_lb, parfit_hb, options, ...
                par_idx, time, time_data, data);
            
            par = ND.get_parvec();
            par(par_idx) = parfit;
            
            %R2 = 1 - norm;
            ND.set_par_fromvec(par);
        end
        
        function plot_combined_dist(ND, time)
            for ir = 1:ND.n_regions
                for ip = 1:ND.n_assembly_paths
                    ND.MO.setpar_ini(ND.par_ini(:, (ir-1)*ND.n_assembly_paths + ip));
                    ND.MO.setpar_kin(ND.par_kin(:,ip));
                    ND.MO.plot_combined_dist(time, ip + ND.n_assembly_paths*(ir-1));
                    figure(ip + ND.n_assembly_paths*(ir-1))
                    title(sprintf('region %d, path %d', ir, ip));
                end
            end
        end
        
        function plot_kinetics(ND, ifig, timemod, timeexp, data, strg_title)
            %%
            [norm, distvec, tot_kin] = ND.dist(timemod, timeexp , data);
            setupFigure(ifig, [0 0 400 300])
            plot(timeexp, data(:,1), 'o', 'Color', 'g')
            plot(timeexp, tot_kin(:,1), 'LineWidth', 2, 'Color', 'g')
            plot(timeexp, data(:,3), 'o', 'Color', 'r')
            plot(timeexp, tot_kin(:,2), 'LineWidth', 2, 'Color', 'r')
            ylim([0 1.2])
            xlim([0 120])
            xlabel('Time (min)')
            ylabel('relative fluorescence')
            title(strg_title)
        end
    end
    methods (Static)
        %% method to load the data
        function [poi_data_all, poi_names] = loaddata(file)
            indata = readtable(file);
            %%
            POIs = unique(indata.POI);
            regions = unique(indata.region);
            for POI = POIs'
                [LIA, locPOI] = ismember(indata.POI, POI{1});
                for region = regions'
                    [LIA, locRegion] = ismember(indata.region, region{1});
                    loc = locPOI.*locRegion;
                    poi_data.(region{1}) = table2array(indata(logical(loc), {'time_min', 'tot_int_mean', 'tot_int_sd'}));
                end
                poi_data_all.(POI{1}) = poi_data;
            end
            poi_names = fieldnames(poi_data_all)
        end
    end
end
