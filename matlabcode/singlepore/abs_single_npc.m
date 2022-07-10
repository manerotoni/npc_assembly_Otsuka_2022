% Copyright MIT-LICENSE 2022 Antonio POLITI
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Abstract class abs_single_npc
% A single pore similation
% There is an initiation dynamics and a single pore kinetics
% par_ini:
% par_kin: are the parameters for the kinetics
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef abs_single_npc < handle
    properties (Access = protected, Abstract = true)
        par_kin_lb; % low boundary for kinetics
        par_kin_hb; % high boundary for kinetics
        par_ini_lb; % low boundary for initiation
        par_ini_hb; % high boundary for initiation
        name;       % a name of the model
    end
    
    properties (Constant, Abstract  = true)
        n_par_kin;  % number of parameters kinetics model
        n_par_ini;  % number of parameters initiation model
    end
    
    properties( Access = protected)
        par_ini;    % are the parameters for the initiation
        par_kin;    % parameters for the kinetics
    end
    methods
        
        function MO = abs_single_npc(varargin)
           if nargin < 1
                par_ini = zeros(1, MO.n_par_ini);
                %display('zero vector for par_ini')
            else
                par_ini = varargin{1};
            end
            if nargin < 2
                par_kin = zeros(1, MO.n_par_kin);
                %display('zero vector for par_kin')
            else
                par_kin = varargin{2};
            end
            %% constructor
            MO.par_ini = par_ini;
            MO.par_kin = par_kin;
            MO.assert_par_kin(par_kin);
            MO.assert_par_ini(par_ini);
        end
        
        function bool = assert_par_ini(MO, par)
            %% check number of parameters initialization function
            bool = (length(par) == MO.n_par_ini);
            if bool == 0
                display('par_ini length = %d.' ,MO.n_par_ini);
            end
        end
        
        function bool = assert_par_kin(MO, par)
            %% check number of parameters initialization function
            bool = (length(par) == MO.n_par_kin);
            if bool == 0
                display('par_kin length = %d.' ,MO.n_par_kin);
            end
        end
        
        function setpar_ini(MO, par)
            %% set parameter for initiation model
            assert(MO.assert_par_ini(par), 'Number of parameters or range for par_ini is not correct');
            MO.par_ini = par;
        end
        
        function setpar_kin(MO, par)
            %% set parameter for kinetic model
            assert(MO.assert_par_kin(par), 'Number of parameters or range for par_kin is not correct');
            MO.par_kin = par;
        end
        
        function X = kin_ini_conv_mom(MO, time, imom)
            %% Convolve kinetics momnent function with initialization
            % imom is the moment of the single pore kinetics.
            % This needs to be derived properly for each model
            % This is X(t) = int_0^\infty g(t') \int_0^\infty f(x, t - t') x^imom d x dt'
            ini_pdf = MO.ini_pdf(time);
            if length(ini_pdf) == 1
                X = 0.*(time < ini_pdf) + ...
                    MO.kin_mom_pdf(time-ini_pdf, imom)/(time(2)-time(1)).*(time >= ini_pdf);
                return
            end
            if (sum(ini_pdf(2:end)) == 0) % no initialization function
                X = MO.kin_mom_pdf(time, imom)/(time(2)-time(1));
            else
                X = conv(MO.ini_pdf(time), MO.kin_mom_pdf(time, imom));
            end
            X = X(1:length(time));
        end
        
        function plot_ini_cdf(MO, time, ifig)
            %% Plot the initialization cumulative density function
            figure(ifig);
            clf;
            plot(time, MO.ini_cdf(time));
            xlabel('Time (s)');
            ylabel('Initiation cdf');
        end
        
        function plot_ini_pdf(MO, time, ifig)
            %% Plot the initialization probability density function
            figure(ifig)
            clf;
            plot(time, MO.ini_pdf(time));
            xlabel('Time (s)');
            ylabel('Initiation pdf');
        end
        
        function plot_kin_mom_pdf(MO, time, imom, ifig)
            %% Plot the initialization probability density function kinetics given the moment
            figure(ifig);
            clf;
            plot(time, MO.kin_mom_pdf(time, imom));
            xlabel('Time (s)');
            ylabel(sprintf('NUP kinietics moment %d', imom));
        end
        
        function plot_combined_dist(MO, time, ifig)
            %% Plot the probability density function of kinetics, initialization and convolved kinetics
            dt = time(2)-time(1);
            figure(ifig);
            clf;
            subplot(3,1,1)
            plot(time, MO.kin_mom_pdf(time, 1));
            xlabel('Time (min)');
            ylabel('NUP kinetics 1st moment');
            
            subplot(3,1,2)
            plot(time, MO.ini_pdf(time));
            xlabel('Time (min)');
            ylabel('Initiation pdf');
            
            subplot(3,1,3)
            kin_ini_conv = MO.kin_ini_conv_mom(time, 1);
            plot(time, kin_ini_conv(1:length(time))*dt);
            xlabel('Time (min)');
            ylabel('Average NPCs occupacy');
        end
        
        function par = get_par_ini(MO)
            %% get parameter of initialization density function
            par = MO.par_ini;
        end
        
        function par = get_par_kin(MO)
            %% get parameters of kinetics probabilities
            par = MO.par_kin;
        end
        function par = get_par_ini_lb(MO)
            %% get parameter of initialization density function
            par = MO.par_ini_lb;
        end
        
        function par = get_par_kin_lb(MO)
            %% get parameters of kinetics probabilities
            par = MO.par_kin_lb;
        end
        
        function par = get_par_ini_hb(MO)
            %% get parameter of initialization density function
            par = MO.par_ini_hb;
        end
        
        function par = get_par_kin_hb(MO)
            %% get parameters of kinetics probabilities
            par = MO.par_kin_hb;
        end
    end
    
    methods(Abstract)
        vec = ini_cdf(MO, time) % int_0^t g(t) dt  cumulative distribution function for the initiation. Returns a vector
        vec = ini_pdf(MO, time) % g(t) probability density function for the initiation. Returns a vector
        mat = kin_pdf(MO, time) % f(x,t) probability density function for the kinetics where x is occupancy. Should return a matrix
        vec = kin_mom_pdf(MO, time, imom) % Moments of the distribution int_0^\infty f(x,t) x^imom. F(t) = kin_mom_pdf(time, 1)
    end
    
    
end