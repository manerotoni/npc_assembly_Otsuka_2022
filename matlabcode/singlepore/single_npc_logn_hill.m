% Copyright MIT-LICENSE 2022 Antonio POLITI<COPYRIGHT HOLDER>
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


%% class single_npc_gamma_hill
% A single npc model with deterministic hill kind of kinetics and
% lognormal initialization probability density function
% Hill function: t^n/(t^n + K^n)
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef single_npc_logn_hill < abs_single_npc
    % the convoluted function is similar to a hill-function for certain
    % parameter combinations
    % par_ini(1) : Mean of distribution of logn
    % par_ini(2) : std of distribution of logn
    % par_kin(1) : the hill-coefficient. t^par_kin(1)./(time^par_kin(1) + par_kin(2)^par_kin(1));
    % par_kin(2) : the inflexion point
    properties (Access = protected)
        % Set of boundaries for the parameters these are default and are
        % Low and high boundary
        par_kin_lb = [0.5, 1];
        par_kin_hb = [10, 200];
        par_ini_lb = [1, 1];
        par_ini_hb = [200, 200];
        name = 'single_npc_logn_hill';
    end
    
    properties (Constant)
        n_par_kin = 2;
        n_par_ini = 2;
    end
    methods
        
        function MO = single_npc_logn_hill(varargin)            
            display('single_npc_logn_hill: Log normal initiation, hill sigmoidal kinetics')
            MO = MO@abs_single_npc(varargin{:});
        end
        
        function vec = ini_pdf(MO, time)
            %% g(t) probability density function for the initiation. Returns a vector
            % convert parmaters to entries for log-normal distribution
            mu = log((MO.par_ini(1)^2)/sqrt(MO.par_ini(2)+MO.par_ini(1)^2));
            sigma = sqrt(log(MO.par_ini(2)/(MO.par_ini(1)^2)+1));
            % check function
            % [M, V]= lognstat(mu, sigma); M == par_ini(1), V == par_ini(2)
            if MO.par_ini(2) == 0 && MO.par_ini(1) > 0 % all pores initiate at the same time
                vec = MO.par_ini(1);
                return
            end
            if MO.par_ini(1) == 0
                vec = time*0;
            else
                vec = lognpdf(time, mu, sigma); % this is faster
            end
        end
        
        function vec = ini_cdf(MO, time)
            %% g(t) probability density function for the initiation. Returns a vector
            
            mu = log((MO.par_ini(1)^2)/sqrt(MO.par_ini(2)+MO.par_ini(1)^2));
            sigma = sqrt(log(MO.par_ini(2)/(MO.par_ini(1)^2)+1));
            % check function
            % [M, V]= lognstat(mu, sigma); M == par_ini(1), V == par_ini(2)
            vec = logncdf(time, mu, sigma);
        end
        
        
        function mat = kin_pdf(MO, time)
            %% NOT IMPLEMENTED f(x,t) probability density function for the kinetics where x is occupancy. Should return a matrix
            %% of probability of being x at time t
            display('Distribution of NUP kinetics has not been defined')
            mat = [];
        end
        
        function vec = kin_mom_pdf(MO, time, imom)
            %% Moments of the distribution int_0^\infty f(x,t) x^imom. F(t) = kin_mom_pdf(time, 1)
            if imom == 1
                vec = time.^MO.par_kin(1)./(time.^MO.par_kin(1) + MO.par_kin(2).^MO.par_kin(1));
            else
                display('Only first moment of kinetics has been defined');
            end
        end
    end
    methods (Static)
        function test()
            %% A test function that can be directly called without instantiating the class
            par_ini = [10; 5];
            par_kin = [2; 40];
            time = [0:0.1:200]';
            amodel = single_npc_logn_hill(par_ini, par_kin);
            amodel.plot_combined_dist(time, 1);
            % find a function that correspond to occupancy
            figure(2)
            clf
            hold
            plot(time, amodel.kin_ini_conv_mom(time, 1)*0.1, 'r')
           [val, pos] = min(abs(0.5 - amodel.kin_ini_conv_mom(time, 1)*0.1));
            par_kin = [3; time(pos)];
            plot(time, time.^par_kin(1)./(time.^par_kin(1) + par_kin(2).^par_kin(1)), 'b');
            
 
        end
    end
    
end