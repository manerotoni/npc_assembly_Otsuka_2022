% Copyright MIT-LICENSE 2022 Antonio POLITI<COPYRIGHT HOLDER>
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% class single_npc_logn_logistic
% A single npc model with deterministic logistic function and logn
% initialization probability denisity function
% The logistic function is the same as used in
% Dultz, Ellenberg 2010, JCB.
% This function is problematic due to the zero value only t=-infty!
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef single_npc_logn_logistic < abs_single_npc
    % par_ini(1) : Mean of distribution of logn
    % par_ini(2) : std of distribution of logn
    properties
        % Set of boundaries for the parameters these are default and are
        % changed with fit
        par_kin_lb = [0.5, 1];
        par_kin_hb = [10, 200];
        par_ini_lb = [1, 2];
        par_ini_hb = [200, 200];
        name = 'single_npc_logn_logistic';
    end
    properties (Constant)
        n_par_kin = 2;
        n_par_ini = 2;
    end
    methods
        function MO = single_npc_logn_logistic(varargin)
            display('single_npc_logn_logistic: Log normal initiation, logistic kinetics')
            MO = MO@abs_single_npc(varargin{:});
        end
        

        
        function vec = ini_pdf(MO, time)
            %% g(t) probability density function for the initiation. Returns a vector
            % convert parmaters to entries for log-normal distribution
            mu = log((MO.par_ini(1)^2)/sqrt(MO.par_ini(2)+MO.par_ini(1)^2));
            sigma = sqrt(log(MO.par_ini(2)/(MO.par_ini(1)^2)+1));
            % [M, V]= lognstat(mu, sigma); M == par_ini(1), V == par_ini(2)
            
            pd = makedist('Lognormal',mu, sigma);  % Create object for distribution
            vec = pdf(pd, time);
        end
        
        function vec = ini_cdf(MO, time)
            %% int_0^t g(t) dt  cumulative distribution function for the initiation. Returns a vector
            % convert parmaters to entries for log-normal distribution
            mu = log((MO.par_ini(1)^2)/sqrt(MO.par_ini(2)+MO.par_ini(1)^2));
            sigma = sqrt(log(MO.par_ini(2)/(MO.par_ini(1)^2)+1));
            
            % Create object for distribution
            pd = makedist('Lognormal',mu, sigma);
            
            vec = cdf(pd, time);
        end
        
        
        function mat = kin_pdf(MO, time)
            %% NOT IMPLEMENTED f(x,t) probability density function for the kinetics where x is occupancy. Should return a matrix
            display('Distribution of NUP kinetics has not been defined')
            mat = [];
        end
        
        function vec = kin_mom_pdf(MO, time, imom)
            %% Moments of the distribution int_0^\infty f(x,t) x^imom. F(t) = kin_mom_pdf(time, 1)
            if imom == 1
                % create offset so that vec(1) at time = 0 is equal par(2)
                tmin = log(MO.par_kin(2)/(1 - MO.par_kin(2)))*MO.par_kin(1); % tmnin = 0 is par(2) = 0.4
                time = time + tmin;  % this is badly defined use logn_hill function instead
                vec = exp(time/MO.par_kin(1))./(1 + exp(time/MO.par_kin(1)));
            else
                display('Only first moment of kinetics has been defined');
            end
        end
    end
end