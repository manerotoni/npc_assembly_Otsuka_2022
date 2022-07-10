% Copyright MIT-LICENSE 2022 Antonio POLITI
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


%% class single_npc_gamma_hill
% A single npc model with deterministic hill kind of kinetics and gamma
% distribution initialization probability density function
% A gamma function can be justified by a multistep process and the time
% distrubution of reaching the final step. For example complete assembly of
% a NPC in order to bind the next pore.
% Zhou, Y., & Zhuang, X. (2007). Kinetic Analysis of Sequential Multistep Reactions Kinetic Analysis of Sequential Multistep Reactions. Journal of Physical Chemistry B, 111(48), 13600�13610. https://doi.org/10.1021/jp073708
% Hill function: t^n/(t^n + K^n)
% Author: Antonio Politi, EMBL Heidelberg, Max Planck Institute for multidisciplinary science, Goettingen.

classdef single_npc_gamma_hill < abs_single_npc
    % par_ini(1) : number of steps (if integer)
    % par_ini(2) : 1/rate ~ tau the average time for each step
    % par_kin(1) : the hill-coefficient. t^par_kin(1)./(time^par_kin(1) + par_kin(2)^par_kin(1));
    % par_kin(2) : the inflexion point
    properties(Access = protected)
        
        par_kin_lb = [0.5, 1];
        par_kin_hb = [10, 200];
        par_ini_lb = [1, 2];
        par_ini_hb = [200, 200];
        name = 'single_npc_gamma_hill';
    end
    
    properties (Constant)
        n_par_kin = 2;
        n_par_ini = 2;
    end
    methods
        function MO = single_npc_gamma_hill(varargin)
            display('single_npc_gamma_hill: Gamma function initiation, hill sigmoidal kinetics')
            MO = MO@abs_single_npc(varargin{:});
        end
        
        function vec = ini_pdf(MO, time)
            %% g(t) probability density function for the initiation. Returns a vector
            % convert parmaters to gamma distribution
            if MO.par_ini(1) == 0
                vec = time*0;
            else
                vec = gampdf(time, MO.par_ini(1), MO.par_ini(2));
            end
        end
        
        function vec = ini_cdf(MO, time)
            %% int_0^t g(t) dt  cumulative distribution function for the initiation. Returns a vector
            % convert parmaters to entries for gamma distribution
            vec = gamcdf(time, MO.par_ini(1), MO.par_ini(2));
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
            par_ini = [0; 5];
            par_kin = [2;10];
            time = [0:0.1:300]';
            amodel = single_npc_gamma_hill(par_ini, par_kin);
            amodel.plot_combined_dist(time, 1);
            % find a function that correspond to occupancy
            figure(2)
            clf
            hold
            plot(time, amodel.kin_ini_conv_mom(time, 1), 'r')
            par_kin = [5;61.5];
            plot(time, time.^par_kin(1)./(time.^par_kin(1) + par_kin(2).^par_kin(1)), 'b');
            
            
        end
    end
end