
classdef SCMState
% classdef SCMState

% Initialize global constants that define the locations of the steady
% states.  

% 10-Feb-2009 Only initialize constants that determine the steady states
%        don't initialize those that only alter stability (tau, gamma, v)
%
% Copyright (c) 2009-2012 Alistair Steyn-Ross: asr(at)waikato.ac.nz
%   School of Engineering, Waikato University, NZ 
% This code forms part of the SOinC (slow oscillations in cortex) project
%
% Feb 2009 modified July 2012
%
% Modified by M. A. Kramer for the Seizing Cortical Model (SCM)
% May 2016.
%
% Modified by E. D. Schlafly, June 2022


%% Methods

    methods
        function S = SCMState(obj, grid_size)
            
            if nargin < 1, return; end
            if nargin < 2, grid_size = []; end
            
            if ~isempty(obj) && ~isa(obj, 'SCMState')
                for ff = string(fieldnames(obj)'), S.(ff) = obj.(ff); end
            end
            
            if ~isempty(grid_size), S.resize(grid_size); end
            
        end

    end

%% Current state (IC)
    properties 
        %% Initial conditions
            % These are the same as in Waikato-Kramer except that {dVe, dVi,
            % Dii, Dee} no longer evolve dynamically (they are direct functions of
            % extracellular potassium). 

            
            K = 0;  % extracellular potassium concentration (cm^2)
            Qe = 0;  % Activity of excitatory population.
            Qi = 0;  % Activity of inhibitory population.
            Ve = -64 % Voltage  of excitatory population.
            Vi = -64 % Voltage of inhibitory population.
            Phi_ee = 0;  % e <--> e synaptic flux
            Phi_ei = 0;  % e <--> i synaptic flux
            Phi_ie = 0;  % i <--> e synaptic flux
            Phi_ii = 0;  % i <--> i synaptic flux
            phi2_ee = 0;  % Wave dynamics
            phi2_ei = 0;
            phi_ee = 0;
            phi_ei = 0;
            F_ee = 0;  % flux dynamics
            F_ei = 0;
            F_ie = 0;
            F_ii = 0;
            
            % Non-dynamic state variables (i.e., these values are overwritten
            % during simulation)
            Dee = 0.008;  % e <--> e gap-junction diffusive-coupling strength (electrodes)
            Dii = 0.8;  % i <--> i gap-junction diffusive-coupling strength in all space (electrodes)
            dVe = -1;  % Excitatory resting potential offset (mV)
            dVi = .1;  % Inhibitory resting potential offset (mV)
            

    end
    
    methods
        
        function S = resize(S, grid_size)
        % Converts scalar values to matrices
        
            for ff = string(properties(S)')
                if all(size(S.(ff)) == grid_size), continue; end
                S.(ff) = S.(ff)(1) .* ones(grid_size); 
            end
        end
        
    end
end




