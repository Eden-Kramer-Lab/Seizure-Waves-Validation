classdef SCM < handle
    
    methods  % constructor and other operators
        function scm = SCM(varargin)
            
            if nargin == 0, scm = SCM('IW'); return; end
            
            switch class(varargin{1})
                case 'SCM'
                    scm = varargin{1};
                    
                case 'char'
                    try
                        mdl = validatestring(varargin{1}, ...
                            {'steyn-ross', 'FS', 'SW', 'IW'});
                        varargin = varargin(2:end);
                        switch mdl
                            
                            case 'steyn-ross'
                                % This will perform similarly
                                scm.label = 'SR';
                                scm = scm.init();
                                scm.v = 140;
                                scm.Nee_sc = 7;
                                scm.Nei_sc = 7;
                                
                                dVe = 1.5;
                                Dii = .35;  % set this to whatever you like to change equilibrium dynamics
                                
                                scm.IZ(:) = dVe;
                                
                                % Adjust sigmoids so that voltage and Dii
                                % are fixed
                                scm.sigmoid_kD = [0 2 * Dii 0];        
                                scm.sigmoid_kdVi = [0 0 0]; 
                                scm.sigmoid_kdVe = [0 0 0];
                                
                                % Add a fixed source to push dynamics into
                                % steady high-firing state
                                scm.duration = 1;
                                center = [5 5]; 
                                radius = [1 1];
                                scm.source = scm.ellipse(center, radius, .5);
                                
                                % No IW
                                scm.map(:) = inf;
                                
                            case 'FS'
                                % Generates a sim with a fixed source. Uses
                                % [10, 10] padding and a duration of 30

                                scm.label = 'FS';
                                scm.save_output = true;
                                scm.visualization_rate = 0;
                                scm.padding = [10 10];
                                scm.duration = 30;
                                scm.dx = .1;
                                
                                % Mess with some other parameters
                                scm.sigmoid_kdVe(2) = .5;  % lower max K+ -> Ve offset (OK)
                                
                                % No IW
                                scm.map(:) = inf;

                                % Add a fixed source
                                center = round( [1.5 1.5] / scm.dx );
                                radius = round( [.25 .25] / scm.dx );
                                scm.source = scm.ellipse(center, radius, .5);

                                % Keep the electrodes in the center so that
                                % it's easy to rotate the source when
                                % generating sims for method comparison
                                scm.dimsNP = [4 4];
                                scm.centerNP = round( scm.grid_size ./ 2 );


                            case 'SW'
                                scm = SCM('FS');
                                scm.label = 'SW';
                                temp = scm.source;
                                scm.rotate(pi/2);
                                scm.source = cat(3, scm.source, temp);
                                scm.SwitchInterval = scm.duration / 2;

                                
                            case 'IW'
                                % This is the sim that shows up in the last
                                % figure
                                                                
                                scm.label = 'IW';
                                scm.basename = 'SCM/IW/IW';
                                scm.dx = 0.1;
                                scm.dimsNP = [4 4];
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);
                                scm.sim_num = 3;
                                scm.save_output = true;
                                                               
                                scm.visualization_rate = 10;
                                scm.depo_block = false;
                                scm.padding = [10 10];
                                scm.duration = 60;
                                                                
                                % Save and visualize some extra fields for testing
                                scm.out_vars = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};

                                % Design the IW
                                scm.expansion_rate = 0.1;

                                % Add a fixed source and IZ
                                scm.IZ = 1.2;
                                scm.set_layout('c');
                                scm.map = scm.generate_contagion;
                                
                                % Move the electrodes off the center
                                scm.centerNP = round( [3.5 3.5] / scm.dx );

                                
                            otherwise
                                error('Input ''%s'' not recognized.', mdl)
                        end
                    catch ME
                        if ~strcmpi(ME.identifier, 'MATLAB:unrecognizedStringChoice')
                            rethrow(ME)
                        end
                    end
                    for ii = 1:2:numel(varargin)
                        ff = validatestring(varargin{ii}, fieldnames(scm));
                        scm.(ff) = varargin{ii + 1};
                    end
            end
            
            scm = scm.init();
            
        end
        
        
        function rotate(scm, theta)
            % Rotate the locations of the source and stim_center about the
            % center of the simulated area by angle <theta>
            
            center = scm.grid_size/2;
            rot =@(xy) round( ...
                    (xy - center) ...  % translate to center
                    * [cos(theta), -sin(theta); sin(theta), cos(theta)] ...  % rotate
                ) + center;  % translate back

            scm.stim_center = rot(scm.stim_center);

            for ii = 1:size(scm.source, 3)
                mat = scm.source(:, :, ii);

                [Y, X] = find(true(size(mat)));
                xnew = rot([X, Y]);
                F = scatteredInterpolant(X, Y, mat(:), 'linear', 'nearest');
                rot_mat = reshape(F(xnew(:, 1), xnew(:, 2)), scm.grid_size);
                scm.source(:, :, ii) = rot_mat;
            end
        end

        
        function str = get.basename(scm)
            str = fullfile(scm.base_dir, scm.label, scm.label);
        end

        
        function P = init(P)
            
            if isempty(P.IC), P.IC = []; end  % set default IC
            P.IC = P.IC.resize(P.grid_size);
            
            % Check some values
            P.dimsNP = arrayfun(@min, P.dimsNP, P.grid_size);
            
            if isempty(P.centerNP)
                P.centerNP = P.grid_size / 2; 
            end
            P.recenterNP;
            
        end

        
        function set_layout(scm, idx)
            
            if nargin < 2, idx = 'a'; end
            iz = scm.IZ;
            [xx, yy] = find(iz > 0);
            C = mean([xx yy]);  % IZ center
            R = max(abs([xx yy] - C));  % IZ radius
            fsr = .25/scm.dx;  % fixed source radius
            
            r = [.3 .1] .* R;
            switch lower(idx)
                case 'a'
                    scm.stim_center = C + [1 -1] .* r;
                    scm.source = scm.ellipse(C + [1 1] .* r, fsr, 1);
                case 'b'
                    scm.stim_center = C + [1 -10] .* r;
                    scm.source = scm.ellipse(C + [1 -8] .* r, fsr, 1);
                case 'c'
                    scm.stim_center = C + [1 -2] .* r;
                    scm.source = scm.ellipse(C + [1 6] .* r, fsr, 1);
                case 'd'
                    scm.stim_center = C + [1 -10] .* r;
                    scm.source = scm.ellipse(C + [1 -2] .* r, fsr, 1);
                otherwise
                    error('Input ''%s'' not recognized', idx)
            end
        end
        
        
        function set_random_layout(scm)
            rng(scm.sim_num+500); % set the seed
            rand(30);  % blow out some numbers
            
            dVe_base = mode(scm.source, 'all');
            dVe_max = max(scm.source, [], 'all');

            scm.source = double(scm.ellipse()) * dVe_base;
            iz_inds = find(scm.ellipse());  % define the IZ
            N = numel(iz_inds);
            gs = scm.grid_size;
            [xx, yy] = meshgrid(1:gs(1), 1:gs(2));

            iw_loc = iz_inds(randi(N));
            fs_loc = iz_inds(randi(N));
            mea_loc = iz_inds(randi(N));

            scm.stim_center = [xx(iw_loc), yy(iw_loc)];  % get IW onset location
            scm.source(scm.ellipse([xx(fs_loc), yy(fs_loc)], 3)) = dVe_max;
            scm.centerNP = [xx(mea_loc) yy(mea_loc)];
            
            scm.init;
                                
        end
        
        
        function [fs, which_source] = fixed_source(scm, t)
            % Get the voltage offset for the fixed source at time t
            % If there are multiple fixed sources, rotate through them
            % every W seconds.
            W = scm.SwitchInterval;
            
            which_source = mod(floor(t / W), size(scm.source, 3)) + 1;
            if t >=0 && t < scm.duration
                fs = scm.source(:, :, which_source);
            else
                which_source = 0;
                fs = zeros(scm.grid_size);
            end
        end
        
        
        function plot_sigmoids(scm, x)
            % Displays the sigmoid functions for the potassium responses
            if nargin < 2, x = linspace(0, 2, 100); end
            ax = axes(figure, 'nextplot', 'add');
            for ff = {'sigmoid_kdVe', 'sigmoid_kdVi', 'sigmoid_kD'}
                y = scm.sigmoid(x, scm.(ff{:}));
                plot(ax, x, y, 'displayname', ff{1}(9:end))
            end
            legend(ax)
        end
        
        
        function clean(scm)
            % Removes previously saved simulations with the same basename
            % and sim number
            delete(sprintf('%s_%d_*.mat', scm.basename, scm.sim_num));
        end
        
        
        function map_ = generate_contagion(scm)
            % pre-generate the contagion style recruitment map
            rng(scm.sim_num+1);  % add 1 in case you use sim_num=0
            
            map_ = inf(scm.grid_size);
            if scm.expansion_rate <= 0, return; end
            dt_ = scm.dt * 10;  % lower precision is ok; worth increase in speed
            
            % Bernoulli => geometric distribution of waiting times
            k = scm.dx / dt_ / scm.expansion_rate;
            p_wavefront = (1/k)/2;  % divide by 2 to account for 2 dimensions of spread
            
            recruited = false(scm.grid_size);  % set state negative to indicate "not recruited"
            recruited(scm.stim_center(1), scm.stim_center(2)) = true;  % stim center is recruited at t = 0;
            map_(recruited) = 0;
            tt = dt_;  % advance tt
            
            % repeat until end of simulation or all nodes are recruited
            while tt < scm.duration && any(isinf(map_), 'all')
                                
                % find non-recruited points with recruited neighbors
                boundary = conv2(recruited, [0 1 0; 1 -4 1; 0 1 0], 'same') > 0;  
                
                % advance wavefront randomly
                dice = rand(size(recruited));  
                wavefront = (dice < p_wavefront) & boundary;
                
                recruited(wavefront) = true;
                
                map_(wavefront) = tt;
                
                tt = tt + dt_;
            end

        end
        
        
        function map_ = generate_map(scm)
            % Generates a linear radial distance from source map with some
            % 2D gaussian noise
            map_ = inf(scm.grid_size);
            if scm.expansion_rate <= 0, return; end
            [xx, yy] = find(ones(scm.grid_size));
            xx = xx - scm.stim_center(1);
            yy = yy - scm.stim_center(2);

            T = reshape( ...
                scm.expansion_rate / scm.dx * sqrt(xx.^2 + yy.^2), ...
                scm.grid_size);

            noise_ = randn(size(T)) * 2;
            f_smooth = round(.5/scm.dx);
            win = gausswin(f_smooth) * gausswin(f_smooth)';
            map_ = T + conv2(noise_, win./sum(win, 'all'), 'same');
            map_(map_ < 0) = 0;

        end
                       
        
        function recenterNP(scm)
            % Makes sure that NP does not go outside of simulation bounds
            ctr = scm.centerNP;
            xx = 1:scm.grid_size(1);
            yy = 1:scm.grid_size(2);
            [~, ext_x] = sort(abs(xx - ctr(1))); ext_x = sort(xx(ext_x(1:scm.dimsNP(1))));
            [~, ext_y] = sort(abs(yy - ctr(2))); ext_y = sort(yy(ext_y(1:scm.dimsNP(2))));
            scm.centerNP = mean([ext_x', ext_y']);
        end
        
        
        function [inds, ext_x, ext_y] = NPinds(scm)
            % Computes the indices extracted as the MEA. 
            
            mea_mask = false(scm.grid_size);
            ctr = scm.centerNP;
            
            xx = 1:scm.grid_size(1);
            yy = 1:scm.grid_size(2);
            [~, ext_x] = sort(abs(xx - ctr(1))); ext_x = sort(xx(ext_x(1:scm.dimsNP(1))));
            [~, ext_y] = sort(abs(yy - ctr(2))); ext_y = sort(yy(ext_y(1:scm.dimsNP(2))));
            mea_mask(ext_x, ext_y) = true;
            inds = find(mea_mask);
    
            
        end
        
        
        function inds = ECinds(scm)
            % Never used this... might need help
            x_offset = ( (1:scm.dimsEC(1)) - floor(scm.dimsEC(1)/2) ) * scm.scaleEC;
            y_offset = ( (1:scm.dimsEC(2)) - floor(scm.dimsEC(2)/2) ) * scm.scaleEC;
            [xx, yy] = ndgrid(x_offset, y_offset);

            inds = sub2ind(scm.grid_size, ...
                scm.centerNP(1) + xx(:), ...
                scm.centerNP(2) + yy(:));
        end
        
        
        function h = show_layout(scm)
            
            
            h = figure('name', 'layout');
            ax = axes(h);
            
            % Set colors
            colors = cool(4);
            cc_ = containers.Map( ...
                ["full", "IZ", "source", "mea"], ...
                arrayfun(@(ii) colors(ii, :), 1:4, 'uni', 0));
            gs = scm.grid_size;
            
            % Show simulation bounds
            patch(ax, [0; gs(1); gs(1); 0], [0; 0; gs(2); gs(2)], ...
                cc_("full"), 'displayname', 'Sim');
            text(ax, gs(1), gs(2), 'Simulated region', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
            hold(ax, 'on')
            
            % Show irritative zone (IZ)
            [x, y] = find(scm.IZ);
            k = boundary(x, y);
            patch(ax, x(k), y(k), cc_("IZ"), 'displayname', 'IZ');
            text(ax, min(x(k)), mean(y(k)), 'IZ', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            
            % Show source(s)
            for d3 = 1:size(scm.source, 3)
                [x, y] = find(scm.source(:, :, d3));
                k = boundary(x, y);
                patch(ax, x(k), y(k), cc_("source"), ...
                    'displayname', sprintf('FS %d', d3));
                text(ax, mean(x(k)), mean(y(k)), sprintf('FS %d', d3), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            
            % Show MEA
            [x, y] = ind2sub(gs, scm.NPinds);
            k = boundary(x, y);
            patch(ax, x(k), y(k), cc_("mea"), 'displayname', 'MEA')
            text(ax, mean(x(k)), mean(y(k)), 'MEA', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            
            % Show IW source
            plot(ax, scm.stim_center(1), scm.stim_center(2), 'r*', ...
                'displayname', 'IW source')
            text(ax, scm.stim_center(1), scm.stim_center(2), 'IW source', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            
            % Cleanup
            hold(ax, 'off');
            set(findobj(ax, 'type', 'Patch'), 'FaceAlpha', .7, 'edgecolor', 'none')
            view(ax, 90, 90);  % match what imagesc shows during visualization
            axis(ax, 'equal')
        end        
        
        
        function angle = theta_FS(scm)
            % Computes the angle between the fixed source and the
            % microelectrode array
            angle = nan(size(scm.source, 3), 1);
            for ii = 1:numel(angle)
                [xx, yy] = find(scm.source(:, :, ii) == max(scm.source(:, :, ii), [], 'all'));
                dy = scm.centerNP(2) - mean(yy);
                dx_ = scm.centerNP(1) - mean(xx);
                angle(ii) = atan2(dy, dx_);
            end
        end
        
        
        function angle = theta_IW(scm)
            % Computes the angle between the IW onset point and the
            % microelectrode array
            dy = scm.centerNP(2) - scm.stim_center(2);
            dx_ = scm.centerNP(1) - scm.stim_center(1);
            angle = atan2(dy, dx_);
        end
        
        
        function Run(self)
            self.CreateDirectory;
            self.RunSimulation;
            self.ConvertToMea;
        end
        

        RunSimulation(self)
        ConvertToMea(self)
        h = Preview(self, show_layout)
    end
    
    methods (Hidden = true)
        function CreateDirectory(scm)

            if ~scm.save_output, return, end

            [base_path, ~, ~] = fileparts(scm.basename);
            if ~exist(base_path, 'dir'), mkdir(base_path), end
            fname = sprintf('%s_%d_info', scm.basename, scm.sim_num);
            save(fname, 'scm');
            fprint('Sim params saved to %s.\n', fname)
        end

    end
    
    properties 
        Qe_movie
        Ve_movie
        K_movie
    end
    
    properties  % updates
        
        % Rotate through sources every SwitchInterval seconds (only
        % applies when there are multiple sources)
        SwitchInterval = 2
        
        % Epileptogenic zone and fixed source. A baseline increase in excitability over a
        % subregion of the sim. Default shapes set in the getters
        IZ = 1.2
        
        % A predetermined IW map. Generate this if the IW is independent of
        % the rest of the sim. It should be a grid of recruitment times. I
        % tried using this with a t:=Cr function (with 2D gaussian noise;
        % see <SCM.generate_map>) and the recruitment pattern looked pretty
        % good, but the smoothness of it made it so no TW formed after the
        % IW until the FS was inside the recruited region
        map  
        
        % Parameters for the inhibitory collapse portion of the IW. The
        % model assumes an increase in inhibition (depression of dVe)
        % followed by inhibitory collapse, modeled here as a drop in the
        % max inhibitory firing rate. 
        %     Used in: SCM.Qi_max_fun
        %     Elements: 
        %       drop: the amount by which the firing rate decreases
        %       offset: time (s) of peak drop relative to duration of dVe depression
        %       width: sd of gaussian used to model the drop 
        % 
        Qi_collapse = [20 1 .5]  % [drop, offset, width]
        
        
        % SIGMOIDS
        % Change dV* and Dii from dynamic functions of [K+]o to direct
        % sigmoid response functions. 
        sigmoid_kdVe = [0.5 .7 10]  % [mid max slope] (width is ~max/slope*4 if you prefer to think of it that way)
        sigmoid_kdVi = [.5 .3 10]
        sigmoid_kD = [1 .3 -4]
        
        
        % For simplicity, remove (for now) the depolarization block
        % function that was added in Martinet. With this parameter set, the
        % voltages never get high enough for it to take effect. Maybe look
        % into this more later.
        depo_block = 0  % apply depolarization block
        
    end
    
    properties  % original-ish
        % meta
        label (1,:) char = 'SCM'
        base_dir (1, :) char = fullfile(pwd, 'SCM')
        basename
        mea_path
        sim_num  
        save_output = true  % Save output
        visualization_rate = 0  % Show this many frames per second
        t_step = 1  % Simulate <t_step> second intervals
        t0_start  % Allows continue from previous sim (IC will use last from file number t0_start-1)
        duration (1,1) double {mustBeNonnegative} = 60
        padding (1,2) = [10 10]  % Padding before and after seizure
        source_drive (1, 1) double = 2.5
        subsample (1,1) double {mustBePositive} = Inf  % Allow downsampling when creating mea
        return_fields (1,:) = {'Qe', 'Ve'}  % Qe required to make mea
        out_vars (1,:) = {'Qe', 'Ve'}  % Define which variables you would like to visualize (can be any in IC)
        source % Define fixed, alternating sources of excitation
        t0  % used to keep track of progress in sims 
        
        
        % model
        stim_center (1,2) = [20 20]
        grid_size (1,2) {mustBePositive} = [50 50] % size of grid to simulate (must be even)
        dt = 2e-4
        dx = .4  % (cm) 
        expansion_rate (1,1) double {mustBeNonnegative} = .1  % in cm^2/s; set to 0 for fixed source
        IC SCMState 
        
        
        % steady states
        tau_e = 0.04  % excit neuron time-constant (/s) [original = 0.04 s, Martinet = 0.02 s]
        tau_i = 0.04  % inhib neuron time-constant (/s) [original = 0.04 s, Martinet = 0.02 s]

        
        % voltage limits
        Ve_rev = 0  % reversal potential (mV)
        Vi_rev = -70
        Ve_rest = -64  % resting potential (mV)
        Vi_rest = -64

        
        % gain per synapse at resting voltage (millivolt.sec)
        rho_e = 1.00e-3
        rho_i = -1.05e-3

        
        % E/IPSP rate constants
        gamma_e = 170  % EPSP decay rate (/s)
        gamma_i = 50  % IPSP decay rate (/s)

        
        % gap-junction diffusive-coupling strength
        D = 0.8  % i <--> i steady state (cm^2) [0, 1]  [Martinet = 1]

        
        % sigmoid characteristics
        Qe_max = 30  % sigmoid maximum (s^-1)
        Qi_max = 60
        theta_e = -58.5     % sigmoid threshold (mV)
        theta_i = -58.5
        sigma_e = 3.0    % sigmoid 'width' (mV)
        sigma_i = 5.0

        
        % connectivities: j-->k convention (dimensionless)            
        Nee_a = 2000  % cortico-cortical
        Nei_a = 2000
        Nee_b = 800
        Nei_b = 800
        Nie_b = 600
        Nii_b = 600
        Nee_sc = 50
        Nei_sc = 50
        

        % axonal conduction velocity (cm/s), 
        v = 280  % [original = 140 cm/s]

        
        % inverse-length scale for connectivity (/cm)
        Lambda = 4.0
        
        
        % potassium
        tau_K = 200    %time-constant (/s).
        k_decay = 0.1  %decay rate (/s).
        kD = 1       %diffusion coefficient (cm^2/s).
        KtoVe = 10     %impact on excitatory population resting voltage.
        KtoVi = 10     %impact on inhibitory population resting voltage.
        KtoD = -50    %impact on inhibitory gap junction strength.
        kR = 0.15   %scale reaction term. 

        
        % electrodes
        centerNP  % defaults to grid center
        dimsNP = [10 10]

        
        % Macroscale
        centerEC  % defaults to grid center
        scaleEC = 4
        dimsEC = [3 3]
        
        
        % noise
        noise_sf (1,1) double {mustBeNonnegative} = 4  % [Martinet = 2]
        noise_sc (1,1) double {mustBeNonnegative} = 0.2;
        
    end

    
    properties (Dependent = true)
        
        % default subcortical fluxes
        % %% ORIGINAL FORMULATION %%
        % [Nee_sc,Nei_sc]= deal(50, 50)  % subcortical  
        % phi_ee_sc = Nee_sc * Qe_max  % original [1500]
        % phi_ei_sc = Nei_sc * Qe_max  % original [1500]

        % %% EDS %%
        phi_ee_sc
        phi_ei_sc
        
        % d/dV derivatives of psi_ij weighting functions
%         d_psi_ee 
%         d_psi_ei 
%         d_psi_ie 
%         d_psi_ii 
%         
%         % Nee and Nie totals for cortico-cortical plus intracortical
%         Nee_ab
%         Nei_ab

    end
    
    
    methods  % Emily parameter functions
        
        function FS = get.source(scm)
            % Default the fixed source to a .25x.25 cm ellipse centered at
            % [2.2, 2.2] cm with the value given in scm.FS (this is in
            % addition to IZ)
            if isempty(scm.source), scm.source = 0; end
            if numel(scm.source) == 1
                center = round( [2.2 2.2] / scm.dx );
                radius = round( [.25 .25] / scm.dx );
                scm.source = double(scm.ellipse(center, radius)) * scm.source;
            end
            FS = scm.source;
        end
        
        function IZ = get.IZ(scm)
            % Default the epileptogenic zone to an ellipse with the value
            % given in scm.IZ
            if numel(scm.IZ) == 1
                scm.IZ = double(scm.ellipse()) * scm.IZ;
            end
            IZ = scm.IZ;
        end
        
        
        function [Qi_max, dVe] = iw_model(scm, time)
            % Model of inhibitory collapse. Max inhibitory firing rate
            % collapses following a gaussian (this was chosen
            % only because it gives a smooth drop and recovery; no proposed
            % relation to mechanism). 
            
            % These are baselines
            Qi_max = scm.Qi_max;
            dVe = 0;
            
            if nargin < 2, return; end  % return now if no time is given
            
            % <state> represents time since IW onset on each node
            state = time - scm.map; 
                    
            % ... and then collapses following a gaussian curve
            % with parameters defined in <scm.Qi_collapse>
            Qi_max = scm.Qi_max ...  
                - scm.Qi_collapse(1) * scm.gaussian( ...
                    state, scm.Qi_collapse(2), scm.Qi_collapse(3));
            
        end
        

    end
    
    methods  % Computed fixed parameters
        function phi_ee_sc = get.phi_ee_sc(self)
            phi_ee_sc = self.Nee_sc * self.Qe_max;
        end
        
        function phi_ei_sc = get.phi_ei_sc(self)
            phi_ei_sc = self.Nei_sc * self.Qe_max;
        end

    end
    
    methods  % functions for convenience/setup
        
        function map = get.map(scm)
            % Default to contagion map
            if isempty(scm.map), scm.map = scm.generate_contagion; end
            map = scm.map;
        end
        
        
        function str = get.mea_path(scm)
            str = fullfile(scm.base_dir, ...
                sprintf('%s_Seizure%d_MEA_%d_%d.mat', ...
                    scm.label, scm.sim_num, scm.padding));
        end
        
        
        function map = ellipse(scm, C, R, value)
            % map = scm.ellipse(center=grid_size/2, radius=grid_size/2-4)
            % Returns a binary map matching scm.grid_size with an ellipse
            % centered at C (1x2) with radius R (1x2)
            
            if nargin < 2 || isempty(C),  C = scm.grid_size ./ 2; end
            if nargin < 3 || isempty(R); R = floor(scm.grid_size / 2 - 4); end
            if nargin < 4 || isempty(value); value = 1; end
            
            if numel(R) == 1, R = [R R]; end
            [ix, iy] = ind2sub(scm.grid_size, 1:prod(scm.grid_size));

            map = zeros(scm.grid_size);
            ellipse = sum( ([ix' iy'] - C).^2 ./ R.^2, 2 ) < 1;
            map(ellipse) = value;
        end
        
        
        function t0s = get.t0_start(p)
            if isempty(p.t0_start) || p.t0_start < -p.padding(1)
                p.t0_start = -p.padding(1);
            end
            t0s = p.t0_start;
        end        
        
        
        function num = get.sim_num(P)
            s = 1;
            while isempty(P.sim_num)
                if isempty(dir([P.basename '_' num2str(s) '*.mat']))
                    P.sim_num = s; 
                end
                s = s+1;
            end
            num = P.sim_num;
        end
        
        
        function grid_size = get.grid_size(p)
            if any(mod(p.grid_size, 2)), p.grid_size = p.grid_size + mod(p.grid_size, 2); end
            grid_size = p.grid_size;
        end

        
        function center = get.centerEC(P)
            if isempty(P.centerEC)
                P.centerEC = round(P.grid_size / 2);
            end
            center = P.centerEC;
        end
        
    end
    

    
    %% Static methods
    methods (Static)
        function y = sigmoid(x, p)
            y = p(2) ./ (1 + exp(-p(3) * (x - p(1))));
        end

        
        function y = gaussian(x, mu, sigma)
            y = exp(-(x - mu).^2 ./ (2 .* sigma.^2));
        end
        
        
        function center = get_center(center, dims, grid_size)
            % (Static) center = get_center(center, dims, grid_size) 
            if isempty(center)
                center = round(grid_size / 2); 
            end
            center = min(max(center, ceil(dims / 2)), grid_size - floor(dims / 2));
        end
    end
    

end

