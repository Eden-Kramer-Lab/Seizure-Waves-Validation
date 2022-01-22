classdef WaveProp < handle & matlab.mixin.Copyable
    
    
	properties 
		
        Patient string
        Seizure string
        Position
        UseLargestCluster = true
        
        GridSize = [10 10]  % Size of the MEA (# of electrodes along each dimension; always [10 10] unless it's a sim)
		Length = 4  % length of the MEA (always 4 mm, unless it's a sim)
        MMpElectrode = .4  % mm per electrode (again, always .4 unless it's a sim)
        
        TOA  % this is the important one to save (especially for delays method since it takes a long time to compute
        TOAcube  % converts TOA into a TxMxM array
        t0

    end
    

	properties (Transient = true)

        ResultOfFit
        NClust  % I don't even know what this does anymore (5/9/21)
		ClustSize
        
        sig = 0.05
        RotateBy = 0
        MinFinite = 30
        Quadratic = false
		Original = false
		Magnitude
        AltMagnitude
		Direction
        DirAlt
		
		scale_quiver = 1
		complexZ
		Inds = []
		Early = [-Inf 25]
		Late = [30 Inf]
		MinDetections = 100
        
        Vx 
		Vy
		p
		Beta
	end
	
	
	properties (Dependent = true)

        Name string
		Z
        Time
		mask
		logp
		NormedMagnitude
% 		Phi_mean_early
% 		Phi_mean_late
		Phi_std_early
		Phi_std_late
		First_detection
		N_detections_early
		N_detections_late
        FName
        Data
        locs
        
        
    end
    
	
	methods  
        %% Getters/setters
        function set.Name(obj, value)
            value = string(value);
            ps = strsplit(value, {'_Seizure', ' ', '_'});  % Assume PP_SeizureXX assignment
            obj.Patient = ps(1);
            obj.Seizure = ps(2);
        end
        
        function name = get.Name(obj)
            name = sprintf("%s_Seizure%s", obj.Patient, obj.Seizure);
        end
        
        function [ci_dir, ci_sp] = direction_ci(M)
            % Estimate a confidence interval around the estimated direction
            % by bootstrapping (as in Martinet, 2017)
            
            res = M.ResultOfFit;
            N = size(res, 1);
            
            [ci_dir, ci_sp] = deal(nan(N, 2));
            
            for ii = 1:N
                if isnan(res.stats(ii).covb(1)), continue; end
                samples = mvnrnd(res.beta(ii, 2:3), res.stats(ii).covb(2:3,2:3), 1000);       % generate samples from model.
                boot_dir = atan2(samples(:,2), samples(:,1));                   % bootstrap direction
                boot_sp = 1./sqrt(samples(:,1).^2 + samples(:,2).^2);           % bootstrap velocity.

                mn = circ_mean(boot_dir);
                ci_dir(ii, :) = ...
                    quantile( ...
                        angle(exp(1j * (boot_dir - mn))), ...
                    [0.025, 0.975]) ...
                    + mn;                    % CI for direction
                ci_sp(ii, :) = quantile(boot_sp, [0.025, 0.975]);  
            end
            ci_dir = angle(exp(1j * (ci_dir - M.RotateBy)));
        end
        
        function out = get.ResultOfFit(obj)
            % If fits haven't been run, run them now (fit a plane to each
            % time)
            
            if isempty(obj.ResultOfFit)
                inds_ = obj.Inds; 
                obj.Inds = [];
                
                % Get the observed TOA
                if obj.UseLargestCluster
                    toa = obj.use_largest_cluster(obj.TOAcube);
                    toa = toa(:, obj.locs);
                else
                    toa = obj.TOA;
                end
                Nt = size(obj.TOA, 1);
                
                % Initialize matrices to store the results and fit a plane
                % at each time point
                V = nan(Nt, 2); p_ = nan(Nt, 1); beta = nan(Nt, 3);
                stats_temp = cell(Nt, 1);
                for ii = 1:Nt
                    [V(ii, :), p_(ii), beta(ii, :), stats_temp{ii}] = ...
                        obj.fit_plane(toa(ii, :), obj.Position);
                end
                
                % Clean up stats object
                mask_ = cellfun(@isstruct, stats_temp);
                tpl = stats_temp{find(mask_, 1)};
                for ff = string(fields(tpl))', tpl.(ff)(:) = nan; end
                stats_temp(~mask_) = {tpl};
                stats = cat(1, stats_temp{:});
                
                % Save the result to a table
                obj.ResultOfFit = table(V, p_, beta, stats);
                
                % Reset inds to whatever they were
                obj.Inds = inds_;
            end
            out = obj.ResultOfFit;
            
        end
        
        function toa = get.TOA(obj)
            if ndims(obj.TOA) == 3
                obj.TOA = obj.TOA(:, obj.locs);
            end
            toa = reshape(obj.TOA, numel(obj.t0), size(obj.Position, 1));
            toa = toa(obj.Inds, :);
        end
        
        function res = get.TOAcube(obj)
            % Convert TOA to 3D (and return .Inds)
            
            m = numel(obj.Time); 
            toa_size = [m, obj.GridSize];  % result should be TIMExXxY
            
            res = nan(toa_size);  % initialize matrix to hold the result

            % Get the indices where the TOA data belongs and generate the
            % resulting TOA
            [tt, xx] = ndgrid(1:m, obj.Position(:, 1));
            [~, yy] = ndgrid(1:m, obj.Position(:, 2));

            locs_ = sub2ind(toa_size, tt(:), xx(:), yy(:));
            res(locs_) = obj.TOA;

        end
        
        function data = get.Data(self)
            % Wrapper for self.TOAcube
            % renamed Data to TOA to be more clear about what's in here
            
            data = self.TOAcube;
        end
        
        function set.Data(self, value)
            % Data is a wrapper for TOA (more clear naming)
            self.TOA = value;
        end
        
        function fname = get.FName(self)
            % Returns the name of the file where the TOA are stored
            fname = sprintf("WaveFits/%s_Seizure%s/", self.Patient, self.Seizure);
            warning('You probably want to call .Name here rather than .FName');
%             fname = sprintf('%s_Seizure%s_fits.mat', self.Patient, self.Seizure);
        end
        
        function pat = get.Patient(self)
            if isempty(self.Patient)
                info = strsplit(self.Name, {'_', 'Seizure'});
                self.Patient = info{1};
                self.Seizure = info{2};
            end
            pat = self.Patient;
        end
        
        function sz = get.Seizure(self)
            if isempty(self.Seizure)
                info = strsplit(self.Name, {'_', 'Seizure'});
                self.Patient = info{1};
                self.Seizure = info{2};
            end
            sz = self.Seizure;
        end
        
		function inds = get.Inds(s)
			inds = s.Inds;
			if isempty(inds)
				inds = (1:numel(s.t0))';
			end
        end
        
        function set.Time(self, value)
            self.t0 = value;
        end
        
		function beta = get.Beta(obj)
			beta = obj.ResultOfFit.beta(obj.Inds, :);
        end
        
        function vx = get.Vx(obj)
            vx = obj.ResultOfFit.V(:, 1);
        end
        
        function vy = get.Vy(obj)
            vy = obj.ResultOfFit.V(:, 2);
        end
        
		function mgn = get.NormedMagnitude(obj)
			mgn = filloutliers(obj.Magnitude, nan, 'median', 'thresholdfactor', 10);
			mgn = normalize(mgn);
			mgn = mgn(obj.Inds);
        end
        
		function logp = get.logp(s)
            % returns XX from expression p=5e-XX
			logp = -log10(s.p/5);
			logp = logp(s.Inds);
        end
        
		function Z = get.complexZ(s)
			Z = complex(s.Vx, s.Vy);
			Z = Z(s.Inds);
        end
        
        function minfinite = get.MinFinite(s)
            if numel(s.MinFinite) > 1
                s.MinFinite = max(s.MinFinite);
            end
            minfinite = s.MinFinite;
        end
        
		function mask = get.mask(s)
            M = sqrt(s.Vx.^2 + s.Vy.^2);
            outs = isoutlier(M, 'thresholdfactor', 100); 
            if mean(outs) < .2   % remove outliers as long as they don't make up more than 20% of the data
                M(outs) = nan;
            end
            
            % require full window (set to .95 bc a few of the windows are short a ms)
            if ismember('T', properties(s))  
                M(s.T < .95*2*s.HalfWin) = nan;
            end
			M = M(s.Inds);
            num_finite = sum(isfinite(s.Data), [2 3]);
            M(num_finite < s.MinFinite(1)) = nan;
			mask = s.p > s.sig | isnan(M);
            
            
            % mask times when VNS is active in P2
            if strcmpi(s.Name, 'P2_seizure1')
                vns_times = s.Time >= 13 & s.Time <= 29;
                mask = mask | vns_times;
            end
            if strcmpi(s.Name, 'P2_seizure2')
                vns_times = s.Time >= 46 & s.Time <= 62;
                mask = mask | vns_times;
            end
        end
        
        function t = time(s)
            t = s.Time;
        end
        
		function time = get.Time(s)
			time = s.t0;
			time = time(s.Inds);
        end
                
		function RB = get.RotateBy(s)
			RB = s.RotateBy .* ones(size(s.Vx));
			RB = RB(s.Inds);
        end
        
		function D = get.Direction(M)
            
            if isempty(M.Direction)  % this should usually be empty so it pulls from (Vx, Vy) but sometimes I want WaveProp functions on a specific set of directions
                D = atan2(M.Vy, M.Vx);

                inds_ = M.Inds;
                M.Inds = [];
                D(M.mask) = nan;

                % Remove points where directions show large shifts within 50
                % ms. (this shouldn't do anything to M and D10 methods, but some
                % experimental methods [e.g. M10] benefit from this)
                D = D - M.RotateBy;
%                 dir_sm = movmean(exp(1j*D), .05, 'samplepoints', M.t0);  
%                 dir_sm(abs(dir_sm) < cos(pi/8)) = nan; % this ensures consecutive angles differ by less than 45°
%                 d2 = movmean(dir_sm, .04, 'omitnan', 'samplepoints', M.Time);  % waves traveling at 100mm/s should be on the MEA ~40 ms
%                 d2(abs(d2) < .9) = nan;  % require strong similarity in directions
%                 D = angle(d2) - M.RotateBy;
                
                
                D = angle(exp(1j * D));  % Keep in range [-pi pi]
                M.Inds = inds_; 
                D = D(inds_);
            else
                D = M.Direction;
            end
        end
        
		function M = get.Magnitude(self)
            % in mm/s
			M = sqrt(self.Vx.^2 + self.Vy.^2) * self.MMpElectrode;  % Fits return electrodes/second. Convert to mm/s
			M = M(self.Inds);
            M(self.mask) = nan;
        end
                
		function psig = get.p(s)
			psig = s.ResultOfFit.p_(s.Inds);
		end
		
        function ll = get.locs(M)
            % convert x-/y-position to indices
            ll = sub2ind(M.GridSize, M.Position(:, 1), M.Position(:, 2));
        end
        
        
        %% Helpful others
        tbl = stable_intervals(s, varargin)
        
        function obj = parse_inputs(obj, varargin)
			for ii = 1:2:numel(varargin)
                ff = validatestring(varargin{ii}, properties(obj));
				obj.(ff) = varargin{ii+1};
			end
        end
        
        function t_inds = time2inds(M, tt)
            [~, t_inds] = min(abs(M.Time(:) - tt(:)'));
        end        
        
        function [DI, N, win] = directionality_index(M, win)
            % dir_index = directionality_index(M, win=[])
            % Wrapper to get directionality index from F.smoothed_direction
            if nargin < 2, win = inf; end
            
            if isinf(win)
                dir = rmmissing(M.Direction);
                DI = circ_r(dir);
                N = numel(dir);
            else
                dir = M.Direction;
                tt = M.Time;
                mask_ = isfinite(dir);
                [DI, N] = circ_mov_stat('r', dir(mask_), tt, win, ...
                    'samplepoints', tt(mask_));
            end
            
        end

        function [mn_est, ci_est] = mov_mean_dir(M, win, times)
            if nargin < 3, times = M.time; end
            if nargin < 2 || isempty(win), win = 10; end
            
            WNG = warning;
            warning('off', 'circ_confmean:requirementsNotMet');
            
            dir = M.Direction;
            mask_ = isfinite(dir);
            
            mn_est = circ_mov_stat('mean', dir(mask_), times, win, ...
                'samplepoints', M.Time(mask_));
            ci_est = circ_mov_stat('confmean', dir(mask_), times, win, ...
                'samplepoints', M.Time(mask_));
            
            warning(WNG);
        end
        
        function M = refit_data(M)
			M.Inds = [];  % Remove time resampling
            M.ResultOfFit = [];
            M.ResultOfFit;
            
        end
                
		function [mn, ci] = mean(s, data)
			if nargin < 2, data = s.Direction; end
            if all(isnan(data)), mn = nan; ci = nan; return; end
			[mn, ul] = circ_mean(data, [], [], 'omitnan');
			ci = ul - mn;
        end
        
        function var_ = var(s, data)
            if nargin < 2, data = s.Direction; end
            if all(isnan(data)), var_ = nan; return; end
            [~, var_] = circ_var(data(~isnan(data)));
        end
        
        function sd = std(s, data)
            if nargin < 2, data = s.Direction; end
            if all(isnan(data)), sd = nan; return; end
            sd = circ_std(data(~isnan(data)));
        end
        
		function d = diff(s, other)
            i0 = other.Inds;
			other = other.resample_t0(s.Time);
			d = angle(exp(1j*(s.Direction - other.Direction)));
            other.Inds = i0;  
        end
        
		function [obj, iq] = resample_t0(obj, t_new)
            
			iq = interp1(obj.t0, 1:numel(obj.t0), t_new, 'nearest', 'extrap');
			obj.Inds = iq(:);
        end	
        
        		
		%% Plotting
        function [ax, mn, ll, ul] = plot_dirs_with_ci(M, ax, win, units, varargin)
            % Plots directions with 95%CI using a <win> second sliding
            % window.
            % [ax, mn, ll, ul] = plot_dirs_with_ci(M, ax=gca, win=10, units="deg" ::plot args::)
            
            if nargin < 2 || isempty(ax), ax = gca; end
            if nargin < 3 || isempty(win), win = 10; end
            if nargin < 4 || isempty(units), units = "deg"; end
            
            if ~isa(ax, 'matlab.graphics.axis.Axes')  % bump everything forward
                varargin = [units varargin];
                units = win;
                win = ax;
                ax = gca;
            end
            
            if ~isnumeric(win)  % if win isn't numeric, assume it's omitted and bump again
                varargin = [units varargin];
                units = win;
                win = 10;
            end
            
            try % see if units actually units or part of varargin
                units = validatestring(units, ["rad" "deg"]);
            catch ME
                if isempty(varargin), rethrow(ME);
                else
                    varargin = [units varargin];
                    units = "deg";
                end
            end
            
            sts = ax.NextPlot;
            if strcmpi(sts, 'replace'), cla(ax); end  % the patches won't clear the axis first
            
            if units == "rad", convert_units = @(x) x; cycle = 2*pi; 
            else, convert_units = @rad2deg; cycle = 360;
            end
            
            switch class(M)
                case 'Delays'
                    mn = M.Direction;
                    ci = M.direction_ci;
                    tt = M.Time;
                    
                    % Unwrap so the plot looks nicer and get the right
                    % units
                    ci = convert_units(angle(exp(1j*(ci - mn))));
                    mn = convert_units(unwrap(mn));
                    ci = ci + mn;
                    ul = ci(:, 1); ll = ci(:, 2); 
                    
                otherwise
                    tt = ( M.Time(1):.1:M.Time(end) )';  % estimate mean and ci at 100 ms intervals
                    [mn, ci] = M.mov_mean_dir(win, tt);
                    mn(isnan(ci) & mn == 0) = nan;
                    
                    % Fill missing ci values, unwrap the directions, and convert to
                    % the requested units.
                    ci = fillmissing(ci, 'constant', pi, 'endvalues', 'none');  % we have no info where ci are nan, so the CI is pi
                    ci = convert_units(ci);
                    mn = convert_units(unwrap(mn));
                    ll = mn - ci;
                    ul = mn + ci;
            end
            
            
            % Plot the mean and 95%CI
            patch_args = varargin;
            patch_args(1:2:end) = strrep(lower(patch_args(1:2:end)), 'color', 'facecolor');
            
            for offset = [-cycle 0 cycle]
                mask_ = isfinite(ll + ul);
                yy = [ll(mask_); flipud(ul(mask_))] + offset;
                xx = [tt(mask_); flipud(tt(mask_))];
                pp = patch(ax, xx, yy, 1, 'displayname', '95%CI', ...
                    'facealpha', .5, 'linestyle', 'none', 'facecolor', [0 0 0], ...
                    patch_args{:}); %#ok<NASGU>
            end
            
            ax.NextPlot = 'add';
            plot(ax, tt, mn + [-cycle 0 cycle], 'k', ...
                'linewidth', 2, 'displayname', 'Mean', varargin{:});
            
            
            % Prettify
            if units == "deg"
                ylim(ax, [-180 180]);
                lbl = cellfun(@(ll) sprintf('%s°', ll), ax.YTickLabel, 'uni', 0);
                ax.YTickLabel = strrep(lbl, '°°', '°');
            else
                ylim(ax, [-pi pi]);
                intv = pi/3;
                tks = -3*pi:intv:3*pi;
                yticks(tks);
                polar_label(ax, 'y');
            end
            title(strrep(M.Name, '_', ' '))
            xlabel('Time [s]')
            ylabel('Direction')
            ax.NextPlot = sts;
            
            
        end
        
		function ax = plot(obj, varargin)
			% ax = plot(obj, type='2D', t0, ax); % if object has only one time point, t0=obj.Time
			directive = '';
			chars = cellfun(@ischar, varargin);
			if any(chars), directive = varargin{chars}; varargin(chars) = []; end
			switch upper(directive)
				case '3D'
					ax = obj.plot3D(varargin{:});
				otherwise
					ax = obj.plot2D(varargin{:});
			end
        end
        
        function ax = plot_residuals(M, varargin)
            [ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
            set(ax.Parent, 'name', 'residuals');
            
            ind = M.time2inds(tt);
            data = squeeze(M.Data(ind, M.locs));
            
            % fit a plane
            xx = M.Position(:, 1); yy = M.Position(:, 2);
            zfit = M.Beta(ind, 1) + M.Beta(ind, 2) * xx + M.Beta(ind, 3) * yy;
            
            % Draw a scatter of the data
            dot_size = M.get_dot_size(ax);
			scatter(ax, xx, yy, dot_size, data(:) - zfit, 'filled', 's'); 
            ax.CLim = [-1 1] * quantile(abs(data(:) - zfit), .95);
            colormap(ax, make_diverging_colormap);
            
			title(ax, ["Residuals" ...
                sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', ...
                tt, M.p(ind), M.Magnitude(ind), 'mm/s')]); 
			hold(ax, 'off');
            cb = colorbar;
            title(cb, 'Resid');
            
            xlim(ax, [0 M.GridSize(1) + 1]);
            ylim(ax, [0 M.GridSize(2) + 1]);
			
        end
        
		function ax = plot3D(M, varargin)
			% ax = plot(obj, t0, ax); % if object has only one time point, t0=obj.Time
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
            t_ind = M.time2inds(tt);
            
            % Get the data and x-/y-axis positions
            data = squeeze(M.Data(t_ind, :, :));
            betas = M.Beta(t_ind, :);
			N = size(data);
            pos = M.Position;
            
            % Draw a scatter of the data
			scatter3(ax, pos(:, 1), pos(:, 2), data(M.locs), [], data(M.locs), 'filled'); 
            
            % Show the fitted plane
			hold(ax, 'on');
			[xx, yy] = ndgrid(1:N(1), 1:N(2));
			zfit = betas(1) + betas(2) * xx + betas(3) * yy;
			im = surf(ax, xx, yy, min(data(:)) * ones(N), data, 'linestyle', 'none'); 
			im.Tag = 'Data';
			surfax = surf(ax, xx, yy, zfit, 'linestyle', 'none', 'facealpha', .5);
			surfax.Tag = 'Fit';
			title(ax, sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', tt, M.p(t_ind), M.Magnitude(t_ind), 'mm/s')); 
			hold(ax, 'off');
            ax.Tag = checkname(fullfile('figs', ...
                sprintf('%s_%s_%d_plot3D', M.Name, class(M), tt)));
			
        end
        
		function [ax, tt] = plot2D(s, varargin)
			% ax = plot2D(obj, tt, ax)
			[ax, tt, directive] = WaveProp.parse_plot_inputs(varargin{:});
            if ndims(s.Data) == 3
				[~, which_t] = min(abs(tt - s.Time));
                tt = s.Time(which_t);
				data = squeeze(s.Data(which_t, :, :));
				vx = s.Vx(which_t);
				vy = s.Vy(which_t);
                pval = s.p(which_t);
                
            else
				data = s.Data;
				vx = s.Vx;
				vy = s.Vy;
                pval = s.p;
            end
            phi = atan2(vy, vx);
            speed = norm([vx vy]) * s.MMpElectrode;
            
            axis(ax, 'square');
            
            [xx, yy] = find(isfinite(data));
            val = data(isfinite(data)) * 1e3;  % convert from s to ms
            units = 'ms';
            
            % Show the time data
			sc = scatter(ax, xx, yy, [], val, ...
                'filled', 's', directive{:});
            set(ax, 'clim', quantile(sc.CData, [.05, .95]));
			set(ax, 'xtick', [], 'ytick', []);
            
            % Add the arrow
            hold(ax, 'on')
            grid_size = size(data);
            sp = sqrt(vx^2 + vy^2) / (.8*max(grid_size)); % Make the arrow 8 units long
            vx = vx/sp; vy = vy/sp;
            grid_center = grid_size/2 + 1;
            X = grid_center - [vx, vy]/2;  % center the arrow on the array
            qq = quiver(ax, X(1), X(2), vx, vy, s.scale_quiver, ...
                'LineWidth', 2, ...
                'MaxHeadSize', .5, 'color', [0 0 0]);
            
            
            hold(ax, 'off')
            ax.Tag = checkname(sprintf("figs/%s_%s_%0.3f_plot2D", s.Name, class(s), tt));
            xlim([0 grid_size(1) + 1]);
            ylim([0 grid_size(2) + 1]);
            sc.SizeData = s.get_dot_size(ax);
            axis(ax, 'equal');
            cb = colorbar;
            title(cb, sprintf('TOA [%s]', units))
            title(ax, sprintf('t = %0.3f s', tt));
            lgd_string = ...
                sprintf('\\phi=%0.0f\\circ\n sp=%0.1f mm/s\n p=%0.4g', rad2deg(phi), speed, pval);
            legend(qq, lgd_string)
        end
                        		
    end
    
		
	methods (Static)
        
        
        
        function dot_size = get_dot_size(ax)
            % Get a good dot size to represent an electrode for a scatter
            % plot.
            uu = ax.Units;
            ax.Units = 'points';
            pos = ax.Position;
            ax.Units = uu;
            gs = range([ax.XLim; ax.YLim]');
            dot_size = min(pos(3:4))/max(gs - 1);
            dot_size = dot_size^2;
        end
        
		function [ax, tt, directive] = parse_plot_inputs(ax, tt, varargin)
            directive = {}; 
            deft_tt = 0;
            deft_ax = @() gca;  % axes(figure);  % just made 70 figures because I forgot to give an axis... probably going to crash matlab
            switch nargin
                case 0
                    ax = deft_ax(); tt = deft_tt; % default to new figure
                case 1
                    if isa(ax, 'matlab.graphics.axis.Axes'), tt = deft_tt; 
                    else, tt = ax; ax = deft_ax();
                    end
                otherwise
                    if isa(ax, 'matlab.graphics.axis.Axes') 
                        if isnumeric(tt), directive = varargin;
                        else, directive = [tt, varargin]; tt = deft_tt; 
                        end
                    elseif isnumeric(ax)
                        if isa(tt, 'matlab.graphics.axis.Axes')
                            directive = varargin;
                            tt_ = tt; tt = ax; ax = tt_;
                        else
                            directive = [tt varargin];
                            tt = deft_tt;
                        end
                    else
                        ax = deft_ax(); tt = deft_tt; directive = varargin;
                    end
            end
			
        end
        
        function [V, p, beta] = fit_data(data, position)
            % Wrapper for .fit_plane()
            [V, p, beta] = WaveProp.fit_plane(data, position);
        end
        
		function [V, p, beta, stats] = fit_plane(data, position)
			% Fit a 2D plane to the data. Return nan if not enough data or
			% beta are zero.
			V = nan(1, 2); p = nan; beta = nan(1, 3); stats = nan;
            if numfinite(data(:)) < 4, return; end
            
            % fit the delay vs two-dimensional positions
            W = warning;
            warning('off', 'stats:statrobustfit:IterationLimit');
            warning('off', 'stats:robustfit:RankDeficient');
            warning('off', 'MATLAB:singularMatrix');
            [beta, stats] = robustfit(position, data, 'fair'); 
            
            % perform F test that last two coefs are both 0.
            H = [0 1 0; 0 0 1];  % These are defaults, but here for clarity
            c = [0 ; 0];
            p = linhyptest(beta, stats.covb, c, H, stats.dfe);  
            
			warning(W);
            
			% beta is invalid if nan or if the slope is 0 in both directions
			invalid = all(beta(2:3).^2 < eps);
			if invalid, return; end			

			V = pinv(beta(2:3));
        end
        
		function [reduced_data, G] = use_largest_cluster(data, position)
			% Use the largest cluster of data points for lfp methods
            % Assumes first dimension is time; if 2D, <position> must be
            % given.
            
            if nargin < 2
                assert(ndims(data) == 3, 'Data must be 3d or electrode position must be given.')
                [px, py] = find(ones(size(data, [2 3])));
                position = [px, py];
            end
            
            reduced_data = nan(size(data));
            G = nan(size(data));
			warning('off', 'stats:clusterdata:MissingDataRemoved');
            
            for ii = 1:size(data, 1)
                temp_dat = data(ii, :);
                X = [position, normalize(temp_dat(:))];
                if numfinite(X(:, end)) < 3, continue; end

                cc = clusterdata(X, ...
                    'distance', 'euclidean', ... % euclidean distance
                    'criterion', 'distance', ...
                    'linkage', 'single', ...   % Nearest neighbor
                    'cutoff', sqrt(2 + .5));  % ... must be ≈connected and values differ by less than .5 sd
                [~, largest] = max(histcounts(cc, max(cc)));
                mask = cc == largest;
%                 reduced_data = data;
                temp_dat(~mask) = nan;
                reduced_data(ii, :) = temp_dat;
                G(ii, :) = cc;
            end
			
        end
        
		function S = resize_obj(s, to_struct)
			if nargin < 2, to_struct = false;  end
			
            if numel(s) == 1; S = s; return; end
             % If all positions match, continue (if not, do not combine
             % these objects)
            temp_position = cat(3, s.Position);
            assert(all(temp_position(:, :, 1) == temp_position(:, :, 1:end), 'all'), ...
                'Field <Position> is not the same for all objects');
            
			S = copy(s(1));
            
            if to_struct
				S = struct;
				for p = string(properties(s(1))')
					S.(p) = []; 
				end
            end
            
            toa_new = cat(1, s.TOA);
            S.TOA = toa_new(:, :);
            S.t0 = cat(1, s.t0);            
            S.Inds = [];
            
            mco = metaclass(s).PropertyList;
            FF = {mco.Name};
            FF(cat(1, mco.Dependent) | cat(1, mco.Hidden) | cat(1, mco.Transient)) = [];
%             out = struct;
            for f = string(FF)
                if ismember(lower(f), ["toa", "t0", "position"]), continue; end
                new_val = cat(1, s.(f));
                siz = size(new_val);
                % Get unique values
                if isvector(new_val), uu = unique(new_val);
                else, uu = unique(new_val(:, :), 'rows');
                end
                
                if size(uu, 1) == 1
                    S.(f) = reshape(uu, [1, siz(2:end)]);
                else
                    S.(f) = new_val;
                end
            end
            

			
        end
                
	end
	
end

