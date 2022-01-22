classdef MEA < handle & matlab.mixin.Copyable

	properties 
        Patient char
        Seizure char
		Data
		Duration
        Padding
		Map
		Position
        Units = '0.25 microvolts'
        AcquisitionSystem
        BadChannels
        
        Path = 'Data/P4_Seizure1.mat'
        
	end
	
	properties (Transient = true)
		locs
		Raw
		AllTime
		MaxDescentData
        IsSim = 0
        GridSize
        Name
        
        SamplingRate = 1000
        Time
        
        FrWindow (1, 1) double {mustBeNonnegative} = 20 % ms
        ArtefactThresh (1, 1) double {mustBeNonnegative} = 16  % stdev
        MinEventSeparation (1, 1) double {mustBeNonnegative} = 1  % ms
        EventThresh (1, 1) double {mustBeNonnegative} = 4  % stdev
        
		SRO
		lfp
		mua
		artefacts
		firing_rate
        iw_firing_rate
		mua_events
		wave_fits
        event_inds
        
        Fits
        IW
        
	end
	
	properties (Access = public, Dependent = true)
		skipfactor
		Center
		Active
    end
    
    
	methods 
        
		function mea = MEA(path, SR)
			
			if nargin == 2, mea.SamplingRate = SR; end
            if nargin < 1, path = mea.Path; end
            if isnumeric(path)
                sz = SeizureInfo;
                path = sz.fname{path};
                disp(path)
            end
            path = dir(path);
            mea.Path = path(1);
			mea.load;
		end
		
		function mea = load(mea)
			temp = load(mea.Path);
			mea.SRO = temp.SamplingRate;
            mea.Raw = temp.Data;
            mea.AllTime = temp.Time;
			
            for f = fieldnames(temp)'
                switch f{:}
                    case {'Data', 'SamplingRate', 'Time', 'Path'}
                        continue
                    case fieldnames(mea)
                        mea.(f{:}) = temp.(f{:});
                    otherwise
                        continue
                end
            end

        end
        
        function mea_struct = save(mea, path)
            % Saves to a struct rather than to an MEA object
            if nargin < 2, path = pwd; end
            
            bad_chan = mea.BadChannels;
            mea.BadChannels = [];
            
            mea_struct = struct( ...
                'SamplingRate', mea.SRO, ...
                'Time', single(mea.AllTime), ...
                'BadChannels', bad_chan, ...
                'Data', mea.Raw ...
                );
            
            for ff = ["Patient", "Seizure", "Duration", "Padding", "Map", ...
                    "Position", "Units", "GridSize", "AcquisitionSystem"]
                mea_struct.(ff) = mea.(ff);
            end
            
            save(fullfile(path, mea.Name), '-v7.3', '-struct', 'mea_struct')
            
        end
        
        
        function [bad_ch, bad_ch0] = get_bad_channels(mea, show_pc_coords)
            % Remove bad channels. Detect channels that are outliers in PC
            % space (using 95% of variance explained). Outliers := > 3 MAD
            % from PC center of mass
            
            if nargin < 2 || isempty(show_pc_coords), show_pc_coords = false; end
            
            S = warning; warning('off');
            bad_ch0 = mea.BadChannels;
            mea.BadChannels = [];
            data = single(mea.Data);
            sd = std(data)';
            outsSD = isoutlier(sd, 'ThresholdFactor', 3);  % remove outliers by std; (threshold = 3 MAD)
            fprintf('outs by sd: ')
            fprintf('%d ', find(outsSD));
            fprintf('\n');
            
            % Get the largest contiguous cluster in PC space
            outsPC = false(size(outsSD));
            [coeff, ~, ~, ~, explained] = pca(data(:, ~outsSD));
            num_pc = min(find(cumsum(explained) >= 95, 1), 3);

            T2 = clusterdata(coeff(:, 1:num_pc), 'criterion', 'distance', ...
                'distance', 'mahalanobis', 'cutoff', sqrt(2));
            if show_pc_coords
                scatter3(coeff(:, 1), coeff(:, 2), coeff(:, 3), [], T2)
            end
            [~, largest_group] = max(histcounts(T2, 'BinMethod', 'integers'));
            outsPC(~outsSD) = T2 ~= largest_group;
            fprintf('outs by pca: ')
            fprintf('%d ', find(outsPC));
            fprintf('\n');
            
            outsF = outsSD | outsPC;
            bad_ch = reshape(find(outsF), [], 1);
            
            warning(S);
            mea.BadChannels = bad_ch0;  % Preserve mea.BadChannels
        end
        
        
        function iw = get.IW(mea)
            if isempty(mea.IW)
                mea.IW = IW(mea); %#ok<CPROP>
            end
            iw = mea.IW;
        end
        
        function gs = get.GridSize(mea)
            gs = size(mea.Map);
        end
        
		function active = get.Active(mea)
			active = (mea.Time > 0) & (mea.Time < mea.Time(end) - mea.Padding(2));
        end
        
        
		function center = get.Center(mea)
			[~, center] = min(sum((mea.Position - mean(mea.Position)).^2, 2));
        end
                
                		
		function raw = get.Raw(mea)
			
            raw = mea.Raw;
			raw(:, mea.BadChannels) = [];
        end
        
        
		function set.Path(mea, value)
            if isstruct(value)  % allow entry of string or dir struct
                file = value;
            else
                file = dir(value);
            end
			mea.Path = fullfile(file.folder, file.name);
        end
		
                
		function time = get.Time(mea)
			time = downsample(mea.AllTime, mea.skipfactor);
        end
        
        
		function time = get.AllTime(mea)
			time = mea.AllTime();  % this is saved as a function to save space
        end
		
        
		function set.SamplingRate(mea, value)
			mea.SamplingRate = value;
        end
        
        
		function SR = get.SamplingRate(mea)
			SR = min(max(mea.SamplingRate, 1), mea.SRO);
        end
        
        
		function set.skipfactor(mea, value)
			if isempty(value), return, end
			mea.skipfactor = max(min(value, mea.SRO), 1);
			mea.SamplingRate = mea.SRO / mea.skipfactor;
        end
		
        
		function skipfactor = get.skipfactor(mea)
			skipfactor = round(mea.SRO / mea.SamplingRate);
        end

        
		function data = get.Data(mea)
			data = single(resample(double(mea.Raw), 1, mea.skipfactor));
        end

        
        function D = make_3d(mea, data)
        % D = mea.make_3d(data=mea.Data)
            if nargin < 2, data = mea.Data; end
            nT = size(data, 1);
            D = nan([nT, max(mea.Position)]);
            temp = nan(max(mea.Position));
            for ii = 1:nT
                temp(mea.locs) = data(ii, :); 
                D(ii, :, :) = temp; 
            end
        end
		
        
		function lfp = get.lfp(mea)
            % mea data bandpass filtered to [1 50] Hz
			lfp = mea.lfp;
			if isempty(lfp)
				fprintf('Filtering lfp ... ')
				lfp = mea.filter(mea.Data, mea.SamplingRate, [1 50]); 
				mea.lfp = lfp;
				fprintf('Done. \n')
			end
        end

        
		function artefacts = get.artefacts(mea)
			if isempty(mea.artefacts), mea.mua; end
			artefacts = mea.artefacts;
        end
        
        
		function mua = get.mua(mea)
			mua = mea.mua;
			if isempty(mua)

                data = single(mea.Raw) - nanmean(mea.Raw, 2);  % perform common average referencing before filtering to reduce sensor noise
				disp('Filtering mua ...')
				mua = mea.filter(data, mea.SRO, [300 3000]);
				
				disp('Done')
				mea.artefacts = mea.get_artefacts(mua, mea.ArtefactThresh);
				mua(mea.artefacts) = nan;
				mea.mua = mua;
			end
        end
        
        
		function EI = get.event_inds(mea)
			EI = mea.event_inds;
			
			if isempty(EI)
                mua_ = filloutliers(mea.mua, 'linear', 'thresholdfactor', 10);  % remove electrode noise
                
                mua_ = normalize(mua_, 2, 'center');  % diminish the effect of simultaneous events/mea noise
                data = normalize(mua_, 'scale');  % Identify activity on each channel
				intervalM = mea.SRO / 1e3 * mea.MinEventSeparation;  % samples per MIN_DIST

				% find events on each channel
				events = false(size(data));
                
				for ch = 1:size(data, 2)
                    temp = fillmissing(data(:, ch), 'linear');
                    if ~any(-temp > mea.EventThresh)
						continue
                    end
					[~, inds] = findpeaks(-temp, ...
						'minpeakdistance', intervalM, ...
						'minpeakheight', mea.EventThresh);  
					events(inds, ch) = true;
				end
				EI = sparse(events);
				mea.event_inds = EI;
			elseif ~issparse(EI)
				sz = size(mea.Raw);
				[i, j] = ind2sub(sz, EI);
				EI = sparse(i, j, true, sz(1), sz(2));
			end
			

        end
		
        
		function events = get.mua_events(mea)
			events = mea.event_inds;
        end
		
		
		function pos = get.Position(mea)
			pos = mea.Position;
			pos(mea.BadChannels, :) = [];
        end
		
   
		function fr = get.firing_rate(mea)
			if isempty(mea.firing_rate)
                if mea.IsSim, fr = -mea.Data; 
                else
                    fr = mea.compute_firing_rate(); 
                    mea.firing_rate = fr;
                end
			else
				fr = mea.firing_rate;
			end
        end
        
        function name = get.Name(mea)
            name = sprintf('%s_Seizure%s', mea.Patient, mea.Seizure);
        end
        
		function [inds, dist] = time2inds(mea, times)
			inds = interp1(mea.Time, 1:length(mea.Time), times, 'nearest', 'extrap');
			dist = times - mea.Time(inds);
        end
        
        
		function locs = get.locs(s)
			locs = sub2ind(max(s.Position), s.Position(:, 1), s.Position(:, 2));
        end
        
        
		function inds = which_t(s, t0)
			[~, inds] = min(abs(s.Time - t0));
		end
		
    end
	
    
	methods % Wrappers for wave fitting and IW detection functions
        
        function [out, M] = get_IW_templates(mea, M, win, max_templates, method)
            % [out, M] = get_IW_templates(mea, M=::saved Miw::, win=4, max_templates=Inf, method)
            % Use max descent type analysis to find IW times and templates
            % Wrapper for method defined in IW.m
            
            iw = IW(mea); %#ok<CPROPLC>
            
            if nargin < 5 || isempty(method), method = iw.Method; end
            if nargin < 4 || isempty(max_templates), max_templates = iw.MaxTemplates; end
            if nargin < 3 || isempty(win), win = iw.W; end
            if nargin < 2 || isempty(M); M = ''; end
            
            iw.MaxTemplates = max_templates;
            iw.W = win;
            iw.Method = method;
            
            [out, M] = iw.get_IW_templates(M);
            
        end
        
        
        function M = max_descent_IW(mea, type, save_bool)
            if nargin < 2 || isempty(type), type = 'new'; end
            if nargin < 3, save_bool = false; end
            iw = IW(mea); %#ok<CPROPLC>
            M = iw.save_IW_fits(type, save_bool);
        end
        
                
		function M = max_descent(mea, times, varargin)
            if nargin < 2 || isempty(times) 
                times = mea.get_discharge_times;
            end
            
            M = MaxDescent(mea, times, varargin{:});
			M.Name = mea.Name;
        end
        
        
		function D = delays(mea, times, varargin)
			
            D = Delays(mea, times, varargin{:});

        end
        
        
        function wt = get_discharge_times(mea, method, min_electrodes, min_peak_height)
			% [wt, peaks] = get_discharge_times(method='lfp_cross', min_electrodes=30, min_peak_height=1)
            % Finds peaks on each channel according to method and then
            % finds where at least <min_electrodes> have peaks in an
            % appropriate window (20 ms in most cases since this is 50Hz,
            % which is the upper bound of the lfp filter)
            
			if nargin < 2, method = 'lfp_cross'; end
            if nargin < 3 || isempty(min_electrodes)
                if mea.IsSim, min_electrodes = 4;
                else, min_electrodes = 30; 
                end
            end
            if nargin < 4 || isempty(min_peak_height), min_peak_height = 1; end
            method = validatestring(method, {'lfp', 'lfp_cross', 'events'});
            
            min_peak_dist = .02;
            switch method
				
				case 'lfp'
					data = -normalize(mea.lfp);
                    
                case 'lfp_cross'
                    data = -normalize(mea.lfp);
                    data = min(data, min_peak_height*1.01); % clip so that findpeaks returns threshold crossing times
                    
				case 'events'
					data = normalize(mea.firing_rate, 'scale');
                    if min_peak_height < 2
                        warning('min_peak_height=%0.1f; suggest using at least 2.', min_peak_height); 
                    end
				otherwise
					error('Method not recognized')
            end
            
            [~, lcs] = arrayfun(@(ii) ...
                findpeaks(data(:, ii), ...
                    'MinPeakHeight', min_peak_height, ...
                    'MinPeakDistance', min_peak_dist * mea.SamplingRate), ...  % require peaks at least 20 ms apart (this gives a max 50hz discharge rate)
                1:size(data, 2), 'uni', 0);
            ch = arrayfun(@(ii) ii * ones(size(lcs{ii})), 1:size(data, 2), 'uni', 0);
            ch = cat(1, ch{:});
            lcs = cat(1, lcs{:});
            [m, n] = size(data);
            N = sum(movmax(sparse(lcs, ch, 1, m, n), .04 * mea.SamplingRate), 2);
            
            [~, wt] = findpeaks(N, mea.Time, ...
                'minpeakheight', min_electrodes, ...
                'minpeakdistance', min_peak_dist);  % again, impose max 50 Hz firing rate

        end
        
        
        function firing_rate = compute_firing_rate(mea)
			
			samplingRateMs = mea.SRO / 1e3;  % samples per ms
			window = mea.FrWindow * samplingRateMs;  % number of samples to use in the window

			% Compute spike rate
			firing_rate = downsample(...
					smoothdata(double(full(mea.mua_events)), 1, 'movmean', window), ...  % mean spikes per sample
                    mea.skipfactor ...
				) * mea.SRO;  % convert to spikes per second 
		end
        
        
        function [S, t, f] = specgram(mea)
            % mostly just here as a reminder of how to use mtspecgramc
            params_.tapers = [3 5 1]; % [T W p]; ntapers = TW-p
            params_.fpass = [0 100];
            params_.Fs = mea.SamplingRate;
            params_.trialave = 1;
            win = 5;
            winstep = 1;
            [S, t, f] = mtspecgramc(mea.Data, [win winstep], params_);
        end
        
        
		function [S, fs] = pspectrum(mea, data, fs, range)
			% [S, fs] = pspectrum(data=Data, fs=SamplingRate, range=[0 100])
			if nargin < 4, range = [0 100]; end
			if nargin < 3, fs = mea.SamplingRate; end
			if nargin < 2, data = mea.Data; end
			
			[S, fs] = pspectrum(data, fs);
			S = S(fs > range(1) & fs < range(2), :);
			fs = fs(fs > range(1) & fs < range(2));
			S = pow2db(S.*fs.^2);
        end        
        
	end
		
	
	methods (Static)
        
        function artefacts = get_artefacts(data, thresh)
			artefacts = (abs(zscore(data)) > thresh); 
        end
        
        
		function [data, band] = filter(data, sampling_rate, band)
            % data = filter(data, sampling_rate, band, use_high_filter_order=false);  % (static)
            % Filters <data> to <band> using 150-order bandpass FIR filter
        
            if band(2) > sampling_rate / 2 / 1.1
				band(2) = (sampling_rate / 2 - 1) / 1.1;
				fprintf('Setting CutoffFrequency2 to %.1f\n', band(2)); 
            end

            if band(1) == 0, band(1) = 1; end  % Neuroports pre-filter to [.3-10k] Hz
			
				bpFilt = designfilt('bandpassfir', ...
                    'FilterOrder',150, ...
					'CutoffFrequency1',band(1), ...
                    'Cutofffrequency2',band(2), ...
					'SampleRate', sampling_rate);
                
			data = single(filtfilt(bpFilt, double(data)));
        end
        
		
    end
	
    
	methods  % Plotting
        
        function mov = preview(mea, times, ax, varargin)
			% mov = preview(times, ax=axes(figure), ::imagesc directives:: )
			if nargin < 3 || isempty(ax), ax = axes(figure); end
			t_inds = mea.time2inds(times);
			N = length(t_inds);
			im = imagesc(ax, nan(max(mea.Position)), varargin{:});
			ttl = title(ax, '0');
			mov(N) = getframe(ax);
			for ii = 1:N
				im.CData(mea.locs) = mea.Data(t_inds(ii), :);
				set(ttl, 'string', num2str(times(ii)));
				drawnow;
				mov(ii) = getframe(ax);
			end
        end

        
		function plot_panels(mea, times, data, layout, h, style)
        % mea.plot_panels(times, data=mea.Data, layout=[], h=gcf, style='')
			if nargin < 6 || isempty(style), style = ''; end
            if nargin < 5 || isempty(h), h = gcf; end
			if nargin < 4 , layout = []; end
            if nargin < 3, data = []; end
			
            style = validatestring(style, {'', 'scatter', 'raw'});
			set(0, 'currentfigure', h);
            clf(h);
            if isempty(layout)
                T = tiledlayout(h, 'flow');
            else
                T = tiledlayout(h, layout(1), layout(2));
            end

            T.Tag = 'panels';
            
			inds = mea.time2inds(times);
			if isempty(data), data = zscore(mea.Data); end
            data = data(inds, :);
            
			temp = nan(max(mea.Position));
            
            % prep plotting params
            ax = nexttile(T, 1);
            axis(ax, 'square');
            units = ax.Units;
            set(ax, 'units', 'points');
            pos = ax.Position;
            ax.Units = units;
            pt_width = max(min(pos([3 4])) / 11, 1);
                
            % prep scatter params
            switch style
                case {'scatter', ''}
                    siz_min = max(.3 * pt_width, 1);  % 1.733);
                    siz_max = (2.5 * pt_width);
                    sd_min = 2;  % show range of 2 sd from mean ...
                    sd_max = max(quantile(abs(data), .999, 'all'), 4);  % to at least 4 sd from mean
                    siz_dat = rescale(abs(data), sqrt(siz_min), sqrt(siz_max), ...
                        'inputmin', sd_min, ...
                        'inputmax', sd_max).^2;
                    col_dat = single(data > 0) - single(data < 0);
                    cmap = lines;
                    clim_ = @(dat) [-2 2];
                    
                    % Don't show very small values
                    mask = abs(data) < sd_min;
                    data(mask) = nan;
                    siz_dat(mask) = nan;
                    col_dat(mask) = nan;
                    
                    
                    
                case 'raw'
                    col_dat = data;
                    siz_dat = pt_width.^2 + zeros(size(data));
                    cmap = gray;
                    cmap(cmap(:, 1) > .85, :) = [];  % don't go all the way to white
                    clim_ = @(dat) quantile(dat, [.1 .95]);
            end
            
            % start plotting
			for ii = 1:length(inds)
                ax = nexttile(T, ii);
                scatter( ...
                    mea.Position(:, 1), mea.Position(:, 2), ...  % x, y
                    siz_dat(ii, :), ...  % size
                    col_dat(ii, :), ...  % color
                    'filled', 's');
                axis square
                colormap(cmap)
                set(gca, 'clim', clim_(col_dat(ii, :)));
                switch style
                    case 'raw0'
                        temp(mea.locs) = data(ii, :);
                        imagesc(temp, [-1 1]);
                        colormap gray
                end
                xticks([]);
                yticks([]);
                
                set(ax, 'xlim', [0 11], 'ylim', [0 11], 'box', 'on');
                lgd = legend(ax, num2str(times(ii)));
                set(lgd, 'visible', 'off');
				if ii == 1
					title(num2str(times(ii)));
				end
			end
			ttl = title(num2str(times(ii)));
			lbl = xlabel(num2str(times(ii)));
			ttl.Position(2:3) = lbl.Position(2:3);
			delete(lbl);
			ttl.VerticalAlignment = 'top';
        end
        
        
		function fr = plot(s, ax)
			if nargin < 2, ax = axes(figure('units', 'normalized', 'position', [0 .5 1 .5])); end
			cmap = ax.ColorOrder;
			set(ax, 'nextplot', 'replaceChildren', 'xlim', [s.Time(1), s.Time(end)]);
			frline = gobjects(3, 1);
			NP = ax.NextPlot;
			
			
			frline(1) = plot(ax, s.Time, zscore(mean(s.lfp, 2)), ...
				'color', (cmap(1, :) + .8)/2, ...
				'linewidth', 2, ...
				'displayname', 'lfp');
            ax.NextPlot = 'add';
			mn_fr = nanmean(s.firing_rate, 2);
			frline(2) = plot(ax, s.Time, mn_fr ./ nanstd(mn_fr) - pi, ...
				'color', (cmap(2, :) + .8)/2, ...
				'linewidth', 2, ...
				'displayname', 'FR');
			
			ylim(ax, [-pi pi]);
			yticks(ax, -pi:pi/4:pi)
			yticklabels(ax, {'-\pi', '', '-\pi/2', '', '0', ...
				'', '\pi/2', '', '\pi'});
			ylabel(ax, {'Direction', 'Normalized Signal'})
			xlabel(ax, 'Time (s)');
			title(ax, [s.Patient ' ' s.Seizure]);
			ax.NextPlot = NP;
			fr.ax = ax;
			fr.frline = frline;
			fr.name = ['figs/' s.Patient '_' s.Seizure '_plot'];
        end
        
        
		function imagesc(s, data, t0, ax)
			if nargin < 4, ax = gca; end
			if nargin < 3 || isempty(t0), t0 = 0; end
			if nargin < 2 || isempty(data), data = mea.Data; end
			temp = nan(max(s.Position));
			ind = s.time2inds(t0);
			temp(s.locs) = data(ind, :);
			imagesc(ax, temp);
		end
	end
	
end



