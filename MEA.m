classdef MEA < handle & matlab.mixin.Copyable

    properties 
        Patient
        Seizure
        LFP
        EventTimes
        Duration
        Padding
        Map
        Position
        Units = '0.25 microvolts'
        GridSize = [10 10]  % unless it's a sim
        AcquisitionSystem
        BadChannels
        SamplingRate
                
    end
    
    properties (Transient = true, Dependent = true)
        locs
        Name
    end
    
    properties (Transient = true)
        
        
        Time
        FiringRate
        
        FrWindow (1, 1) double {mustBeNonnegative} = 20 % ms
        
        MaxDescentData
        IsSim = false
        
        
    end
    
    
    methods % Constructor and getters
        
        function mea = MEA(path)
            
            % Load the requested struct
            if nargin < 1 || isempty(path)
                path = 'Data/P4_Seizure1.mat'; 
            end
            struct = load(path);
            
            % Pass structure values to mea object
            for ff = string(fieldnames(struct)')
                if ismember(ff, {'Name'}), continue; end
                mea.(ff) = struct.(ff);
            end
            
            % Flag sims
            if ismember(mea.Patient, {'SIM', 'SCM', 'FS', 'SW'})
                mea.IsSim = true;
            end
            
            % Compute time vector from sample rate and padding
            mea.Time = (0:size(mea.LFP, 1) - 1) / mea.SamplingRate - mea.Padding(1);
            
        end

        
        function data = get.LFP(mea)
            data = mea.LFP;
            data(:, mea.BadChannels) = [];
        end
        
        
        function et = get.EventTimes(mea)
            et = mea.EventTimes;
            et(mea.BadChannels) = [];
        end

        
        function pos = get.Position(mea)
            pos = mea.Position;
            pos(mea.BadChannels, :) = [];
        end
        
   
        function fr = get.FiringRate(mea)
            % wrapper to compute firing rate based on whether data is
            % simulated or in vivo
            if isempty(mea.FiringRate)
                if mea.IsSim
                    load(mea.Path, 'firingRate'); fr = firingRate;
                    mea.FiringRate = fr;
                else
                    fr = mea.compute_firing_rate(); 
                    mea.FiringRate = fr;
                end
            else
                fr = mea.FiringRate;
            end
        end
        
    end
    
    
    methods % Convenience
        
        
        function D = make_3d(mea, data)
            % Convenience: transform data 3-D 
            % D = mea.make_3d(data=mea.LFP)
            if nargin < 2, data = mea.LFP; end
            nT = size(data, 1);
            D = nan([nT, max(mea.Position)]);
            temp = nan(max(mea.Position));
            for ii = 1:nT
                temp(mea.locs) = data(ii, :); 
                D(ii, :, :) = temp; 
            end
        end

        
        function locs = get.locs(s)
            % convenience: converts x,y-coords to indices
            locs = sub2ind(max(s.Position), s.Position(:, 1), s.Position(:, 2));
        end
        
        
        function name = get.Name(mea)
            name = sprintf('%s_Seizure%s', mea.Patient, mea.Seizure);
        end
        
        
        function inds = time2inds(mea, times)
            % Get indices of requested time points
            inds = interp1(mea.Time, 1:length(mea.Time), times, 'nearest');
        end
        
        
    end
    
    
    methods  % Computed properties
        
        
        function wt = get_discharge_times(mea, min_electrodes, min_peak_height)
            % [wt, peaks] = get_discharge_times(min_electrodes=30, min_peak_height=1)
            % Finds lfp crossing times on each channel and then
            % find where at least <min_electrodes> have crossing times in a 40 ms
            % window
            
            if nargin < 2 || isempty(min_electrodes)
                if mea.IsSim, min_electrodes = 4;
                else, min_electrodes = 30; 
                end
            end
            if nargin < 3 || isempty(min_peak_height), min_peak_height = 1; end
            
            min_peak_dist = .02;  % require discharges to be at least 20 ms apart
            
            % Clip data so that crossing times become peaks
            data = -normalize(mea.LFP);
            data = min(data, min_peak_height*1.01); % clip so that findpeaks returns threshold crossing times
                
            % Identify crossing times on each channel
            [~, lcs] = arrayfun(@(ii) ...
                findpeaks(data(:, ii), ...
                    'MinPeakHeight', min_peak_height, ...
                    'MinPeakDistance', min_peak_dist * mea.SamplingRate), ...  
                1:size(data, 2), 'uni', 0);
            
            % Convert crossing times to sparse matrix and count the number
            % of channels with crossing times in sliding 40 ms windows
            ch = arrayfun(@(ii) ii * ones(size(lcs{ii})), 1:size(data, 2), 'uni', 0);
            ch = cat(1, ch{:});
            lcs = cat(1, lcs{:});
            [m, n] = size(data);
            N = sum(movmax(sparse(lcs, ch, 1, m, n), .04 * mea.SamplingRate), 2);
            
            % Return time points where crossing times are detected on at least
            % <min_electrodes> and at least 20 ms apart
            [~, wt] = findpeaks(N, mea.Time, ...
                'minpeakheight', min_electrodes, ...
                'minpeakdistance', min_peak_dist);  % again, impose max discharge rate

        end
        

        function firing_rate = compute_firing_rate(mea)
            % Convert event times to firing rate traces
            
            % Convert event times on each channel to time indices
            event_inds = mea.EventTimes;
            for ii = 1:numel(event_inds)
                event_inds{ii} = ...
                    round((event_inds{ii} + mea.Padding(1)) * mea.SamplingRate) + 1;
            end
            
            % Create a binary matrix from event indices
            firing_rate = zeros(size(mea.LFP));
            for ii = 1:size(firing_rate, 2)
                firing_rate(event_inds{ii}, ii) = 1;
            end
            
            % Compute sliding window firing rate
            firing_rate = smoothdata(firing_rate, 1, 'movmean', mea.FrWindow);
            firing_rate = firing_rate * mea.SamplingRate;  % convert to spikes per second
            
        end

    end
    
    
    methods % Wrappers to wave fitting functions
        
        function [out, M] = get_IW_templates(mea, M, win, max_templates, method)
            % [out, M] = get_IW_templates(mea, M=::saved Miw::, win=4, max_templates=Inf, method)
            % Use max descent type analysis to find IW times and templates
            % Wrapper for method defined in IW.m
            
            iw = IW(mea); 
            
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
            iw = IW(mea); 
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
        
    end
        
    
    methods  % Plotting
        
        function mov = preview(mea, times, ax, varargin)
            % mov = preview(times, ax=axes(figure), ::imagesc directives:: )
            % Does a drawnow playback of LFP frames
            if nargin < 3 || isempty(ax), ax = axes(figure); end
            t_inds = mea.time2inds(times);
            N = length(t_inds);
            im = imagesc(ax, nan(max(mea.Position)), varargin{:});
            axis(ax, 'equal')
            ttl = title(ax, '0');
            mov(N) = getframe(ax);
            for ii = 1:N
                im.CData(mea.locs) = mea.LFP(t_inds(ii), :);
                set(ttl, 'string', num2str(times(ii)));
                drawnow;
                mov(ii) = getframe(ax);
            end
        end

        
        function plot_panels(mea, times, data, layout, h, style)
        % mea.plot_panels(times, data=mea.LFP, layout=[], h=gcf, style='')
        % Plots the LFP at <times> as a series of frames
        %   <style> can be 'scatter' or 'raw'
        
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
            if isempty(data), data = zscore(mea.LFP); end
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
                
        
        function im = imagesc(mea, t0, data, ax, varargin)
            % imagesc(mea, t0=0, data=mea.lfp, ax=gca, ::imagesc directives::)
            % Shows the data at frame at <t0>
            
            if nargin < 4, ax = gca; end
            if nargin < 3 || isempty(data), data = mea.LFP; end
            if nargin < 2 || isempty(t0), t0 = 0; end
            
            temp = nan(max(mea.Position));
            ind = mea.time2inds(t0);
            temp(mea.locs) = data(ind, :);
            im = imagesc(ax, temp, varargin{:});
        end
        
    end
    
    
    methods (Static)
        
        
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

    
end



