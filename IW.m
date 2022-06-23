classdef IW < handle
    
    properties
        
        name
        wave = 1
        MinPeakHeight = 2  % Defines the threshold for candidate IW crossing times in standard deviations
        MinPeakFr = 30
        DiffsOrPeaks = 'peaks'
        FiringRateWin = 20  % in ms
        SmoothingArgs = {'movmean', 1}  % alternately, use {'gaussian', .2}; {method, seconds}
        
        MaxTemplates = inf
        W  % Big smoothing window for looking for IW; empty means autodetect
        MinElectrodes = 10
        
        iw_templates
    end
    
    properties (Transient = true)
        ManualOuts
        all_pks   % Used in get_stats plots
        all_locs
        all_durs

        mea
        locs
        fr
        fr_movmean
        fr_gaussian
        fr_smooth
        fr_at_peak
        iw_fwhm
    end
    
    properties (Dependent = true)
        time
        nch
        pks
        t_mask
        V
        speed
        phi
        num_waves
        center
        range
        show
        position
        
        onsets
        durs
        outliers
        GridSize
    end
    
    methods
        
        function self = IW(mea, varargin)
            for ii = 1:2:nargin-1
                field = validatestring(varargin{ii}, properties(self));
                self.(field) = varargin{ii + 1};
            end
            self.mea = mea;
            self.name = mea.Name;
            if ~isempty(self.FiringRateWin)
                self.mea.FrWindow = self.FiringRateWin;
            end

        end
        
        % basic mea properties
        function nch = get.nch(self)
            nch = size(self.mea.Raw, 2); 
        end
        
        function locs = get.locs(self) 
            locs = find(squeeze(any(self.iw_templates.template, 1)));
%             if isempty(self.locs), self.locs = self.mea.locs; end
%             locs = self.locs;
        end
        
        function gs = get.GridSize(self)
            gs = size(self.iw_templates.template, [2 3]);
        end
        
        function set.GridSize(self, value)
            self.mea.GridSize = value;
        end
                
        function mea = get.mea(self)
            if isempty(self.mea)
                self.mea = MEA(sprintf('Data/%s.mat', ...
                    self.name));
                self.mea.FrWindow = self.FiringRateWin;
            end
            mea = self.mea;
        end
        
        function fr = get.fr(self)
            if isempty(self.fr) || self.FiringRateWin ~= self.mea.FrWindow
                self.mea.FiringRate = [];
                self.mea.FrWindow = self.FiringRateWin;
                
                % Identify where firing rate is above the mean across the
                % array
                self.fr = self.mea.FiringRate;
                self.fr_smooth = [];
            end
            fr = self.fr;
        end
        
        function t = get.time(self); t = self.mea.Time; end
        
                
        %% Get IW templates
        
        function tpl = get.iw_templates(self)
            
            RECOMPUTE = ...
                isempty(self.iw_templates) ...
                || ~strcmpi(self.DiffsOrPeaks, self.iw_templates.mdorpeaks) ...
                ;
            
            if RECOMPUTE
                self.compute_IW_templates;
            end
            tpl = self.iw_templates;
        end
        
        function ind = main_wave(self)
            % Choose the wave with the highest firing. Set
            % non-participating electrodes to 0 and use the mean to compute
            % the firing.
            
            % use sum instead of mean to balance high firing and higher
            % numbers of electrodes
            fr_ = nansum(self.iw_templates.firing_rate, [2 3]);  
            [~, ind] = max(fr_);
        end
        
        function [M, fwhm] = max_descent_IW(iw, times, varargin)
            % Use max descent method to compute IW-like wave propagation
            
            mea_ = iw.mea;
            time_ = iw.time;
            
            data = iw.fr_smoothC;
            dataN = normalize(data, 'scale');
            
            S = warning; warning off;
            D0 = mea_.MaxDescentData;
            mea_.MaxDescentData = -dataN;
            M = MaxDescent(mea_, times, ...
                'halfwin', 1, 'diffsorpeaks', iw.DiffsOrPeaks, ...  % defaults
                varargin{:});

            M.TOA = M.use_largest_cluster(M.TOA, M.Position);
            M.MinFinite = 10;  % Set to 10 to accommodate slow movement relative to window
			M.Name = ['IW_' mea_.Name];
            
            warning(S);
            mea_.MaxDescentData = D0;
            
            % Get the IW width (duration) of the IW on each electrode
            dat = M.TOAcube + M.t0 - M.HalfWin;
            fwhm = zeros(size(dat));
            
            [~, locs_w, ww] = arrayfun(@(ii) ...
                findpeaks(dataN(:, ii), time_), 1:size(dataN, 2), 'uni', 0);
            
            for ii = 1:numel(ww)
                if isempty(ww{ii}), continue; end  % this happens when the sd of an electrode is 0
                li = mea_.locs(ii);
                temp_tq = dat(:, li);
                temp_t_all = locs_w{ii}(:);
                [~, nearest] = min(abs(temp_tq - temp_t_all'), [], 2);
                temp_w = ww{ii}(nearest);
                
                % don't add data where there isn't any, or where the times
                % don't match
                t_diff = temp_tq - temp_t_all(nearest);
                mask = isfinite(t_diff) & abs(t_diff) < 3e-3;  % within 3 ms
                temp_w(~mask) = nan;
                
                fwhm(:, li) = temp_w;
            end

            
            
        end

        function [M, s] = save_IW_fits(iw, save_bool)
            % Creates the IW files in WaveFits
            if nargin < 2, save_bool = true; end
            s = [];
            
            M = iw.max_descent_IW(iw.fr_smooth, ...
                'diffsorpeaks', 'peaks');
            
            
            Mname = 'Miw';
            if save_bool
                s = mkdir(['WaveFits/' iw.mea.Name]);
                save(['WaveFits/' iw.mea.Name filesep Mname], 'M');
            end
        end
        
        function win = get_W(iw, show)
            
            if nargin < 2, show = false; end
            
            % Auto-detect min peak window length: Use the width of the
            % largest peak as the window to differentiate between peaks
            if isempty(iw.W)  
                mea_ = iw.mea;
                time_ = mea_.Time;
                fr_ = iw.fr_smooth;
                frC = iw.fr_smoothC;


                % Normalize the firing rate based on std; set to zero any
                % electrodes with low firing rates
                frN = normalize(frC, 'zscore'); 
                frN(fr_ < iw.MinPeakFr) = 0;
                
                % Look for where firing rates are 2 sd above mean (the low
                % threshold makes the windows longer so there is less risk
                % of over segmenting IW events) using a 2 s smoothing
                % window (i.e. assume channels recruited within 2 s of each
                % other are part of the same IW event).
                win0 = 0;
                win_new = 2;
                while win_new / win0 > 1.05
                    win0 = win_new;

                    num_hi_fr = smoothdata(sum(frN > iw.MinPeakHeight, 2), ...
                        'gaussian', mea_.SamplingRate*win_new);
                    
                    peak_fun = @() findpeaks(normalize(num_hi_fr, 'scale'), time_, ...
                        'SortStr', 'descend', 'Npeaks', 1, ...
                        'WidthReference', 'halfprom', 'Annotate', 'extents');  % halfprom is what you want here... don't change this! (halfheight sounds like what you want but it truncates; prom is what you want)
                    [~, ~, win_new] = peak_fun();
                    if show
                        fprintf('win = %0.2f (%0.2f)\n', win_new, win_new / win0)
                    end
                    
                end
                win = win_new;
                if show, peak_fun(); end
                if isempty(win), win = 0; end
                win = max(min(win, 60), 2);
            else
                win = iw.W;
            end

        end
        
        function out = compute_IW_templates(iw, varargin)
            % [out, M] = compute_IW_templates(iw, method, varargin)
            
            max_templates = iw.MaxTemplates;  % Return at most this many templates
            win = iw.get_W;  % Use this as the min peak window (otherwise, autodetect)
            mea_ = iw.mea;
            time_ = mea_.Time(:)';
            MIN_FR = iw.MinPeakFr;  % channels must hit this firing rate to be considered in IW

            
            % Get the firing rate according to the indicated method
            fr_ = iw.fr_smooth;
            frC = iw.fr_smoothC;


            % Normalize the firing rate based on std; set to zero any
            % electrodes with low firing rates
            frN = normalize(frC, 'zscore');  
            

            % Find time of peaks on each channel
            wng = warning; warning('off');
            [~, locs_t] = arrayfun( ...
                @(ii) findpeaks(frN(:, ii), time_, ...
                'minpeakheight', iw.MinPeakHeight, 'minpeakdistance', win), ...
                1:size(fr_, 2), 'uni', 0);
            ch = arrayfun(@(ii) ii*ones(size(locs_t{ii})), 1:size(fr_, 2), 'uni', 0);
            warning(wng);

            
            % Convert peak times to locs
            locs_t = cat(2, locs_t{:});
            ch = cat(2, ch{:});
            [~, locs_i] = min(abs(time_(:) - locs_t));

            
            % Create a matrix of peaks and determine how many channels have
            % a peak in each sliding window
            [m, n] = size(fr_);
            temp = movmax(full(sparse(locs_i, ch, 1, m, n)), win * mea_.SamplingRate);
            N = sum(temp, 2);
            
            
            % Find time points where many electrodes have a peak
            [~, times, durs_, proms_] = findpeaks(N, time_, ...
                'MinPeakDistance', win, 'MinPeakHeight', iw.MinElectrodes);

            
            % Compute max descent (or min peak) on time points
            [M, fwhm] = iw.max_descent_IW(times, 'halfwin', win/2, varargin{:});
            dat = WaveProp.use_largest_cluster(M.Data) + M.time - M.HalfWin;
            [~, inds] = min(abs(time_(:) - times));
            fr_max = movmax(fr_, win * mea_.SamplingRate); 
            
            
            fr_max = fr_max(inds, :);
            fr_max(isnan(dat(:, mea_.locs)) | fr_max < MIN_FR) = nan;
            fr_maxN = normalize(fr_max, 2);
            
            
            fr_max(fr_maxN < -1) = nan;  % exclude any channels with low firing rates
            fr_max = mea_.make_3d(fr_max);
            dat(isnan(fr_max)) = nan;
            fr_max(isnan(dat)) = nan;
            fwhm(isnan(dat)) = nan;
            
            
            % Eliminate times where not enough electrodes are active
            nchan_ = sum(isfinite(dat), [2 3]);
            mask = nchan_ >= iw.MinElectrodes;

            
            % Return the highest firing templates
            if max_templates < sum(mask)
                fr_max(~mask, :, :) = 0;
                [~, so] = sort(sum(fillmissing(fr_max, 'constant', 0), [2 3]), 'descend');
                mask(so(max_templates+1:end)) = false;
            end
            
            
            out = struct( ...
                'template', dat(mask, :, :), ...
                'firing_rate', fr_max(mask, :, :), ...
                'fwhm', fwhm(mask, :, :), ...
                'time', M.Time(mask), ...
                'win', win, ...
                'nchan', nchan_(mask), ...
                'durs', durs_(mask), ...
                'proms', proms_(mask), ...
                'mdorpeaks', iw.DiffsOrPeaks, ...
                'M', M);
            
            iw.iw_templates = out;

        end
        

        %% IW stats
        
        function iw_fwhm_ = get.iw_fwhm(self)
            ind = self.wave;
            tpl = self.iw_templates;
            iw_fwhm_ = tpl.fwhm(ind, self.locs);
        end
        
        function fr = get.fr_at_peak(self)
            ind = self.wave;
            tpl = self.iw_templates;
            fr = tpl.firing_rate(ind, self.locs);
        end
        
        function onsets = get.onsets(self)
            onsets = self.iw_templates.template(self.wave, self.locs);
        end
        
        function durs = get.durs(self)
            % This used to be the duration of the peak on each channel; now
            % it is the duration of the iw wave as seen on all channels 
            durs = self.iw_templates.durs(self.wave);
        end
        
        function set.outliers(self, value)
            self.ManualOuts = value;
        end
        
        function outliers = get.outliers(self)
            outliers = ~isfinite(self.onsets);
        end
        
        function V = get.V(self)
            [V, ~] = self.regress;
        end
        
        function speed = get.speed(self)
            fit = self.wave_fit;
            speed = fit.speed;
        end
        
        function phi = get.phi(self)
            phi = atan2(self.V(2), self.V(1));
        end
        
        function num_waves = get.num_waves(self)
            num_waves = numel(self.iw_templates.time);
        end
        
        function center = get.center(self)
            center = nanmedian(self.onsets(~self.outliers));
        end
        
        function range = get.range(self)
            range = quantile(self.onsets(~self.outliers), [0 1]);
            
        end
        
        function pos = get.position(self)
            [p1, p2] = ind2sub(self.GridSize, self.locs);
            pos = [p1, p2];
        end
        
        function ap = get.all_pks(self)
            if isempty(self.all_pks)
                [pks_, locs_, durs_] = self.findpeaks([-Inf Inf]);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Visual inspection for multiple waves
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                self.all_pks = cat(1, pks_{:});
                self.all_locs = cat(2, locs_{:});
                self.all_durs = cat(2, durs_{:});
            end
            ap = self.all_pks;
        end
        
        function al = get.all_locs(self)
            if isempty(self.all_locs), self.all_pks; end
            al = self.all_locs;
        end
        
        function ad = get.all_durs(self)
            if isempty(self.all_durs), self.all_pks; end
            ad = self.all_durs;
        end
                
        function frS = fr_smoothC(self)
            % remove spatial mean (mean along dimension 2)
            frS = normalize(max(self.fr_smooth, self.MinPeakFr), 2, 'center');
        end
        
        function frS = get.fr_smooth(self)
            % Smooth firing rate according to method
            if isempty(self.fr_smooth)
                mtd = self.SmoothingArgs{1};
                win = self.SmoothingArgs{2};
                frS = smoothdata(self.fr, mtd, win, ...
                        'SamplePoints', self.time);
                self.fr_smooth = frS;
            end
            frS = self.fr_smooth;

        end
        
        function [pks, locs, fwhm] = findpeaks(self, t_bounds_)
            if nargin < 2, t_bounds_ = self.t_bounds; end
            S = warning;
            
            data = self.fr_smooth;
            tt = self.time;
            data(tt < 0 | tt > self.seizure_dur, :) = nan;
            data = normalize(data);
            thresh = self.MinPeakHeight;

            t_mask_ = tt >= t_bounds_(1) & tt <= t_bounds_(2);
            warning('off')
            [pks, locs, fwhm] = arrayfun( ...
                        @(ich) findpeaks( ...
                        data(t_mask_, ich), tt(t_mask_), ...
                        'minpeakheight', thresh, ...
                        'SortStr', 'descend'), ...
                        1:self.nch, 'uni', 0);
            warning(S)
        end
        
        function fr_max = peak_fr(self)
            
            t_on = self.onsets ...
                - self.durs/2;
            t_off = self.onsets ...
                + self.durs/2;

            fr_ = self.fr;
            
            [~, ind_on] = min(abs(self.time' - t_on));
            [~, ind_off] = min(abs(self.time' - t_off));
            
            N = size(fr_, 2);
            fr_max = nan(1, N);
            for ich = 1:N
                fr_max(ich) = max(fr_(ind_on(ich):ind_off(ich), ich));
            end

        end
        
        function fr_mn = mean_fr(self)
            
            t_on = self.onsets ...
                - self.durs/2;
            t_off = self.onsets ...
                + self.durs/2;

            fr_ = self.fr;
            
            [~, ind_on] = min(abs(self.time' - t_on));
            [~, ind_off] = min(abs(self.time' - t_off));
            
            N = size(fr_, 2);
            fr_mn = nan(1, N);
            for ich = 1:N
                fr_mn(ich) = mean(fr_(ind_on(ich):ind_off(ich), ich));
            end

        end
                
        function [V, p, v0, phi, speed] = regress(self, field)
            if nargin < 2 || isempty(field), field = 'onsets'; end
            y = self.(field)(~self.outliers)';
            
            X = [ones(size(y(:))) self.position(~self.outliers, :)];
            [b, ~, ~, ~, stats] = regress(y(:), X);
            p = stats(3);
            v0 = b(1);
            V = pinv(b(2:3));
            phi = atan2(V(2), V(1));
            speed = norm(V);
            
        end
        
        function fit = wave_fit(self)
            
            [V_, p, v0] = self.regress;
            fit.V = V_;
            fit.speed = norm(V_) * .4;  % in mm / s  
            fit.phi = atan2(V_(2), V_(1));
            fit.p = p;
            fit.v0 = v0;
        end
        
        function s_dur = seizure_dur(self)
            p = strsplit(self.name, '_');
            if strcmpi(p{1}, 'P4') 
                s_dur = 34;
            else
                load([p{1} filesep ...
                    self.name '_Neuroport_10_10.mat'], ...
                    'Duration');
                s_dur = Duration;
            end
        end
        
        function ctr = crossing_time_rel(self)
            ctr = diff(self.range) / self.seizure_dur;
        end
        
        function er = ending_rel(self)
            er = self.range(2) / self.seizure_dur;
        end
        
        function show = get.show(self)
            show = sum(~self.outliers) >= self.MinElectrodes;  
        end

        function ax = plot2D(self, type, ax)
            if nargin < 3, ax = gca; end
            if nargin < 2 || isempty(type), type = 'onsets'; end
            if isa(type, 'matlab.graphics.axis.Axes'), ax = type; type = 'onsets'; end
            type = validatestring(type, { ...
                'onsets', 'firingrate', 'duration' ...
                });
            
            switch type
                case 'onsets'
                    axis(ax, 'square')
                    set(ax, 'nextplot', 'replacechildren', 'units', 'points');
                    gs = self.GridSize;
                    pt_width = min(ax.Position([3 4])) / max(gs);
                    
                    
                    fit = self.wave_fit;
                    ons = self.onsets(~self.outliers);
                    pos = self.position(~self.outliers, :);
                    sz = self.fr_at_peak(~self.outliers);
                    szR = rescale(sz, (pt_width / 4).^2, (1.1*pt_width).^2);
                    scatter(ax, pos(:, 1), pos(:, 2), szR, ons, 'filled', 's');
                    
                    xlim([0 gs(1)+1]);
                    ylim([0 gs(2)+1]);
                    colorbar(ax);
                    title(ax, {strrep(self.name, '_', ' '); 'IW onset time (s)'});
                    
                    
                    % Add direction arrow
                    V_ = .8 * gs(1) * fit.V./vecnorm(fit.V);
                    x = V_/2;
                    hold on
                    quiver(ax, (gs(1) + 1)/2 - x(1), (gs(2) + 1)/2 - x(2), ...
                        V_(1), V_(2), 0, ...
                        'color', [0 0 0], 'linewidth', 2, 'maxheadsize', .5)
                    hold off
                    xstring = strsplit(sprintf('Speed (mm/s): %0.3f, phi: %0.1f, p=%0.5g', ...
                        fit.speed, atan2(V_(2), V_(1)) / pi * 180, fit.p), ',')';
                    xlabel(ax, xstring); 
                    
                case 'firingrate'
                    % Peak power during IW
                    data = self.fr_at_peak;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Firing rate at peak (spikes/s)'});
                    
                    axis(ax, 'square')
                    
                case 'duration'
                    % Peak power during IW
                    
                    data = self.durs;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Duration of IW [s]'});
                    
                    axis(ax, 'square')

            end
        end        
        
        function plot(self, style, ax, show_onsets)
            % self.plot(style='raster');  
            % Valid styles are {'raster', 'lowpass'}. 
            %   lowpass: full (normalized) low pass trace on each channel.
            %   raster: shows suprathreshold (self.MinPeakHeight) intervals
            %       on each channel with red dot indicating wave onset time
            
            if nargin < 2 || isempty(style), style = 'raster'; end  
            if nargin < 3 || isempty(ax), ax = gca; end
            if nargin < 4, show_onsets = false; end
            
            style = validatestring(style, ...
                {'raster', 'lowpass', 'raster_fr', 'fr', 'fr_pow'});
            switch style
                case 'fr_pow'
                    fr_mn = nanmean(self.fr, 2);
                    plo_mn = nanmean(self.p_lo, 2);
                    yyaxis(ax, 'left')
                    plot(ax, self.time, fr_mn, ...
                        'linewidth', .5, 'color', .8 * [1 1 1])
                    ylim(ax, [0 1.1*max(fr_mn)]);
                    ylabel(ax, 'Firing rate (spikes/s)');
                    yyaxis(ax, 'right')
                    plot(ax, self.t, plo_mn, ...
                        'linewidth', 2, 'color', .15 * [1 1 1])
                    ylim(ax, [0 1.1*max(plo_mn)]);
                    ylabel(ax, {'Power [0-2] Hz'; '(normalized)'})
                    
                    for ii = 1:self.num_waves
                        self.wave = ii;
                        if sum(isfinite(self.onsets)) < 10, continue; end
                        t0 = nanmedian(self.onsets);
%                         wave_range = ...
%                             [self.onsets(:) - self.durs(:)/2, ...
%                             self.onsets(:) + self.durs(:)/2];
%                         x_range = [min(wave_range(~self.outliers, 1)), ...
%                             max(wave_range(~self.outliers, 2))];
                        x_range = self.range;
                        y_range = ylim(ax);
                        xline(ax, t0);
                        
                        [v1, v2] = ndgrid(x_range, y_range);
                        pp = patch(ax, 'vertices', [v1(:) v2(:)], ...
                            'faces', [1 2 4 3], ...
                            'facecolor', .5*[1 1 1], ...
                            'facealpha', .5, ...
                            'linestyle', 'none'); %#ok<NASGU>
                    end

                    axis(ax, 'tight')
                    title(ax, {strrep(self.name, '_', ' ')});
                case 'fr'
                    
                    fr_mn = mean(self.fr, 2);
                    fr_mn1 = mean(self.fr_smooth, 2);
                    ff = self.vmr;
                    ff = ff.trace;

                    plot(ax, self.time, rescale(fr_mn, 0, max(fr_mn1(:))), ...
                        'linewidth', .5, 'color', .8 * [1 1 1], 'tag', 'fr')
                    ylabel(ax, {'Firing rate'; '(spikes/s)'});
                    hold(ax, 'on')
                    plot(ax, self.time, fr_mn1, ...
                        'linewidth', 2, 'color', .15 * [1 1 1], 'tag', 'frS');
                    
                    plot(ax, self.time, ff, ...
                        'linewidth', 2, 'color', .4 * [1 1 1], 'tag', 'VMR');
                    hold(ax, 'off')
                    
                    ylabel(ax, {'Firing rate'; '(spikes/s)'})
                    
                    for ii = 1:self.num_waves
                        self.wave = ii;
                        if ~self.show, continue; end
                        t0 = self.center;
                        x_range = self.range;
                        y_range = ylim(ax);
                        xline(ax, t0);
                        
                        [v1, v2] = ndgrid(x_range, y_range);
                        pp = patch(ax, 'vertices', [v1(:) v2(:)], ...
                            'faces', [1 2 4 3], ...
                            'facecolor', .5*[1 1 1], ...
                            'facealpha', .5, ...
                            'linestyle', 'none'); %#ok<NASGU>
                    end

                    axis(ax, 'tight')
                    dat = [fr_mn1, ff];
                    dat = dat(self.time > 0 & self.time < self.seizure_dur, :);
                    ylim(ax, [0 1.1*max(dat(:))])
                    title(ax, {strrep(self.name, '_', ' ')});
                case 'raster_fr'
                    frS = self.fr_smooth;
                    frS(self.time < 0 | self.time > self.seizure_dur, :) = nan;
                    frN = normalize(frS);
%                     [n_t, n_ch] = find(frN > self.MinPeakHeight & frS > self.MinPeakFr);
%                     plot(ax, self.time(n_t), n_ch, '.');
%                     hold on
%                     plot(self.onsets, 1:self.nch, 'r.');
%                     hold off
                    mat = frN > self.MinPeakHeight & frS > self.MinPeakFr;
                    [ons, so] = sort(self.onsets);
                    mat = mat(:, so);
%                     mat = frS > self.MinPeakFr;
%                     mat = false(size(frS));
%                     tt = self.time;
%                     for ii = 1:size(mat, 2)
%                         bounds = self.onsets(ii) + [-.5 .5] * self.durs(ii);
%                         mask = tt >= bounds(1) & tt <= bounds(2);
%                         mat(mask, ii) = true;
%                     end
                    imagesc(ax, self.time, 1:numel(self.onsets), mat');
                    if show_onsets
                        hold on
                        plot(ons,  1:self.nch, 'r.');
                        hold off;
                    end
                    colormap(ax, 1 - gray(2)); axis(ax, 'xy')
                    axis(ax, 'tight');
                    xlim(ax, [self.time(1) self.time(end)])
                    title(ax, {strrep(self.name, '_', ' '); ...
                        sprintf('Mean duration = %0.3f s', nanmean(self.durs))});
                    ylabel(ax, 'Channel')
                case 'raster'
%                     pwr = self.p_lo;
%                     inds = interp1(self.time, 1:numel(self.time), self.t, 'nearest');
%                     frS = self.fr_smooth(inds, :);
%                     [n_t, n_ch] = find(pwr > self.MinPeakHeight & frS > self.MinPeakFr);
%                     plot(ax, self.t(n_t), n_ch, '.');
%                     mat = pwr > self.MinPeakHeight & frS > self.MinPeakFr;
                    frS = self.fr_smooth;
                    mat = frS > self.MinPeakFr & normalize(frS) > self.MinPeakHeight;
                    imagesc(ax, self.time, 1:numel(self.onsets), mat');
                    if show_onsets
                        hold on
                        plot(self.onsets,  1:self.nch, 'r.');
                        hold off;
                    end
                    colormap(ax, 1-gray(2)); axis(ax, 'xy')
%                     hold on
%                     plot(ax, self.onsets, 1:self.nch, 'r.');
%                     hold off
                    axis(ax, 'tight');
                    xlim(ax, [self.time(1) self.time(end)])
                    title(ax, {strrep(self.name, '_', ' '); ...
                        sprintf('Mean duration = %0.3f s', nanmean(self.durs))});
                    ylabel(ax, 'Channel')
                case 'lowpass'
                    plot(ax, self.time, normalize(self.fr_smooth)/4 + (1:self.nch), 'k') 
                    axis(ax, 'tight')
                    title(ax, strrep(self.name, '_', ' '));
                    ylabel(ax, 'Channel')
            end
            xlabel(ax, 'Time (s)');
        end
        

    end
    
end



