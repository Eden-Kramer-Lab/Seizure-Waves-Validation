function [data, P] = stable_intervals(M, varargin)
    % Computes the TW direction in sliding P.Win=S.Win windows (every S.Step s)
    % with low variance.
    %
    %
    % Returns:
    %   shifts: the change in direction between phases [rad]
    %   duration_of_phase: the duration of each phase [s]
    %   
    % 
    
    %% Parse inputs
    P = struct(...
        'Win', 5, ...  % Determine stability on 5 s windows
        'Step', .1, ...   % ... every 100 ms.
        'GapTolerance', 2, ...  % Allow gaps of 2 s without splitting the interval.
        'SplitDphiThresh', deg2rad(30), ...  % Split intervals if the direction trend changes by more than 30°
        'MinDuration', 2, ...  % Ignore/omit intervals with duration shorter than 2 s.
        'MinSamples', 10, ...  % Require at least 10 direction estimates in WIN for a time to be classed as stable
        'StableThresh', [], ...  % Require DI be at least .5◊.97 for a time to be classed as stable
        'AngStat', 'circ_modeS');  % Define the direction trend as the distribution mode and smooth direction trends
    
    for ii = 1:2:numel(varargin)
        key = validatestring(varargin{ii}, fieldnames(P));
        val = varargin{ii + 1};
        P.(key) = val;
    end
    
    if isempty(P.StableThresh)
        if class(M) == "MaxDescent"
            P.StableThresh = 0.5;
        elseif class(M) == "Delays"
            P.StableThresh = .97;
        end
    end
    
    % Function defining where to split intervals
    issplit_ = @(data) ...
        data.dt > P.GapTolerance ...
        | abs(data.dphi) > P.SplitDphiThresh ...
        | (abs(data.dphi./data.dt) > 3 * P.SplitDphiThresh ) ...  % & abs(data.dphi) > deg2rad(5)
        ;

    
    %% Load TW direction data and compute sliding window stats
    WNG = warning;

    % Load the data
    data = table(M.Time, M.Direction, 'VariableNames', {'time', 'dir'});
    
    % ... filter: only use samples where directions are estimated
    data = data(isfinite(data.dir), :);
    
    % Compute the range of times with direction estimates
    data.seizure_dur(:) = range(data.time);
    
    % Initialize additional variables
    for ff = ["circ_r", "circ_confmean", "circ_mean", "circ_mode", ...
            "dir_trend", "test_times", "dt", "dphi", "phase_num"]
        data.(ff)(:) = nan;
    end
    data.phase_num(:) = 0;  % nan's don't register as unique so use 0
    

    %% Get stable times
    di = circ_mov_stat('r', data.dir, data.time, P.Win, ...
                    'samplepoints', data.time);
    N = arrayfun(@(t) sum(abs(t - data.time) < P.Win/2), data.time);
    isstable = di > P.StableThresh & N >= P.MinSamples;

    
    if ~any(isstable)  % return early if there are no stable windows
        data = data(1, :);
        
        % Adjust some fields 
        [data.time, data.test_times] = deal(0);
        
        warning(WNG);
    	return;
    end

    %% Remove non-stable times
    data = data(isstable, :);

    %% Get directions at stable times
    data.test_times = round(data.time/P.Step) * P.Step;
            
    if class(M) == "MaxDescent"
        % Local functions to help compute sliding stats
        mask_ = @(tt) abs(data.time - tt) <= P.Win/2;
        sliding_stat_ = @(stat_fun, args) ...
            arrayfun(@(tt) stat_fun(args(mask_(tt), :)), data.test_times);

        % Compute statistics about the observations (TW directions) in each window
        warning('off', 'circ_confmean:requirementsNotMet');

        data.circ_r = sliding_stat_(@circ_r, data.dir);
        data.circ_confmean = sliding_stat_(@circ_confmean, data.dir);
        data.circ_mean = sliding_stat_(@circ_mean, data.dir);
        data.circ_mode = sliding_stat_(@circ_mode, data.dir);
        
        if contains(P.AngStat, "mean")
            data.dir_trend = data.circ_mean;
        elseif contains(P.AngStat, "mode")
            data.dir_trend = data.circ_mode;
        end
    
    elseif class(M) == "Delays" && M.HalfWin >= P.Win/2
        data.dir_trend = data.dir;
    end
    

    %% Remove any duplicate test times
    [~, uinds] = unique(data.test_times, 'stable');
    data = data(uinds, :);
    
    %% Smooth stats using a 1 s moving median filter
    if ismember(P.AngStat, ["circ_modeS", "circ_meanS"])
        data.dir_trend = fix_angle( ...
            movmedian(unwrap(data.dir_trend), P.Win, ...
            'SamplePoints', data.test_times));
    end

    %% Split data into intervals of stability
    data.dt = [P.Step; diff(data.time)];
    data.dphi = [0; diff_phase(data.dir_trend)];
    data.phase_num = cumsum(issplit_(data)) + 1;
    
    %% Remove phases that are too short
    dur = sortrows( ...
        groupsummary(data(:, ["test_times", "phase_num"]), "phase_num", 'range'), ...
        ['range_test_times', "GroupCount"], 'desc');
    keepers = dur.phase_num(dur.range_test_times >= P.MinDuration & dur.GroupCount >= 5);
    data = data(ismember(data.phase_num, keepers), :);
    
    
    %% Re-apply is_split_ 
    % This should appropriately merge any gaps that may have occurred if a
    % short disruption was removed.
    data.phase_num = cumsum(issplit_(data)) + 1;
    
    %% Again, remove phases that are too short and fix phase numbers
    dur = sortrows( ...
        groupsummary(data(:, ["test_times", "phase_num"]), "phase_num", 'range'), ...
        ['range_test_times', "GroupCount"], 'desc');
    keepers = dur.phase_num(dur.range_test_times >= P.MinDuration & dur.GroupCount >= 5);
    data = data(ismember(data.phase_num, keepers), :);
    data.phase_num = findgroups(data.phase_num);

    
    %% Return
    warning(WNG);
    
end
    