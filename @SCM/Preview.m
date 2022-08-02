function h = Preview(sim, show_layout)
    % Generate panels showing the state of the sim (Qe and K) at each save
    % point.
    
    % Parse inputs
    if nargin < 2 || isempty(show_layout), show_layout = false; end
    
    % Load Qe and K movies if not already done
    if isempty(sim.Qe_movie)
        load(sprintf('%s_%d_info.mat', sim.basename, sim.sim_num), 'Qe_movie', 'K_movie');
        aspect_ratio = sim.grid_size(2) / sim.grid_size(1);
        for ii = 1:length(Qe_movie)
            Qe_movie(ii).cdata = ...
                imresize(Qe_movie(ii).cdata, [500 500 * aspect_ratio]);
            K_movie(ii).cdata = ...
                imresize(K_movie(ii).cdata, [500 500 * aspect_ratio]);
        end
        sim.Qe_movie = Qe_movie;
        sim.K_movie = K_movie;
    end

    % Make the figures
    fig_names = ["Qe_movie" "K_movie"];
    h = arrayfun(@(nn) figure('units', 'inches', ...
        'position', [0    0.8472   10.6250   10.2083], ...
        'name', nn), fig_names); 
    dims = size(sim.Qe_movie(1).cdata);
    
    % Plot the data
    for ff = fig_names
        T = tiledlayout(h(ff == fig_names), 10, 10, 'tilespacing', 'compact');
        for ii = 1:min(numel(sim.(ff)), 100)
            ax = nexttile(T, ii); 
            imagesc(ax, sim.(ff)(ii).cdata); 
            colormap(ax, sim.(ff)(ii).colormap);
            xticks(ax, []); yticks(ax, []);
            axis(ax, 'square');
            if show_layout
                show_layout_(sim, ax);
            end
            text(dims(1), dims(2), sprintf('t = %d s ', ii - sim.padding(1)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end

    end
    

end

%% Local fun

function show_layout_(sim, ax)
    
    
    [mea_x, mea_y] = ind2sub(sim.grid_size, sim.NPinds);
    k_mea = boundary(mea_x, mea_y);
    
    Nsrcs = size(sim.source, 3);
    [src_x, src_y] = arrayfun(@(ii) find(sim.source(:, :, ii) > 0), ...
        1:Nsrcs, 'uni', 0);
    src_x = cellfun(@mean, src_x);
    src_y = cellfun(@mean, src_y);
    dx = sim.dx;
    
    % Show the IZ
    [iz_x, iz_y] = find(sim.IZ ~= 0);
    k_iz = boundary(iz_x, iz_y);
    % ... smooth the IZ boundary
    ctr = [mean(iz_x) mean(iz_y)]/dx;
    rad = mean(vecnorm(([iz_x(k_iz)/dx, iz_y(k_iz)/dx] - ctr), 2, 2));
    theta = linspace(-pi, pi, 200);
    iz_x = rad*(cos(theta)) + ctr(1);
    iz_y = rad*(sin(theta)) + ctr(2);
    
    hold(ax, 'on')
    pp = fill(ax, mea_y(k_mea)/dx, mea_x(k_mea)/dx, 'r');
    set(pp, 'facealpha', .5, 'edgecolor', 'none');
    plot(ax, src_y/dx, src_x/dx, 'r.', 'markersize', 10);
    plot(ax, iz_y, iz_x, 'r:');
    hold(ax, 'off')
    
    
end

