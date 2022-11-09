function [fig] = plotIonoDelay(t, I, legend_str, fig )
t = (t - t(1)) / 60; % format time vector
idx = find(t>0.15);

if (strcmp(legend_str, 'carrier'))
    I(idx) = I(idx) - I(idx(1));
end

if ishandle(fig)
    fig = figure(fig);
    hold on;
else
    fig = figure();
end
plot(t(idx), I(idx), 'LineWidth', 2, 'DisplayName', legend_str)
xlabel('Time [min]')
ylabel('Ionospheric delay [meters]')
legend()

end

