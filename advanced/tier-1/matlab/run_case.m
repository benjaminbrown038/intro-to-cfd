function run_case()
% Solve u_t = alpha*u_xx, zero endpoint values, sine initial condition.
% Open this directory in MATLAB and run: run_case
base = fileparts(fileparts(mfilename('fullpath')));
c = jsondecode(fileread(fullfile(base, 'cases', 'sine_decay.json')));
out = fullfile(base, 'results', 'matlab');
if ~exist(out, 'dir'), mkdir(out); end
assert(c.diffusivity > 0 && c.length > 0 && c.final_time > 0);
assert(c.diffusion_number > 0 && c.diffusion_number <= 0.5);
counts = c.node_counts;
metrics = zeros(numel(counts), 9);
f1 = figure('Visible', 'off'); hold on;
for k = 1:numel(counts)
    n = counts(k);
    assert(n >= 3 && n == floor(n));
    x = linspace(0, c.length, n);
    dx = c.length/(n-1);
    steps = max(1, ceil(c.final_time/(c.diffusion_number*dx^2/c.diffusivity)));
    dt = c.final_time/steps;
    r = c.diffusivity*dt/dx^2;
    u = sin(pi*x/c.length); u([1,end]) = 0;
    timer = tic;
    for j = 1:steps
        u(2:end-1) = u(2:end-1) + r*(u(3:end)-2*u(2:end-1)+u(1:end-2));
    end
    elapsed = toc(timer);
    exact = sin(pi*x/c.length)*exp(-c.diffusivity*(pi/c.length)^2*c.final_time);
    exact([1,end]) = 0;
    err = u-exact;
    l2 = sqrt(dx*sum(err.^2));
    order = NaN;
    if k > 1
        order = log(metrics(k-1,6)/l2)/log(metrics(k-1,2)/dx);
    end
    metrics(k,:) = [n, dx, dt, steps, r, l2, max(abs(err)), elapsed, order];
    data = table(x,u,exact,err,'VariableNames',{'x','numerical','analytical','error'});
    writetable(data, fullfile(out,sprintf('solution_%d.csv',n)));
    plot(x,u,'DisplayName',sprintf('%d nodes',n));
end
plot(x,exact,'k--','DisplayName','Analytical');
xlabel('x (m)'); ylabel('Normalized scalar'); title('1D diffusion: final profile');
legend('show'); print(f1,fullfile(out,'solution.png'),'-dpng','-r160'); close(f1);
t = array2table(metrics,'VariableNames',{'nodes','dx','dt','steps','diffusion_number','l2_error','max_error','integration_seconds','observed_order'});
writetable(t,fullfile(out,'convergence.csv'));
f2 = figure('Visible','off'); loglog(metrics(:,2),metrics(:,6),'o-'); grid on;
xlabel('Grid spacing (m)'); ylabel('Discrete L2 error');
title('Refinement with dt proportional to dx squared');
print(f2,fullfile(out,'convergence.png'),'-dpng','-r160'); close(f2);
matlab_version = version;
save(fullfile(out,'run.mat'),'c','metrics','matlab_version');
disp(t);
end
