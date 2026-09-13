function run_all()
% Run all supplied MATLAB examples. MATLAB execution is not verified here.
root=fileparts(mfilename('fullpath')); original=pwd;
cleanup=onCleanup(@()cd(original));
for tier=1:5
    cd(fullfile(root,'advanced',sprintf('tier-%d',tier),'matlab'));
    clear run_case;
    fprintf('Running Tier %d\n',tier);
    run_case;
end
end
