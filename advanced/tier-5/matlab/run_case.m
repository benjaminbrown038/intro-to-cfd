function run_case()
% Tier 5: see ../README.md for scope and run instructions.
base=fileparts(fileparts(mfilename('fullpath')));
root=fileparts(fileparts(base));
addpath(fullfile(root,'common','matlab'));
cfd_run_cases(5,base);
end
