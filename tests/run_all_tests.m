function fails = run_all_tests
%RUN_ALL_TESTS Run all three GraphWaveletOptimization verification harnesses.
%   Returns the total number of failed test groups and raises an error
%   if any check fails, so that CI exits nonzero on regression.
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end
fails = 0;
fprintf('==== Module 1: run_tests ====\n');   fails = fails + run_tests();
fprintf('\n==== Module 2: run_tests2 ====\n'); fails = fails + run_tests2();
fprintf('\n==== Module 3: run_tests3 ====\n'); fails = fails + run_tests3();
fprintf('\n==== TOTAL: %d failing test group(s) ====\n', fails);
if fails > 0
    error('GraphWaveletOptimization:testFailure', '%d test group(s) failed.', fails);
end
end
