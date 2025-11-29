function run_tests()
    addpath(genpath(pwd));
    testFolder = fullfile(pwd,'tests');
    if ~exist(testFolder,'dir')
        error('Tests folder not found: %s', testFolder);
    end

    % Run tests and display results
    results = runtests(testFolder);
    disp(results);

    % Save out a summary to tests/
    outFile = fullfile(pwd,'tests','test_results.mat');
    save(outFile,'results');
    fprintf('Saved test results to %s\n', outFile);

    % Print failures in readable form
    failed = results([results.Failed] == 1);
    if ~isempty(failed)
        fprintf('Some tests failed (%d). See details above and test_results.mat\n', numel(failed));
    else
        fprintf('All tests passed.\n');
    end
end
