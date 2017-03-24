function [C,lags] = xcorrPeriodic(A,B,maxlags)
    % xcorrPeriodic - Computes the cross correlation of two periodic signals
    %
    %   The built-in cross-correlation function is limited because in order to
    %   take the cross-correlation of periodic signals, you have to repmat each
    %   of them, and the resulting correlations are HIGHLY dependent upon how
    %   many repeats you choose to do. Also, you have to worry about edge bias.
    %
    %   This implementation shifts them in time and loops the signals around
    %   when necessary, thus producing a reliable and unbiased
    %   cross-correlation
    %
    %   The output is strictly the correlation coefficient because really all
    %   other outputs are basically garbage
    %
    %   A is the stationary signal, while the lags correspond to B
    %
    %   NOTE:   If you specify a maxlag > M/2, you will obtain two maximum
    %           correlations in your result. This is because the signals shift
    %           past eachother twice and so you get maximum correlations
    %           shifted by 2*pi.
    %
    % USAGE:
    %   [corr,lags] = xcorrPeriodic(A,B,maxlags)
    %
    % INPUTS:
    %   A:          [M x 1] Vector, Periodic reference signal
    %   B:          [M x 1] Vector, Periodic signal that is shited in time
    %   maxlags:    Integer, specifies the maximum delay (either positive or
    %               negative) that we will compute. Defaults to (M/2).
    %
    % OUTPUTS:
    %   corr:       [1 x maxlags*2+1] Array, Correlation values computed for
    %               each time shift
    %   lags:       [1 x maglags*2+1] Array, Phase lag corresponding to each
    %               of the CORR values.
    %
    % Last Modified: 03-14-2012
    % Modified By: Suever (suever@gmail.com)
    assert(numel(A)==numel(B),'xcorrPeriodic:UnEqualInputs',...
        'Inputs must be the same length since they have the same period');

    % De-mean and normalize each of the signals
    A = (A(:) - mean(A(:))) / std(A(:));
    B = (B(:) - mean(B(:))) / std(B(:));

    % By default use the Matlab maxlags
    if ~exist('maxlags','var')
        maxlags = ceil(numel(A)/2);
    end

    if numel(maxlags) == 2
        lags = maxlags(1):maxlags(2);
    else
        lags = -maxlags:maxlags;
    end

    C = lags;                     % initialize C

    % Pre-compute these guys for speed
    nPoints     = numel(A);
    originds    = 1:nPoints;

    for i = 1:numel(lags)
        % Shifted indices for B
        inds = mod(originds - lags(i) - 1,nPoints) + 1;
        C(i) = A' * B(inds);
    end

    C = C / (nPoints - 1);
end
