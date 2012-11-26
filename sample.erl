-module(sample).

-export([
        %% Descriptive functions
          range/1

        %% Statistics of location
        , mean/1
        , mean_weighted/1
        , harmonic_mean/1
        , geometric_mean/1

        %% Statistics of dispersion

        %% Functions over central moments
        , central_moment/2
        , central_moments/3
        , skewness/1
        , kurtosis/1

        %% Multi-pass functions (numerically robust)
        , variance/1
        , variance_unbiased/1
        , mean_variance/1
        , mean_variance_unb/1
        , std_dev/1
        , variance_weighted/1

        %% Single-pass functions (faster, less safe)
        , fast_variance/1
        , fast_variance_unbiased/1
        , fast_std_dev/1
        ]).

%% Range. The difference between the largest and smallest
%% elements of a sample.
-spec range([float()]) -> float().
range(XS) ->
    {Min, Max} = function:min_max(XS),
    Max - Min.

%% Arithmetic mean. This uses Welford's algorithm to provide
%% numerical stability, using a single pass over the sample data.
-spec mean([float()]) -> float().
mean(XS) ->
    Go = fun (X, {M, N}) ->
                 N2 = N + 1,
                 M2 = M + (X - M) / N2,
                 {M2, N2}
         end,
    {A, _} = lists:foldl(Go, {0, 0}, XS),
    A.

%% Arithmetic mean for weighted sample. It uses a single-pass
%% algorithm analogous to the one used by 'mean'.
-spec mean_weighted([{float(), float()}]) -> float().
mean_weighted(XS) ->
    Go = fun ({X, XW}, {M, W}) ->
                 W2 = W + XW,
                 M2 = case W2 == 0 of
                          true  -> 0;
                          false -> M + XW * (X - M) / W2
                      end,
                 {M2, W2}
         end,
    {A, _} = lists:foldl(Go, {0, 0}, XS),
    A.

%% Harmonic mean.  This algorithm performs a single pass over
%% the sample.
-spec harmonic_mean([float()]) -> float().
harmonic_mean(XS) ->
    Go = fun (N, {X, Y}) ->
                 {X + 1 / N, Y + 1}
         end,
    {A, B} = lists:foldl(Go, {0, 0}, XS),
    B / A.

%% Geometric mean of a sample containing no negative values.
-spec geometric_mean([float()]) -> float().
geometric_mean(XS) ->
    Go = fun (A, {P, N}) ->
                 {P * A, N + 1}
         end,
    {P, N} = lists:foldl(Go, {1, 0}, XS),
    math:pow(P, 1 / N).

%% Compute the k-th central moment of a sample.  The central moment
%% is also known as the moment about the mean.
%%
%% For samples containing many values very close to the mean, this
%% function is subject to inaccuracy due to catastrophic cancellation.
%%
%% Example from WolframAlpha:
%%   CentralMoment[{4, 36, 45, 50, 75}, 3] == -3722.4
-spec central_moment(integer(), [float()]) -> float().
central_moment(A, _) when A < 0 ->
    exit({?MODULE, central_moment, "negative input"});
central_moment(A, _) when A == 0 -> 1.0;
central_moment(A, _) when A == 1 -> 0.0;
central_moment(A, XS) ->
    Go = fun (X) ->
                 M = mean(XS),
                 math:pow(X - M, A)
         end,
    lists:sum(lists:map(Go, XS)) / length(XS).

%% Compute the k-th and j-th central moments of a sample.
%%
%% For samples containing many values very close to the mean, this
%% function is subject to inaccuracy due to catastrophic cancellation.
-spec central_moments(integer(), integer(), [float()]) -> {float(), float()}.
central_moments(A, B, XS) when (A < 2) or (B < 2) ->
    {central_moment(A, XS), central_moment(B, XS)};
central_moments(A, B, XS) ->
    Go = fun (X, {I, J}) ->
                 M = mean(XS),
                 D = X - M,
                 {I + math:pow(D, A), J + math:pow(D, B)}
         end,
    {I, J} = lists:foldl(Go, {0, 0}, XS),
    N = length(XS),
    {I / N, J / N}.

%% Compute the skewness of a sample. This is a measure of the
%% asymmetry of its distribution.
%%
%% A sample with negative skew is said to be "left-skewed". Most of
%% its mass is on the right of the distribution, with the tail on the
%% left.
%%
%% 1> skewness([1,100,101,102,103]).
%% -1.497681449918257607
%%
%% A sample with positive skew is said to be "right-skewed".
%%
%% 2> skewness([1,2,3,4,100]).
%% 1.4975367033335198
%%
%% A sample's skewness is not defined if its 'variance' is zero.
%%
%% For samples containing many values very close to the mean, this
%% function is subject to inaccuracy due to catastrophic cancellation.
-spec skewness([float()]) -> float().
skewness(XS) ->
    {C3, C2} = central_moments(3, 2, XS),
    C3 * math:pow(C2, -1.5).

%% Compute the excess kurtosis of a sample. This is a measure of
%% the "peakedness" of its distribution. A high kurtosis indicates
%% that more of the sample's variance is due to infrequent severe
%% deviations, rather than more frequent modest deviations.
%%
%% A sample's excess kurtosis is not defined if its 'variance' is
%% zero.
%%
%% For samples containing many values very close to the mean, this
%% function is subject to inaccuracy due to catastrophic cancellation.
-spec kurtosis([float()]) -> float().
kurtosis(XS) ->
    {C4, C2} = central_moments(4, 2, XS),
    C4 / (C2 * C2) -3.


%% The 'variance'and hence the 'standard deviation' of a sample
%% of fewer than two elements are both defined to be zero.
%%
%% These functions use the compensated summation algorithm of Chan et
%% al. for numerical robustness, but require two passes over the
%% sample data as a result.

-spec robust_sum_var(float(), [float()]) -> float().
robust_sum_var(M, Samp) ->
    lists:sum(lists:map(fun (I) ->
                                math:pow(I-M, 2)
                        end, Samp)).

%% Maximum likelihood estimate of a sample's variance. Also known
%% as the population variance, where the denominator is 'n'.
-spec variance([float()]) -> float().
variance(Samp) ->
    N = length(Samp),
    case N > 1 of
        true  -> robust_sum_var(mean(Samp), Samp) / N;
        false -> 0
    end.

%% Unbiased estimate of a sample's variance.  Also known as the
%% sample variance, where the denominator is 'n'-1.
-spec variance_unbiased([float()]) -> float().
variance_unbiased(Samp) ->
    N = length(Samp),
    case N > 1 of
        true  -> robust_sum_var(mean(Samp), Samp) / (N - 1);
        false -> 0
    end.

%% Calculate mean and maximum likelihood estimate of variance. This
%% function should be used if both mean and variance are required
%% since it will calculate mean only once.
-spec mean_variance([float()]) -> {float(), float()}.
mean_variance(Samp) ->
    N = length(Samp),
    M = mean(Samp),
    case N > 1 of
        true  -> {M, robust_sum_var(M, Samp)/ N};
        false -> {M, 0}
    end.

%% Calculate mean and unbiased estimate of variance. This
%% function should be used if both mean and variance are required
%% since it will calculate mean only once.
-spec mean_variance_unb([float()]) -> {float(), float()}.
mean_variance_unb(Samp) ->
    N = length(Samp),
    M = mean(Samp),
    case N > 1 of
        true  -> {M, robust_sum_var(M, Samp)/ (N - 1)};
        false -> {M, 0}
    end.

%% Standard deviation.  This is simply the square root of the
%% unbiased estimate of the variance.
-spec std_dev([float()]) -> float().
std_dev(Samp) ->
    math:sqrt(variance_unbiased(Samp)).

-spec robust_sum_var_weighted([{float(), float()}]) -> {float(), float()}.
robust_sum_var_weighted(Samp) ->
    M = mean_weighted(Samp),
    Go = fun ({X, XW}, {S, W}) ->
                 D = X - M,
                 {S + XW * D * D, W + XW}
         end,
    lists:foldl(Go, {0, 0}, Samp).

%% Weighted variance. This is biased estimation.
-spec variance_weighted([{float(), float()}]) -> float().
variance_weighted(Samp) when length(Samp) > 1 ->
    {S, W} = robust_sum_var_weighted(Samp),
    S / W;
variance_weighted(_) -> 0.


%% The functions prefixed with the name 'fast' below perform a single
%% pass over the sample data using Knuth's algorithm. They usually
%% work well, but see below for caveats.
%%
%% Note: in cases where most sample data is close to the sample's
%% mean, Knuth's algorithm gives inaccurate results due to
%% catastrophic cancellation.

-spec fast_var([float()]) -> {float(), float(), float()}.
fast_var(Samp) ->
    Go = fun (X, {N, M, S}) ->
                 D  = X - M,
                 N2 = N + 1,
                 M2 = M + D / N2,
                 S2 = S + D * (X - M2),
                 {N2, M2, S2}
         end,
    lists:foldl(Go, {0, 0, 0}, Samp).

%% Maximum likelihood estimate of a sample's variance.
-spec fast_variance([float()]) -> float().
fast_variance(Samp) ->
    {N, _M, S} = fast_var(Samp),
    case N > 1 of
        true  -> S / N;
        false -> 0
    end.

%% Unbiased estimate of a sample's variance.
-spec fast_variance_unbiased([float()]) -> float().
fast_variance_unbiased(Samp) ->
    {N, _M, S} = fast_var(Samp),
    case N > 1 of
        true  -> S / (N - 1);
        false -> 0
    end.

%% Standard deviation. This is simply the square root of the
%% maximum likelihood estimate of the variance.
-spec fast_std_dev([float()]) -> float().
fast_std_dev(Samp) ->
    math:sqrt(fast_variance(Samp)).


%% References:
%%
%% * Chan, T. F.; Golub, G.H.; LeVeque, R.J. (1979) Updating formulae
%%   and a pairwise algorithm for computing sample
%%   variances. Technical Report STAN-CS-79-773, Department of
%%   Computer Science, Stanford
%%   University. <ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/79/773/CS-TR-79-773.pdf>
%%
%% * Knuth, D.E. (1998) The art of computer programming, volume 2:
%%   seminumerical algorithms, 3rd ed., p. 232.
%%
%% * Welford, B.P. (1962) Note on a method for calculating corrected
%%   sums of squares and products. /Technometrics/
%%   4(3):419&#8211;420. <http://www.jstor.org/stable/1266577>
%%
%% * West, D.H.D. (1979) Updating mean and variance estimates: an
%%   improved method. /Communications of the ACM/
%%   22(9):532&#8211;535. <http://doi.acm.org/10.1145/359146.359153>
