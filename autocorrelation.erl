-module(autocorrelation).

-export([
          autocovariance/1
        , autocorrelation/1
        ]).

%% Functions for computing autocovariance and autocorrelation of a sample.


%% Compute the autocovariance of a sample, i.e. the covariance of the
%% sample against a shifted version of itself.
-spec autocovariance([float()]) -> [float()].
autocovariance(A) ->
    L = length(A),
    Mean = sample:mean(A),
    C = lists:map(fun (X) -> X - Mean end, A),
    F = fun (K) ->
                lists:sum(lists:zipwith(fun (X, Y) -> X * Y end,
                                        lists:sublist(C, L-K),
                                        lists:sublist(C, K+1, L))) / L
        end,
    lists:map(F, lists:seq(0, L-2)).

%% Compute the autocorrelation function of a sample, and the upper
%% and lower bounds of confidence intervals for each element.
%%
%% Note: The calculation of the 95% confidence interval assumes a
%% stationary Gaussian process.
-spec autocorrelation([float()]) -> {[float()], [float()], [float()]}.
autocorrelation(A) ->
    C = autocovariance(A),
    H = hd(C),
    R = lists:map(fun (X) -> X / H end, C),
    L = length(A),
    F0 = fun (V) -> 1.96 * math:sqrt((V * 2 + 1) / L) end,
    DLLSE = lists:map(F0, scanl1(fun (X, Y) -> X + Y end,
                                 lists:map(fun (X) -> X * X end, R))),
    CI = fun (F) -> [1 | tl(lists:map(fun (X) -> F(-1 / L, X) end, DLLSE))] end,
    {R, CI(fun (X, Y) -> X - Y end), CI(fun (X, Y) -> X + Y end)}.

%% Utility functions

%%-type scan_fun(A, B) :: fun ((A) -> fun ((B) -> A)).
%%-spec scanl(scan_fun(A, B), A, [B]) -> [A].
scanl(_, Q, [])     -> [Q];
scanl(F, Q, [X|XS]) -> [Q | scanl(F, F(Q, X), XS)].

%%-spec scanl1(scan_fun(A, A), [A]) -> [A].
scanl1(F, [X|XS]) -> scanl(F, X, XS);
scanl1(_, [])     -> [].
