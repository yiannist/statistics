-module(quantile).

-compile(export_all).

%% Functions for approximating quantiles, i.e. points taken at regular
%% intervals from the cumulative distribution function of a random
%% variable.
%%
%% The number of quantiles is described below by the variable 'q', so
%% with q=4, a 4-quantile (also known as a "quartile") has 4
%% intervals, and contains 5 points. The parameter 'k' describes the
%% desired point, where 0 <= k <= q.


%% Estimate the /k/th /q/-quantile of a sample, using the weighted
%% average method.
-spec weightedAvg(integer(), integer(), [float()]) -> float().
weightedAvg(K, Q, X) ->
    N   = length(X),
    case N == 1 of
        true  -> hd(X);
        false ->
            %% Assertions:
            true = Q >= 2,
            true = K >= 0,
            true = K <  Q,
            IDX = (N - 1) * K / Q,
            J   = sample_powers:floor(IDX),
            SX  = lists:sort(X), %STUB: should be partial_sort!
            XJ  = lists:nth(J + 1, SX),
            XJ1 = lists:nth(J + 2, SX),
            G   = IDX - J,
            XJ + G * (XJ1 - XJ)
    end.

%% Estimate the k-th q-quantile of a sample 'x', using the continuous
%% sample method with the given parameters. This is the method used by
%% most statistical software, such as R, Mathematica, SPSS, and S.
-spec continuous_by({float(), float()}, integer(), integer(), [float()]) -> float().
continuous_by({A, B}, K, Q, X) ->
    true = Q >= 2,
    true = K >= 0,
    true = K =< Q,
    N = length(X),
    Eps = sample_histogram:m_epsilon() * 4,
    P = K / Q,
    T = A + P * (N + 1 - A - B),
    J = sample_powers:floor(T + Eps),
    R = T - J,
    H = case abs(R) < Eps of
            true  -> 0;
            false -> R
        end,
    Bracket = fun (M) -> min(max(M, 0), N - 1) end,
    SX = lists:sort(X), %STUB: should be partial_sort!
    Item = fun (I) -> lists:nth(Bracket(I) + 1, SX) end,
    (1 - H) * Item(J - 1) + H * Item(J).

%% Estimate the range between q-quantiles 1 and q-1 of a sample x,
%% using the continuous sample method with the given parameters.
%%
%% For instance, the interquartile range (IQR) can be estimated as
%% follows:
%%
%% > midspread medianUnbiased 4 (U.fromList [1,1,2,2,3])
%% > ==> 1.333333
-spec midspread({float(), float()}, integer(), [float()]) -> float().
midspread({A, B}, K, X) ->
    true = K > 0,
    N = length(X),
    Eps = sample_histogram:m_epsilon() * 4,
    T = fun (I) -> A + I * (N + 1 - A - B) end,
    J = fun (I) -> sample_powers:floor(T(I) + Eps) end,
    R = fun (I) -> T(I) - J(I) end,
    H = fun (I) ->
                Ri = R(I),
                case abs(Ri) < Eps of
                    true  -> 0;
                    false -> Ri
                end
        end,
    Bracket = fun (M) -> min(max(M, 0), N - 1) end,
    SX = lists:sort(X), %STUB: should be partial_sort!
    Item = fun (I) -> lists:nth(Bracket(I) + 1, SX) end,
    Frac = 1 / K,
    Quantile = fun (I) -> (1 - H(I)) * Item(J(I - 1)) + H(I) * Item(J(I)) end,
    Quantile(1 - Frac) - Quantile(Frac).
