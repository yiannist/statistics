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
-spec continuous_by({float(), float()}, integer(), integer(), [float()])
                   -> float().
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
    Frac = 1 / K,
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
    Quantile = fun (I) -> (1 - H(I)) * Item(J(I) - 1) + H(I) * Item(J(I)) end,
    Quantile(1 - Frac) - Quantile(Frac).

%% California Department of Public Works definition, a=0, b=1.
%% Gives a linear interpolation of the empirical CDF.  This
%% corresponds to method 4 in R and Mathematica.
-spec cadpw() -> {float(), float()}.
cadpw() -> {0.0, 1.0}.

%% Hazen's definition, a=0.5, b=0.5. This is claimed to be
%% popular among hydrologists. This corresponds to method 5 in R and
%% Mathematica.
-spec hazen() -> {float(), float()}.
hazen() -> {0.5, 0.5}.

%% Definition used by the SPSS statistics application, with a=0,
%% b=0 (also known as Weibull's definition). This corresponds to
%% method 6 in R and Mathematica.
-spec spss() -> {float(), float()}.
spss() -> {0.0, 0.0}.

%% Definition used by the S statistics application, with a=1,
%% b=1. The interpolation points divide the sample range into n-1
%% intervals. This corresponds to method 7 in R and Mathematica.
-spec s() -> {float(), float()}.
s() -> {1.0, 1.0}.

%% Median unbiased definition, a=1/3, b=1/3. The resulting
%% quantile estimates are approximately median unbiased regardless of
%% the distribution of x. This corresponds to method 8 in R and
%% Mathematica.
-spec medianUnbiased() -> {float(), float()}.
medianUnbiased() -> {1/3, 1/3}.


%% Normal unbiased definition, a=3/8, b=3/8. An approximately
%% unbiased estimate if the empirical distribution approximates the
%% normal distribution.  This corresponds to method 9 in R and
%% Mathematica.
-spec normalUnbiased() -> {float(), float()}.
normalUnbiased() -> {3/8, 3/8}.


%% References:
%%
%% * Weisstein, E.W. Quantile. /MathWorld/.
%%   <http://mathworld.wolfram.com/Quantile.html>
%%
%% * Hyndman, R.J.; Fan, Y. (1996) Sample quantiles in statistical
%%   packages. /American Statistician/
%%   50(4):361&#8211;365. <http://www.jstor.org/stable/2684934>
