-module(sample_powers).

-export([
        %% Constructor
          powers/2
        %% Descriptive functions
        , order/1
        , count/1
        , sum/1
        %% Statistics of location
        , mean/1
        %% Statistics of dispersion
        , variance/1
        , std_dev/1
        , variance_unbiased/1
        %% Functions over central moments
        , central_moment/2
        , skewness/1
        , kurtosis/1
        %% Utility
        , floor/1
        ]).

-export_type([
             %% Types
               powers/0
             ]).

%% Fast statistics over simple powers of a sample. These can all be
%% computed efficiently (in just a single (?) pass over a sample).
%%
%% The tradeoff is that some of these functions are less numerically
%% robust than their counterparts in the 'statistics_sample' module.
%% Where this is the case, the alternatives are noted.

-type powers() :: {'powers', [float()]}.

%% Functions computed over a sample's simple powers require at least a
%% certain number (or order) of powers to be collected.
%%
%% * To compute the k-th 'centralMoment', at least k simple powers
%%   must be collected.
%%
%% * For the 'variance', at least 2 simple powers are needed.
%%
%% * For 'skewness', we need at least 3 simple powers.
%%
%% * For 'kurtosis', at least 4 simple powers are required.
-spec powers(integer(), [float()]) -> powers().
powers(K, _) when K < 2 -> exit({?MODULE, powers, "too few powers"});
powers(K, XS) ->
    L = K + 1,
    T = ets:new(tab, [named_table]),
    true = ets:insert(T, lists:zip(lists:seq(0, L-1),
                                   lists:duplicate(L, 0.0))), % mutable.
    Go = fun (X, MS) ->
                 ok = powers__loop({MS, X, L}, 0, 1),
                 MS
         end,
    T = lists:foldl(Go, T, XS),
    {_, PS} = lists:unzip(lists:sort(fun ({A1, _}, {A2, _}) -> A1 < A2 end,
                                     ets:tab2list(T))), % freeze.
    true = ets:delete(T),
    {'powers', lists:sort(PS)}.

-spec powers__loop({ets:tab(), float(), integer()}, integer(), integer())
                  -> ok.
powers__loop({_, _, L}, I, _) when I == L -> ok; % MSTab is "global".
powers__loop({MSTab, X, L}, I, XK) ->
    [{I, M}] = ets:lookup(MSTab, I), % 'I' intentionally the same
    true = ets:insert(MSTab, [{I, M + XK}]),
    powers__loop({MSTab, X, L}, I + 1, XK * X).

%% The order (number) of simple powers collected from a 'sample'.
-spec order(powers()) -> integer().
order({'powers', XS}) ->
    length(XS) - 1.

%% Compute the k-th central moment of a sample. The central moment
%% is also known as the moment about the mean.
%%STUB: sample_powers:central_moment(40, sample_powers:powers(52, lists:seq(1, 42000))).
%% WRONG result.
-spec central_moment(integer(), powers()) -> float().
central_moment(0, _) -> 1;
central_moment(K, P={'powers', PA}) ->
    case ((K < 0) or (K > order(P))) of
        true  -> error({?MODULE, central_moment, "invalid argument"});
        false -> N = hd(PA),
                 M = mean(P),
                 Go = fun ({I, E}) ->
                              spec_functions:choose(K, I) * math:pow(-M, K-I) * E
                      end,
                 Indexed = function:indexed(lists:sublist(PA, K + 1)),
                 lists:sum(lists:map(Go, Indexed)) / N
    end.

%% Maximum likelihood estimate of a sample's variance.  Also known
%% as the population variance, where the denominator is n. This is
%% the second central moment of the sample.
%%
%% This is less numerically robust than the variance function in the
%% 'Statistics.Sample' module, but the number is essentially free to
%% compute if you have already collected a sample's simple powers.
%%
%% Requires 'Powers' with 'order' at least 2.
-spec variance(powers()) -> float().
variance(P) ->
    central_moment(2, P).

%% Standard deviation. This is simply the square root of the
%% maximum likelihood estimate of the variance.
-spec std_dev(powers()) -> float().
std_dev(P) ->
    math:sqrt(variance(P)).

%% Unbiased estimate of a sample's variance. Also known as the
%% sample variance, where the denominator is n-1.
%%
%% Requires 'Powers' with 'order' at least 2.
-spec variance_unbiased(powers()) -> float().
variance_unbiased(P={'powers', PA}) when hd(PA) > 1 ->
    N = hd(PA),
    variance(P) * N / (N-1);
variance_unbiased(_) -> 0.0.

%% Compute the skewness of a sample. This is a measure of the
%% asymmetry of its distribution.
%%
%% A sample with negative skew is said to be "left-skewed". Most of
%% its mass is on the right of the distribution, with the tail on the
%% left.
%%
%% 1> skewness(3, powers([1,100,101,102,103])).
%% -1.497681449918257
%%
%% A sample with positive skew is said to be "right-skewed".
%%
%% 2> skewness(3, powers([1,2,3,4,100])).
%% 1.4975367033335198
%%
%% A sample's skewness is not defined if its 'variance' is zero.
%%
%% Requires 'Powers' with 'order' at least 3.
-spec skewness(powers()) -> float().
skewness(P) ->
    central_moment(3, P) * math:pow(variance(P), -1.5).

%% Compute the excess kurtosis of a sample. This is a measure of
%% the "peakedness" of its distribution. A high kurtosis indicates
%% that the sample's variance is due more to infrequent severe
%% deviations than to frequent modest deviations.
%%
%% A sample's excess kurtosis is not defined if its 'variance' is
%% zero.
%%
%% Requires 'Powers' with 'order' at least 4.
-spec kurtosis(powers()) -> float().
kurtosis(P) ->
    V = variance(P),
    central_moment(4, P) / (V * V) - 3.

%% The number of elements in the original 'Sample'.  This is the
%% sample's zeroth simple power.
-spec count(powers()) -> integer().
count({'powers', PA}) ->
    floor(hd(PA)).

-spec floor(float()) -> integer().
floor(X) ->
    T = erlang:trunc(X),
    case (X - T) < 0 of
        true  -> T - 1;
        false -> T
    end.

%% The sum of elements in the original 'Sample'.  This is the
%% sample's first simple power.
-spec sum(powers()) -> float().
sum({'powers', PA}) ->
    lists:nth(2, PA).

%% The arithmetic mean of elements in the original 'Sample'.
%%
%% This is less numerically robust than the mean function in the
%% 'sample' module, but the number is essentially free to
%% compute if you have already collected a sample's simple powers.
-spec mean(powers()) -> float().
mean({'powers', PA}) when hd(PA) == 0 -> 0.0;
mean(P={'powers', PA}) ->
    N = hd(PA),
    sum(P) / N.


%% References:
%%
%% * Besset, D.H. (2000) Elements of statistics.
%%   /Object-oriented implementation of numerical methods/
%%   ch. 9, pp. 311&#8211;331.
%%   <http://www.elsevier.com/wps/product/cws_home/677916>
%%
%% * Anderson, G. (2009) Compute k-th central moments in one
%%   pass. /quantblog/. <http://quantblog.wordpress.com/2009/02/07/compute-kth-central-moments-in-one-pass/>
