-module(sample_histogram).

-export([
          histogram/2
        %% Building blocks
        , histogram_/4
        , range/2
        ]).

%% Compute a histogram over a data set.
%%
%% The result consists of a pair of lists:
%%
%% * The lower bound of each interval.
%% * The number of samples within the interval.
%%
%% Interval (bin) sizes are uniform, and the upper and lower bounds
%% are chosen automatically using the 'range' function.  To specify
%% these parameters directly, use the 'histogram_' function.
-spec histogram(integer(), [float()]) -> {[float()], [float()]}.
histogram(NumBins, XS) ->
    {Lo, Hi} = range(NumBins, XS),
    D = (Hi - Lo) / NumBins,
    Step = fun (I) ->
                   Lo + D * I
           end,
    {[Step(I) || I <- lists:seq(0, NumBins - 1)]
      , histogram_(NumBins, Lo, Hi, XS)}.

%% Compute a histogram over a data set.
%%
%% Interval (bin) sizes are uniform, based on the supplied upper
%% and lower bounds.
-spec histogram_(integer(), float(), float(), [float()]) -> [float()].
histogram_(NumBins, Lo, Hi, XS0) ->
    Bin = fun (XS, Bins) ->
                  D = ((Hi - Lo) * (1 + m_epsilon())) / NumBins,
                  histogram__go(XS, Lo, D, length(XS), 0, Bins)
          end,
    Bins = ets:new(bins, [named_table]),
    true = ets:insert(Bins, [{I, 0.0} || I <- lists:seq(0, NumBins - 1)]),
    Bins = Bin(XS0, Bins),        % Fills 'Bins' ETS table through side-effects!
    {_, Elems} = lists:unzip(ets:tab2list(Bins)),
    true = ets:delete(Bins),
    Elems.

-spec histogram__go([float()], float(), float(), integer(), integer()
                    , ets:tab()) -> ets:tab().
histogram__go(_,  _,  _, Len, I, Bins) when I >= Len -> Bins;
histogram__go(XS, Lo, D, Len, I, Bins) ->
    X = lists:nth(I+1, XS),
    B = trunc((X - Lo) / D),
    [{_, Belem}] = ets:lookup(Bins, B), % Matches the 'elem', ignores the index.
    true = ets:insert(Bins, [{B, Belem + 1}]),
    histogram__go(XS, Lo, D, Len, I + 1, Bins).

%% Compute decent defaults for the lower and upper bounds of a histogram,
%% based on the desired number of bins and the range of the sample data.
%%
%% The upper and lower bounds used are @(lo-d, hi+d)@, where
%%
%% D = (maximum sample - minimum sample) / ((bins - 1) * 2)
-spec range(integer(), [float()]) -> {float(), float()}.
range(NumBins, _) when NumBins < 1 ->
    exit({?MODULE, range, "invalid bin count"});
range(_, []) ->
    exit({?MODULE, range, "empty sample"});
range(NumBins, XS) ->
    {Lo, Hi} = function:min_max(XS),
    D = case NumBins == 1 of
            true  -> 0;
            false -> (Hi - Lo) / ((NumBins - 1) * 2)
        end,
    {Lo - D, Hi + D}.

%% The smallest 'float' such that 1 + E =/= 1.0.
%% STUB: found by thorough testing...
%%
%% See for yourself:
%%
%% 1> 1 + 1.1102230246251566637e-16 =/= 1.0.
%% true
%%
%% 2> 1 + 1.1102230246251566636999999999999999e-16 =/= 1.0.
%% false
-spec m_epsilon() -> float().
m_epsilon() -> 1.1102230246251566637e-16.
