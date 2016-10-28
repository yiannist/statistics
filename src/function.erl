-module(function).

-export([
          min_max/1
        , sort/1
        , sort_by/2
        , indexed/1
        , indices/1
        , next_highest_power_of_two/1
        , within/3
        ]).

-spec sort([T]) -> [T].
sort(L) ->
    lists:sort(L).

-type comp_func(T) :: fun((T, T) -> boolean()).

-spec sort_by(comp_func(T), [T]) -> [T].
sort_by(F, L) ->
    lists:sort(F, L).

-spec min_max([float()]) -> {float(), float()}.
min_max([H|T]) ->
    Go = fun(K, {Min, Max}) ->
                 {min(Min, K), max(Max, K)}
         end,
    lists:foldl(Go, {H, H}, T).

-spec indexed([T]) -> [{integer(), T}].
indexed(L) ->
    lists:zip(indices(L), L).

-spec indices([_]) -> [integer()].
indices(L) ->
    lists:seq(0, length(L) - 1).

-spec next_highest_power_of_two(integer()) -> integer().
next_highest_power_of_two(N) ->
    Go = fun (I, M) ->
                 M bor (M bsr I)
         end,
    1 + lists:foldl(Go, N - 1, [1, 2, 4, 8, 16, 32]).

%% Compare two 'Double' values for approximate equality, using
%% Dawson's method.
%%
%% The required accuracy is specified in ULPs (units of least
%% precision).  If the two numbers differ by the given number of ULPs
%% or less, this function returns 'true'.
-spec within(integer(), float(), float()) -> boolean().
within(Ulps, A, B) ->
    %% XXX: Translate a float to integer! Don't try this at home...
    <<Ai0:64/integer-signed>> = <<A/float>>,
    <<Bi0:64/integer-signed>> = <<B/float>>,
    Big = 16#8000000000000000,
    Ai = case Ai0 < 0 of
             true  -> Big - Ai0;
             false -> Ai0
         end,
    Bi = case Bi0 < 0 of
             true  -> Big - Bi0;
             false -> Bi0
         end,
    abs(Ai - Bi) =< Ulps.


%% References:
%%
%% Approximate floating point comparison, based on Bruce Dawson's
%% "Comparing floating point numbers":
%% <http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm>
