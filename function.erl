-module(function).

-export([
          min_max/1
        , sort/1
        , sort_by/2
        , indexed/1
        , indices/1
        , next_highest_power_of_two/1
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
