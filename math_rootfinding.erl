-module(math_rootfinding).

-export([
          from_root/2
        , ridders/3
        ]).

-export_type([
               root/1
             ]).

%% Erlang functions for finding the roots of mathematical functions.

%% The result of searching for a root of a mathematical function.
-type root(T) :: 'notBracketed'   %% The function does not have opposite signs
                                  %% when evaluated at the lower and upper
                                  %% bounds of the search.
               | 'searchFailed'   %% The search failed to converge to within the
                                  %% given error tolerance after the given
                                  %% number of iterations.
               | {'root', T}.     %% A root was successfully found.


%% Returns either the result of a search for a root, or the default
%% value if the search failed.
-spec from_root(T, root(T)) -> T.
from_root(_, {'root', A}) -> A;
from_root(A, _)           -> A.

%% Use the method of Ridders to compute a root of a function.
%%
%% The function must have opposite signs when evaluated at the lower
%% and upper bounds of the search (i.e. the root must be bracketed).
-type root_func() :: fun((float()) -> float()).

-spec ridders(float(), {float(), float()}, root_func()) -> root(float()).
ridders(Tol, {Lo, Hi}, F) ->
    Flo = F(Lo),
    Fhi = F(Hi),
    case Flo == 0 of
        true  -> {'root', Lo};
        false ->
            case Fhi == 0 of
                true  -> {'root', Hi};
                false ->
                    case Flo * Fhi > 0 of
                        true  -> 'notBracketed';
                        false -> ridders__wrap_go(Tol, F, Lo, Flo, Hi, Fhi, 0)
                    end
            end
    end.

-spec ridders__wrap_go(float(), root_func(), float(), float(), float(), float()
                       , integer()) -> root(float()).
ridders__wrap_go(Tol, F, A, FA, B, FB, I) ->
    D  = abs(B - A),
    DM = (B - A) * 0.5,
    M  = A + DM,
    FM = F(M),
    Signum = fun(Num) when Num >= 0 -> +1;
                (_) -> -1
             end,
    DN = Signum(FB - FA) * DM * FM / math:sqrt(FM * FM - FA * FB),
    N  = M - Signum(DN) * min(abs(DN), abs(DM) - 0.5 * Tol),
    FN = F(N),
    Within = function:within(1, A, B),
    ridders__go(Tol, F, {D, DM, M, FM, Signum, DN, N, FN, Within}, A, FA, B, FB, I).

-type sign_fun() :: fun((number()) -> -1 | 1).

-spec ridders__go(float(), root_func(), {float(), float(), float(), integer()
                  , sign_fun(), float(), float(), float(), boolean()}, float()
                  , float(), float(), float(), integer()) -> root(float()).
%% Root is bracketed within 1 ulp. No improvement could be made
ridders__go(_, _, {_, _, _, _, _, _, _, _, Within}, A, _, _, _, _) when Within ->
    {'root', A};
%% Root is found. Check that F(M) == 0 is nessesary to ensure that root is never
%% passed to 'ridders__go'.
ridders__go(_, _, {_, _, M, FM, _, _, _, _, _}, _, _, _, _, _) when FM == 0 ->
    {'root', M};
ridders__go(Tol, _, {D, _, _, _, _, _, N, FN, _}, _, _, _, _, _)
  when (FN == 0) or (D < Tol) ->
    {'root', N};
%% Too many iterations performed. Fail
ridders__go(_, _, _, _, _, _, _, I) when I > 100 ->
    'searchFailed';
%% Ridder's approximation coincide with one of old bounds. Revert to bisection.
ridders__go(Tol, F, {_, _, M, FM, _, _, N, _, _}, A, FA, B, _, I)
  when (((N == A) or (N == B)) and (FM * FA < 0)) ->
    ridders__wrap_go(Tol, F, A, FA, M, FM, I + 1);
ridders__go(Tol, F, {_, _, M, FM, _, _, N, _, _}, A, _, B, FB, I)
  when ((N == A) or (N == B)) ->
    ridders__wrap_go(Tol, F, M, FM, B, FB, I + 1);
%% Proceed as usual
ridders__go(Tol,F,{_, _, M, FM, _, _, N, FN, _},_,_,_,_,I) when FN*FM < 0 ->
    ridders__wrap_go(Tol, F, N, FN, M, FM, I + 1);
ridders__go(Tol,F,{_, _, _, _, _, _, N, FN, _},A,FA,_,_,I) when FN*FA < 0 ->
    ridders__wrap_go(Tol, F, A, FA, N, FN, I + 1);
ridders__go(Tol, F, {_, _, _, _, _, _, N, FN, _}, _, _, B, FB, I) ->
    ridders__wrap_go(Tol, F, N, FN, B, FB, I + 1).


%% References:
%%
%% * Ridders, C.F.J. (1979) A new algorithm for computing a single
%%   root of a real continuous function.
%%   /IEEE Transactions on Circuits and Systems/ 26:979&#8211;980.
