-module(sample_kerneldensity).

-export([
        %% Estimation functions
          kde/2
        , kde_/4
         , kde___go/3
        ]).

%% Kernel density estimation. This module provides a fast, robust,
%% non-parametric way to estimate the probability density function of
%% a sample.
%%
%% This estimator does not use the commonly employed "Gaussian rule
%% of thumb".  As a result, it outperforms many plug-in methods on
%% multimodal samples with widely separated modes.


%% Gaussian kernel density estimator for one-dimensional data, using
%% the method of Botev et al.
%%
%% The result is a pair of vectors, containing:
%%
%% * The coordinates of each mesh point.  The mesh interval is chosen
%%   to be 20% larger than the range of the sample.  (To specify the
%%   mesh interval, use 'kde_'.)
%%
%% * Density estimates at each mesh point.
-spec kde(integer(), [float()]) -> {[float()], [float()]}.
kde(N0, XS) ->
    {Lo, Hi} = function:min_max(XS),
    Range = case length(XS) =< 1 of
                true  -> 1;      % Unreasonable guess
                false -> Hi - Lo
            end,
    kde_(N0, Lo - Range / 10, Hi + Range / 10, XS).

%% Gaussian kernel density estimator for one-dimensional data, using
%% the method of Botev et al.
%%
%% The result is a pair of vectors, containing:
%%
%% * The coordinates of each mesh point.
%%
%% * Density estimates at each mesh point.
-spec kde_(integer(), float(), float(), [float()]) -> {[float()], [float()]}.
kde_(_,  _, _, []) ->
    exit({?MODULE, kde_, "empty sample"});
kde_(N0, _, _, _) when N0 < 1 ->
    exit({?MODULE, kde_, "invalid number of points"});
kde_(N0, Min, Max, XS) ->
    N = function:next_highest_power_of_two(N0),
    R = Max - Min,
    Len = length(XS),
    Sqr = fun (X) when is_number(X) -> math:pow(X, 2) end,
    %% Compute mesh
    Gen  = fun (Z) -> Min + (R / (N - 1) * Z) end,
    Mesh = [Gen(I) || I <- lists:seq(0, N - 1)],
    %% Compute density
    H = lists:map(fun (X) -> X / Len end,
                  sample_histogram:histogram_(N, Min, Max, XS)),
    Sum = lists:sum(H),
    A   = transform:dct(lists:map(fun (X) -> X / Sum end, H)),
    A2V = lists:map(fun (X) -> Sqr(0.5 * X) end, tl(A)),
    IV  = lists:map(Sqr, lists:seq(1, N-1)),
    F0  = fun (Q, T) ->
                 G = fun (I, A2) ->
                             math:pow(I, Q) * A2
                                 * math:exp((-I) * Sqr(math:pi()) * T)
                     end,
                 2 * math:pow(math:pi(), Q * 2)
                     * lists:sum(lists:zipwith(G, IV, A2V))
         end,
    F1 = fun (X) ->
                 X - math:pow(Len * 2 * math:sqrt(math:pi())
                              * kde___go({F0, Len}, 6, F0(7, X)), -0.4)
         end,
    T_Star = math_rootfinding:from_root(0.28 * math:pow(Len, -0.4),
             math_rootfinding:ridders(1.0e-14, {0.0, 0.1}, F1)),
    F2 = fun (B, Z) ->
                 %% b * exp (sqr z * sqr pi * t_star * (-0.5))
                 B * math:exp(Sqr(Z) * Sqr(math:pi()) * T_Star * (-0.5))
         end,
    Density = lists:map(fun (Z) -> Z / (2 * R) end,
                        transform:idct(lists:zipwith(F2, A, lists:seq(0, N-1)))),
    {Mesh, Density}.

-spec kde___go({fun ((number(), number()) -> number()), integer()},
               pos_integer(), number()) -> number().
kde___go(_, 1, H) -> H;
kde___go({F, Len}, S, H) ->
    Const = (1 + math:pow(0.5, S+0.5)) / 3,
    Product = fun (XS) ->
                      lists:foldl(fun (X, Acc) -> X * Acc end, 1.0, XS)
              end,
    K0 = Product(lists:seq(1, 2*S-1, 2)) / math:sqrt(2 * math:pi()),
    Time = math:pow(2 * Const * K0 / Len / H, 2 / (3 + 2 * S)),
    kde___go({F, Len}, S-1, F(S, Time)).


%% References:
%%
%% Botev. Z.I., Grotowski J.F., Kroese D.P. (2010). Kernel density
%% estimation via diffusion. /Annals of Statistics/
%% 38(5):2916&#8211;2957. <http://arxiv.org/pdf/1011.2602>
