-module(transform).

-export([
        %% Discrete cosine transform
          dct/1
        , dct_/1
        , idct/1
        , idct_/1
        %% Fast Fourier transform
        , fft/1
        , ifft/1
        ]).

-export_type([
             %% Type synonyms
              cd/0
             ]).

%% Fourier-related transformations of mathematical functions.
%%
%% These functions are written for simplicity and correctness, not
%% speed.  If you need a fast FFT implementation for your application,
%% you should strongly consider using a library of FFTW bindings
%% instead.


-type cd() :: complex:complex(). % Complex-type alias.

%% Discrete cosine transform (DCT-II)
-spec dct([float()]) -> [float()].
dct(XS) ->
    dct_worker(lists:map(fun (X) -> complex:make(X, 0) end, XS)).

%% Discrete cosine transform (DCT-II). Only real part of vector is
%% transformed, imaginary part is ignored.
-spec dct_([cd()]) -> [float()].
dct_(XS) ->
    dct_worker(lists:map(fun (X) -> complex:make(complex:real(X), 0) end, XS)).

-spec dct_worker([cd()]) -> [float()].
dct_worker(XS) ->
    Len = length(XS),
    Mul = fun complex:multiply/2,
    Exp = fun complex:exp/1,
    WeightsGen = fun (X) ->
                         %% 2 * exp ((0:+(-1)) * (x+1) * pi/(2*n))
                         Mul(Exp(Mul(complex:make(0,-1),
                                     (X+1) * math:pi() / (2*Len))), 2)
                 end,
    Weights = [complex:make(2) | [WeightsGen(X) || X <- lists:seq(0, Len-2)]],
    Backpermute = fun (IS) ->
                          lists:map(fun (I) -> lists:nth(I+1, XS) end, IS)
                  end,
    Interleaved = Backpermute(lists:seq(0, Len-2, 2) ++ lists:seq(Len-1, 1, -2)),
    lists:map(fun complex:real/1,
	      lists:zipwith(Mul, Weights, fft(Interleaved))).



%% Inverse discrete cosine transform (DCT-III). It's inverse of
%% 'dct' only up to scale parameter:
%%
%% > (idct . dct) x = (* lenngth x)
-spec idct([float()]) -> [float()].
idct(XS) ->
    idct_worker(lists:map(fun (X) -> complex:make(X, 0) end, XS)).

%% Inverse discrete cosine transform (DCT-III). Only real part of vector is
%% transformed, imaginary part is ignored.
-spec idct_([cd()]) -> [float()].
idct_(XS) ->
    idct_worker(lists:map(fun (X) -> complex:make(complex:real(X), 0) end, XS)).

-spec idct_worker([cd()]) -> [float()].
idct_worker(XS) ->
    Len = length(XS),
    Mul = fun complex:multiply/2,
    Exp = fun complex:exp/1,
    WeightsGen = fun (X) ->
                         %% 2 * n * exp ((0:+1) * (x+1) * pi/(2*n))
                         Mul(Exp(Mul(complex:make(0,1),
                                     (X+1) * math:pi() / (2*Len))), 2 * Len)
                 end,
    Weights = [complex:make(Len) | [WeightsGen(X) || X <- lists:seq(0, Len-2)]],
    IFFTed = ifft(lists:zipwith(Mul, Weights, XS)),
    Vals = lists:map(fun complex:real/1, IFFTed),
    Interleave = fun (Z) when Z rem 2 == 0 -> lists:nth(halve(Z) + 1,   Vals);
                     (Z)                   -> lists:nth(Len - halve(Z), Vals)
                 end,
    [ Interleave(I) || I <- lists:seq(0, Len - 1) ].

%%STUB: ifft . fft in Erlang =/= in Haskell

%% Inverse fast Fourier transform.
-spec ifft([cd()]) -> [cd()].
ifft(XS) ->
    FFTed = fft(lists:map(fun complex:conjugate/1, XS)),
    Len = length(XS),
    lists:map(fun (Z) -> complex:divide(complex:conjugate(Z), Len) end, FFTed).

%% Radix-2 decimation-in-time fast Fourier transform.
-spec fft([cd()]) -> [cd()].
fft(V) ->
    Len = length(V),
    T = ets:new(tab, [named_table]),
    true = ets:insert(T, lists:zip(lists:seq(0, Len-1), V)), % mutable.
    ok = mfft({T, Len}),
    {_, VS} = lists:unzip(lists:sort(fun ({A1, _}, {A2, _}) -> A1 < A2 end,
				     ets:tab2list(T))), % freeze.
    true = ets:delete(T),
    VS.

-spec mfft({ets:tab(), integer()}) -> ok.
mfft({VecTab, Len}) ->
    M = spec_functions:log2(Len),
    case (1 bsl M) =/= Len of
	true  -> exit({?MODULE, mfft, "bad vector size"});
	false -> mfft__bitReverse({VecTab, Len, M}, 0, 0)
    end.

-spec mfft__bitReverse({ets:tab(), integer(), integer()},
                       integer(), integer()) -> ok.
mfft__bitReverse({VecTab, Len, M}, I, _) when I == Len - 1 ->
    mfft__stage({VecTab, Len, M}, 0, 1);
mfft__bitReverse({VecTab, Len, M}, I, J)                   ->
    ok = case I < J of
             true  -> [{I, IElem}] = ets:lookup(VecTab, I),
                      [{J, JElem}] = ets:lookup(VecTab, J),
                      true = ets:insert(VecTab, [{J,IElem}, {I,JElem}]),
                      ok; % swap.
             false -> ok
         end,
    mfft__inner({VecTab, Len, M, I}, Len bsr 1, J).

-spec mfft__inner({ets:tab(), integer(), integer(), integer()},
                  integer(), integer()) -> ok.
mfft__inner({VecTab, Len, M, I}, K, L) when K =< L ->
    mfft__inner({VecTab, Len, M, I}, K bsr 1, L - K);
mfft__inner({VecTab, Len, M, I}, K, L)             ->
    mfft__bitReverse({VecTab, Len, M}, I+1, L+K).

-spec mfft__stage({ets:tab(), integer(), integer()}, integer(), integer())
                 -> ok.
mfft__stage({_, _, M}, M, _)  -> ok;
mfft__stage({VecTab, Len, M}, L, L1) ->
    L2 = L1 bsl 1,
    E = -6.283185307179586 / L2,
    mfft__flight({VecTab, Len, M, L, L1, L2, E}, 0, 0.0).

-spec mfft__flight({ets:tab(), integer(), integer(), integer(), integer(),
                    integer(), float()}, integer(), float()) -> ok.
mfft__flight({VecTab, Len, M, L, L1, L2, _}, L1, _) ->
    mfft__stage({VecTab, Len, M}, L+1, L2);
mfft__flight({VecTab, Len, M, L, L1, L2, E}, J, A) ->
    mfft__butterfly({VecTab, Len, M, A, E, J, L, L1, L2}, J).

-spec mfft__butterfly({ets:tab(), integer(), integer(), float(), float(),
                       integer(), integer(), integer(), integer()}, integer())
                     -> ok.
mfft__butterfly({VecTab, Len, M, A, E, J, L, L1, L2}, I) when I >= Len ->
   mfft__flight({VecTab, Len, M, L, L1, L2, E}, J+1, A+E);
mfft__butterfly({VecTab, Len, M, A, E, J, L, L1, L2}, I) ->
    I1 = I + L1,
    [{I1, Z}] = ets:lookup(VecTab, I1),
    {Xi1, Yi1} = {complex:real(Z), complex:imag(Z)},
    {C, S} = {math:cos(A), math:sin(A)},
    D = complex:make(C * Xi1 - S * Yi1, S * Xi1 + C * Yi1),
    [{I, Ci}] = ets:lookup(VecTab, I),
    true = ets:insert(VecTab, [{I1, complex:subtract(Ci, D)},
                               {I, complex:add(Ci, D)}]),
    mfft__butterfly({VecTab, Len, M, A, E, J, L, L1, L2}, I + L2).

-spec halve(integer()) -> integer().
halve(I) ->
    I bsr 1.
