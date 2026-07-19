# Integral electronic proof supplement: the SU(N) adjoint Q_3 branch

This document is part of the proof package for `paper.tex`, not a progress
note.  It records the replacement for the invalid two-moment Choquet route
and supplies the endpoint arguments used by Theorem `thm:su-ratio`.  The
original offset-by-offset transition formulas in Propositions 24--62 are
retained for audit provenance but have been superseded in the load-bearing
chain by the uniform transition-strip lemma in `paper.tex`.  The temporary
named ratio assumptions below are local
interfaces: Propositions 11, 15, 19 and Corollary 69 discharge all of them,
and Theorem 71 contains no unproved assumption.  The clean-room source
manifest authenticates this exact supplement.

## Proposition 1 (failure of the two-moment cone route)

Let
\[
\nu=\frac{32}{297}\delta_{1/40}
  +\frac{37}{63}\delta_{2/5}
  +\frac{634}{2079}\delta_{5/2}.
\]
Then \(\nu\) is a probability measure on \([0,4]\) satisfying
\[
\int x\,d\nu(x)=1,\qquad \int x^2\,d\nu(x)=2,
\]
but
\[
\iint (x-y)^2(x+y-2)^9\,d\nu(x)d\nu(y)
=-\frac{587425016412531}{1310720000000000}<0.
\]

Proof.
Direct rational arithmetic gives
\[
\frac{32}{297}+\frac{37}{63}+\frac{634}{2079}=1,
\]
\[
\frac{32}{297}\frac1{40}
 +\frac{37}{63}\frac25
 +\frac{634}{2079}\frac52=1,
\]
and
\[
\frac{32}{297}\frac1{40^2}
 +\frac{37}{63}\frac{4}{25}
 +\frac{634}{2079}\frac{25}{4}=2.
\]
Expanding the finite double sum over the nine atom-pairs gives the displayed
negative value. Therefore the claim that first two moments alone imply
Q_3-nonnegativity is false. \(\square\)

## Lemma 2 (stable trace moments)

Let \(U\) be Haar-distributed in \(SU(N)\), and put \(X=|\operatorname{tr}U|^2\).
For every integer \(k\) with \(0\le k\le N\),
\[
E[X^k]=k!.
\]

Proof.
The function \(U\mapsto |\operatorname{tr}U|^{2k}\) is invariant under scalar
phases, so its \(SU(N)\) and \(U(N)\) Haar integrals agree. Schur-Weyl duality
gives
\[
(\mathbb C^N)^{\otimes k}
 \cong \bigoplus_{\lambda\vdash k,\ \ell(\lambda)\le N}
 S_\lambda(\mathbb C^N)\otimes \operatorname{Sp}_\lambda .
\]
Hence
\[
E_{U(N)}[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le N}
(\dim\operatorname{Sp}_\lambda)^2.
\]
If \(k\le N\), every partition of \(k\) has length at most \(N\). The
regular representation identity
\[
\sum_{\lambda\vdash k}(\dim\operatorname{Sp}_\lambda)^2=k!
\]
then gives the claim. \(\square\)

## Lemma 3 (stable adjoint moments)

For \(0\le r\le N\),
\[
E[(\chi_{\rm adj}(U))^r]=!r,
\]
where \(!r=\sum_{j=0}^{r}(-1)^{r-j}\binom rj j!\) is the derangement number.

Proof.
For \(SU(N)\),
\[
\chi_{\rm adj}(U)=|\operatorname{tr}U|^2-1=X-1.
\]
Using Lemma 2 for every \(j\le r\le N\),
\[
E[(X-1)^r]
=\sum_{j=0}^{r}(-1)^{r-j}\binom rj E[X^j]
=\sum_{j=0}^{r}(-1)^{r-j}\binom rj j!
=!r.
\]
\(\square\)

## Lemma 4 (positive stable Q_3 coefficients)

Let \(Y\) have the exponential distribution of mean \(1\), and set \(A=Y-1\).
Define
\[
q_n=\frac12 E[(A-B)^2(A+B)^n],
\]
where \(A,B\) are independent copies. Then \(q_n>0\) for every \(n\ge0\).

Proof.
The moment generating function of \(A\) is
\[
M(z)=E[e^{zA}]=\frac{e^{-z}}{1-z}.
\]
For every \(n\ge0\),
\[
q_n=n![z^n]\{M''(z)M(z)-M'(z)^2\}.
\]
Since
\[
M''M-(M')^2=M^2(\log M)''
=\frac{e^{-2z}}{(1-z)^4},
\]
write
\[
\frac{e^{-2z}}{(1-z)^4}=\sum_{n\ge0}q_n\frac{z^n}{n!}.
\]
Verification. Differentiating \(F(z)=e^{-2z}(1-z)^{-4}\) gives
\[
(1-z)F'(z)=(2+2z)F(z).
\]
Comparing exponential-generating coefficients gives, for \(n\ge1\),
\[
q_{n+1}=(n+2)q_n+2nq_{n-1},
\]
with \(q_0=1\) and \(q_1=2\). The recurrence has positive coefficients
and positive initial values, so \(q_n>0\) for every \(n\ge0\). \(\square\)

## Theorem 5 (stable-rank SU(N) adjoint Q_3 closure)

For \(SU(N)\), if \(N\ge n+2\), then
\[
Q_3^{SU(N),{\rm adj}}(n)\ge0.
\]
In fact \(Q_3^{SU(N),{\rm adj}}(n)/2=q_n>0\), with \(q_n\) from Lemma 4.

Proof.
The polynomial formula
\[
\frac12 Q_3(n)=
\sum_{k=0}^{n}\binom nk
\left(m_{k+2}m_{n-k}-m_{k+1}m_{n-k+1}\right)
\]
uses only the adjoint moments \(m_0,\ldots,m_{n+2}\). If \(N\ge n+2\),
Lemma 3 identifies all of these moments with the derangement moments
\(!0,\ldots,!(n+2)\), equivalently the moments of \(Y-1\) for
\(Y\sim{\rm Exp}(1)\). Lemma 4 gives the required positivity. \(\square\)

## Theorem 6 (SU(2) adjoint Q_3 closure)

For \(SU(2)\),
\[
Q_3^{SU(2),{\rm adj}}(n)\ge0
\]
for every \(n\ge0\).

Proof.
Parametrize \(SU(2)\) conjugacy classes by
\(\alpha\in[0,\pi]\), with normalized Haar density
\((2/\pi)\sin^2\alpha\,d\alpha\). The fundamental trace is
\(\tau(\alpha)=2\cos\alpha\), and the adjoint character is
\[
\chi_{\rm adj}(\alpha)=\tau(\alpha)^2-1=1+2\cos(2\alpha).
\]
For two independent angles \(\alpha,\beta\), put
\[
\theta=\alpha+\beta,\qquad \phi=\alpha-\beta .
\]
Then
\[
\chi_{\rm adj}(\alpha)+\chi_{\rm adj}(\beta)
=2+4\cos\theta\cos\phi,
\]
and
\[
(\chi_{\rm adj}(\alpha)-\chi_{\rm adj}(\beta))^2
=16\sin^2\theta\sin^2\phi .
\]
Also
\[
4\sin^2\alpha\sin^2\beta=(\cos\theta-\cos\phi)^2,
\qquad d\alpha\,d\beta=\frac12\,d\theta\,d\phi .
\]
Therefore the measure
\[
\frac12(\chi_{\rm adj}(\alpha)-\chi_{\rm adj}(\beta))^2
\,d\mu_{SU(2)}(\alpha)\,d\mu_{SU(2)}(\beta)
\]
pushes forward, after folding the \((\theta,\phi)\)-diamond by the
\(C_2\) Weyl symmetries, to
\[
\frac{8}{\pi^2}
(\cos\theta-\cos\phi)^2\sin^2\theta\sin^2\phi\,d\theta\,d\phi
\qquad (0\le\theta,\phi\le\pi).
\]
This is the normalized \(USp(4)\) Weyl conjugacy density from the Weyl
integration formula (Bröcker-tom Dieck, IV.(1.11)). The normalization is
fixed by the case \(n=0\), since
\[
\frac12E[(\chi_{\rm adj}(g)-\chi_{\rm adj}(h))^2]
=E[\chi_{\rm adj}^2]-E[\chi_{\rm adj}]^2=1
\]
by Schur orthogonality.

For \(A\in USp(4)\) with eigenangles \(\theta,\phi\), the coefficient
\(a_2(A)\) of the characteristic polynomial is
\[
a_2(A)=2+4\cos\theta\cos\phi
=\operatorname{tr}(\Lambda^2 V)(A),
\]
where \(V\) is the four-dimensional defining representation. Hence
\[
\frac12Q_3^{SU(2),{\rm adj}}(n)
=\int_{USp(4)} a_2(A)^n\,dA
=\int_{USp(4)} \operatorname{tr}((\Lambda^2 V)^{\otimes n})(A)\,dA .
\]
By Schur orthogonality (Bröcker-tom Dieck, Ch. II §4), the last integral is
the multiplicity of the trivial representation in
\((\Lambda^2V)^{\otimes n}\), hence is a nonnegative integer. This proves
the claim. \(\square\)

## Remaining SU(N) obligation

Theorem 5 closes the stable-rank region \(N\ge n+2\), and Theorem 6
closes \(N=2\). The remaining non-stable region is
\[
3\le N<n+2.
\]
The next hard estimate is therefore:

**Estimate SU-tail(N).** For each fixed \(N\ge3\), prove
\[
Q_3^{SU(N),{\rm adj}}(n)\ge0
\]
for all \(n>N-2\), using the exact Schur-Weyl partition formula
\[
E[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le N}
(\dim\operatorname{Sp}_\lambda)^2.
\]

## Proposition 7 (SU(3) finite exact tail through \(n=41\))

For \(SU(3)\),
\[
Q_3^{SU(3),{\rm adj}}(n)>0
\]
for every odd \(1\le n\le 41\). Together with pointwise positivity for even
\(n\), this closes \(SU(3)\) through \(n=41\).

Proof.
For \(SU(3)\), Schur-Weyl duality gives
\[
E[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le 3}
(\dim\operatorname{Sp}_\lambda)^2.
\]
The Specht dimensions are computed by the hook-length formula. Therefore
each adjoint moment
\[
m_r=\sum_{j=0}^r(-1)^{r-j}\binom rj
\sum_{\lambda\vdash j,\ \ell(\lambda)\le 3}
(\dim\operatorname{Sp}_\lambda)^2
\]
is an exact integer. Substituting these integers into
\[
\frac12 Q_3(n)=
\sum_{k=0}^{n}\binom nk
\left(m_{k+2}m_{n-k}-m_{k+1}m_{n-k+1}\right)
\]
gives the following exact positive values of \(Q_3(n)/2\):
\[
\begin{array}{c|rrrrrrrrrrr}
n&1&3&5&7&9&11&13&15&17&19&21\\
\hline
Q_3(n)/2&
2&28&654&21312&902860&47963656&3102699474&
237017650352&20745270011904&2026586247524336&
216415850589952420
\end{array}
\]
and
\[
\begin{array}{c|rrrrrrrrrr}
n&23&25&27&29&31&33&35&37&39&41\\
\hline
Q_3(n)/2&
24867505378433342912&
3038207135550489061624&
391094519170131227840848&
52665481345568936284426666&
7376887964890704207857676528&
1069805588210196920738159847080&
160011133433530908540032170124752&
24604097145845889887526308074133504&
3878710568820053596239667141694589952&
625419666608024802966499311190595760896
\end{array}
\]
All displayed integers are strictly positive. \(\square\)

## Estimate SU3-tail

After Theorem 5 and Proposition 7, the remaining SU(3) tail estimate is:
\[
Q_3^{SU(3),{\rm adj}}(n)\ge0
\qquad\text{for every odd }n\ge43.
\]

The next theorem supplies the trace recurrence used by the tail argument.

## Theorem 8 (SU(3) trace recurrence)

Let
\[
E_k:=E_{SU(3)}[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le 3}
(\dim\operatorname{Sp}_\lambda)^2.
\]
For every \(k\ge1\),
\[
(k+3)^2E_{k+1}=(10k^2+22k+9)E_k-9k^2E_{k-1}.
\]

Proof.
By Schur-Weyl duality,
\[
E_k=\sum_{\lambda\vdash k,\ \ell(\lambda)\le3}
(\dim\operatorname{Sp}_\lambda)^2.
\]
The Robinson-Schensted-Knuth correspondence identifies this sum with the
number \(u_3(k)\) of permutations of \(\{1,\ldots,k\}\) whose longest
increasing subsequence has length at most \(3\).

Gessel, *Symmetric functions and P-recursiveness*, J. Combin. Theory Ser. A
53 (1990), §7, Theorem 16 and the formula following (26), gives
\[
U_3(X):=\sum_{k\ge0}u_3(k)\frac{X^{2k}}{(k!)^2}
=\det(b_{|i-j|}(X))_{1\le i,j\le3},
\]
where \(b_j(X)=I_j(2X)\) and \(I_j\) is the modified Bessel function.
Thus, with \(t=X^2\) and
\[
W(t):=U_3(\sqrt t)=\sum_{k\ge0}E_k\frac{t^k}{(k!)^2},
\]
we may compute \(W\) from this determinant.

The determinant is
\[
U_3=b_0^3-2b_0b_1^2+2b_1^2b_2-b_0b_2^2.
\]
Let
\[
y(t):=b_0(\sqrt t)=I_0(2\sqrt t).
\]
The Bessel recurrence \(Xb_2=Xb_0-b_1\) and the derivative identity
\(\frac{d}{dX}b_0=2b_1\) give
\[
b_1=Xy'(t),\qquad b_2=y(t)-y'(t).
\]
Therefore
\[
W=y'(2y^2-yy'-2t(y')^2).
\]
The series \(y=\sum_{j\ge0}t^j/(j!)^2\) satisfies
\[
ty''+y'-y=0.
\]

Verification. Differentiating the displayed formula for \(W\), and using
\(ty''=y-y'\) and its derivatives to eliminate \(y'',y^{(3)},y^{(4)},y^{(5)}\),
gives
\[
\bigl[(\theta+1)(\theta+3)^2\partial_t
 -(10\theta^2+22\theta+9)+9t\bigr]W=0,
\qquad \theta=t\partial_t .
\]
Equivalently,
\[
t^3W^{(4)}+10t^2W^{(3)}+(23t-10t^2)W''
 +(9-32t)W'+(9t-9)W=0.
\]
Taking the coefficient of \(t^k/(k!)^2\) for \(k\ge1\) gives
\[
(k+3)^2E_{k+1}-(10k^2+22k+9)E_k+9k^2E_{k-1}=0,
\]
as required. \(\square\)

## Proposition 9 (transfer to adjoint moments)

The SU(3) adjoint moments
\[
m_k:=E_{SU(3)}[\chi_{\rm adj}(U)^k]
\]
satisfy, for every \(k\ge1\),
\[
(k+3)^2m_{k+1}=k(7k+11)m_k+8k(k+1)m_{k-1},
\qquad m_0=1,\quad m_1=0.
\]

Proof.
Let
\[
A(z)=\sum_{k\ge0}E_k\frac{z^k}{k!}.
\]
Multiplying Theorem 8 by \(z^k/k!\), summing over \(k\ge1\),
and using \(E_0=E_1=1\), gives the differential equation
\[
z^2A'''+(7z-10z^2)A''+(9-32z+9z^2)A'+(-9+9z)A=0.
\]
Since \(\chi_{\rm adj}=|\operatorname{tr}U|^2-1\), the adjoint-moment EGF is
\[
M(z)=\sum_{k\ge0}m_k\frac{z^k}{k!}=e^{-z}A(z).
\]
Substituting \(A=e^zM\) into the displayed equation gives
\[
z^2M'''+(7z-7z^2)M''+(9-18z-8z^2)M'-16zM=0.
\]
Verification. The coefficient of \(z^k/k!\) in this equation is
\[
\bigl(k(k-1)+7k+9\bigr)m_{k+1}
 -\bigl(7k(k-1)+18k\bigr)m_k
 -\bigl(8k(k-1)+16k\bigr)m_{k-1}=0.
\]
This simplifies to
\[
(k+3)^2m_{k+1}=k(7k+11)m_k+8k(k+1)m_{k-1}.
\]
The initial values are \(m_0=1\) and \(m_1=E_1-E_0=0\). \(\square\)

## Assumption SU3-RatioMono

Let \(m_k\) be the SU(3) adjoint moments and
\[
\rho_k:=m_{k+1}/m_k .
\]
Assume
\[
\rho_3=\rho_4=4,\qquad 4\le \rho_k\le \rho_{k+1}
\quad\text{for every }k\ge3.
\]

## Proposition 10 (conditional boundary-layer positivity)

Under Proposition 9 and Assumption SU3-RatioMono,
\[
Q_3^{SU(3),{\rm adj}}(n)>0
\]
for every odd \(n\ge 11\).

Proof.
For \(0\le j\le n\), put
\[
B_j=m_{j+2}m_{n-j}+m_jm_{n+2-j}-2m_{j+1}m_{n+1-j}.
\]
Then
\[
\frac12Q_3(n)=\frac12\sum_{j=0}^{n}\binom nj B_j .
\]
For \(3\le j\le n-3\), Assumption SU3-RatioMono gives
\[
\frac{B_j}{m_{j+1}m_{n+1-j}}
=\frac{\rho_{j+1}}{\rho_{n-j}}
 +\frac{\rho_{n+1-j}}{\rho_j}-2\ge0
\]
by AM-GM, since
\[
\frac{\rho_{j+1}\rho_{n+1-j}}{\rho_{n-j}\rho_j}\ge1.
\]
It remains to check the depth-2 boundary contribution
\[
L_2:=B_0+nB_1+\binom n2 B_2 .
\]
Here
\[
B_0=m_n+m_{n+2},\quad
B_1=2m_{n-1}-2m_n,\quad
B_2=8m_{n-2}+m_n-4m_{n-1}.
\]
Set
\[
S:=\rho_{n-2}=m_{n-1}/m_{n-2},\qquad
R:=\rho_{n-1}=m_n/m_{n-1}.
\]
Assumption SU3-RatioMono gives \(R\ge S\ge4\), so write
\[
S=4+s,\qquad R=4+s+t,\qquad s,t\ge0.
\]
Using the recurrence of Proposition 9 to eliminate \(m_{n+2}\), a direct
calculation gives
\[
L_2=
\frac{m_{n-2}}{2(n+3)^2(n+4)^2}P(n,s,t),
\]
where
\[
\begin{aligned}
P(n,s,t)=&
n^6s^2+n^6st+4n^6s+4n^6t+8n^6\\
&+9n^5s^2+9n^5st+24n^5s+36n^5t+56n^5\\
&+119n^4s^2+119n^4st+884n^4s+476n^4t+2104n^4\\
&+479n^3s^2+479n^3st+4256n^3s+1916n^3t+10120n^3\\
&+716n^2s^2+716n^2st+7184n^2s+2864n^2t+17088n^2\\
&+636ns^2+636nst+6528ns+2544nt+14784n\\
&+576s^2+576st+4608s+2304t+9216 .
\end{aligned}
\]
Verification. Expanding \(B_0+nB_1+\binom n2B_2\), substituting
\[
m_{n-1}=S m_{n-2},\qquad m_n=RS m_{n-2},
\]
and using
\[
m_{n+1}=\frac{n(7n+11)m_n+8n(n+1)m_{n-1}}{(n+3)^2},
\]
\[
m_{n+2}=
\frac{(n+1)(7n+18)m_{n+1}+8(n+1)(n+2)m_n}{(n+4)^2}
\]
produces exactly the displayed \(P(n,s,t)\).
All coefficients of \(P\) are nonnegative and the constant term is positive
for \(n\ge1\). Hence \(L_2>0\). The interior terms are nonnegative, so
\(Q_3(n)/2\ge L_2>0\). \(\square\)

## Proposition 11 (ratio monotonicity)

Assumption SU3-RatioMono follows from the recurrence of Proposition 9.

Proof.
Set
\[
U_k:=8-\rho_k .
\]
The recurrence of Proposition 9 gives, for \(k\ge2\),
\[
U_k=
\frac{k^2+37k+72-\dfrac{8k(k+1)}{8-U_{k-1}}}{(k+3)^2}.
\]
Define
\[
G_k:=\frac{32k+100}{(k+4)^2},\qquad
H_k:=\frac{32(k+11/5)}{(k+3)^2}.
\]
We prove
\[
G_k\le U_k\le H_k\qquad(k\ge5).
\]
The base case is exact:
\[
U_5=\frac{111}{32},\qquad
G_5=\frac{260}{81},\qquad
H_5=\frac{18}{5}.
\]
Assume \(G_k\le U_k\le H_k\) for some \(k\ge5\). The map
\[
\Phi_k(u):=
\frac{k^2+37k+72-\dfrac{8k(k+1)}{8-u}}{(k+3)^2}
\]
is decreasing for \(u<8\). Also \(H_k<4<8\), since
\[
4-H_k=\frac{4(5k^2-10k-43)}{5(k+3)^2}>0
\qquad(k\ge5).
\]
Therefore
\[
U_{k+1}=\Phi_{k+1}(U_k)\ge \Phi_{k+1}(H_k),
\]
and a direct expansion gives
\[
\Phi_{k+1}(H_k)-G_{k+1}
=
\frac{4(34k^4+288k^3+683k^2+278k-403)}
{(k+4)^2(k+5)^2(5k^2+10k+1)}
>0 .
\]
Similarly,
\[
U_{k+1}=\Phi_{k+1}(U_k)\le \Phi_{k+1}(G_k),
\]
and another direct expansion gives
\[
H_{k+1}-\Phi_{k+1}(G_k)
=
\frac{29k^2+91k+54}
{5(k+4)^2(2k^2+8k+7)}
>0 .
\]
This proves \(G_k\le U_k\le H_k\) for every \(k\ge5\).

It remains to compare \(U_k\) with the fixed point threshold for ratio
monotonicity. For \(k\ge3\), the inequality \(\rho_{k+1}\ge\rho_k\) is
equivalent to
\[
F_k(\rho_k)\le0,
\]
where
\[
F_k(r):=(k+4)^2r^2-(k+1)(7k+18)r-8(k+1)(k+2).
\]
Let \(r_k^*\) be the positive root of \(F_k\), and put
\(\epsilon_k^*:=8-r_k^*\). Substituting \(r=8-\epsilon\) gives
\[
(k+4)^2\epsilon^2-(9k^2+103k+238)\epsilon+(288k+864)=0.
\]
Verification. Its discriminant satisfies the exact identity
\[
\begin{aligned}
 &(9k^2+103k+238)^2-4(k+4)^2(288k+864)\\
 &\qquad=(9k^2+39k+38)^2+16(k-3)(k+2).
\end{aligned}
\]
For \(k\ge3\), the square root is at least \(9k^2+39k+38\). Hence
\[
\epsilon_k^*
\le
\frac{(9k^2+103k+238)-(9k^2+39k+38)}{2(k+4)^2}
=G_k .
\]
For \(k\ge5\), \(U_k\ge G_k\ge\epsilon_k^*\), so
\(\rho_k\le r_k^*\), equivalently \(\rho_{k+1}\ge\rho_k\). Also
\(U_k\le H_k<4\), so \(\rho_k\ge4\) for every \(k\ge5\). The remaining
initial values are
\[
\rho_3=\rho_4=4,\qquad \rho_5=\frac{145}{32}>4.
\]
Thus Assumption SU3-RatioMono holds. \(\square\)

## Corollary 12 (SU(3) tail closure)

For SU(3),
\[
Q_3^{SU(3),{\rm adj}}(n)\ge0
\]
for every \(n\ge0\).

Proof.
Theorem 8 gives Proposition 9. Proposition 11 then
discharges Assumption SU3-RatioMono, so Proposition 10 gives strict
positivity for every odd \(n\ge11\). Proposition 7 gives the odd cases
\(1\le n\le9\). Even \(n\) are pointwise nonnegative. \(\square\)

Thus the SU(3) tail branch of Estimate SU-tail(N) is closed. The remaining
finite-rank SU(N) branch is
\[
4\le N<n+2.
\]

## Estimate SU4-tail

The first unresolved branch after Corollary 12 is
\[
Q_3^{SU(4),{\rm adj}}(n)\ge0
\qquad\text{for every odd }n\ge5.
\]
The cases \(n=1,3\) are exact:
\[
Q_3^{SU(4),{\rm adj}}(1)/2=2,\qquad
Q_3^{SU(4),{\rm adj}}(3)/2=39.
\]

## Theorem 13 (SU(4) adjoint moment recurrence)

Let
\[
m_k:=E_{SU(4)}[\chi_{\rm adj}(U)^k].
\]
For every \(k\ge0\),
\[
\begin{aligned}
&(k^3+22k^2+161k+392)m_{k+4}\\
&\quad-(16k^3+230k^2+1054k+1524)m_{k+3}\\
&\quad+(10k^3-340k-750)m_{k+2}\\
&\quad+(72k^3+522k^2+1242k+972)m_{k+1}\\
&\quad+(45k^3+270k^2+495k+270)m_k=0.
\end{aligned}
\]

Proof.
Let
\[
E_k:=E_{SU(4)}[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le4}
(\dim\operatorname{Sp}_\lambda)^2 .
\]
By RSK, \(E_k\) is the number of permutations of \(\{1,\ldots,k\}\) whose
longest increasing subsequence has length at most \(4\). Gessel,
*Symmetric functions and P-recursiveness*, J. Combin. Theory Ser. A 53
(1990), §7, Theorem 16, gives
\[
U_4(X):=\sum_{k\ge0}E_k\frac{X^{2k}}{(k!)^2}
=\det(b_{|i-j|}(X))_{1\le i,j\le4},
\]
where \(b_j(X)=I_j(2X)\).

Put \(t=X^2\), \(W(t)=U_4(\sqrt t)\), and
\[
y(t):=b_0(\sqrt t)=I_0(2\sqrt t).
\]
The Bessel recurrence \(Xb_j=(1-j)b_{j-1}+Xb_{j-2}\) gives
\[
b_1=Xy',\qquad b_2=y-y',\qquad
b_3=Xy'-\frac{2(y-y')}{X}.
\]
The determinant then simplifies to
\[
W=
-\frac{
4t^2(y')^4-8ty^2(y')^2+8ty(y')^3-t(y')^4
+4y^4-8y^3y'+4y^2(y')^2}{t}.
\]
Since \(y=\sum_{j\ge0}t^j/(j!)^2\), it satisfies \(ty''+y'-y=0\).

Verification. Differentiating the displayed expression for \(W\), and
using \(ty''=y-y'\) and its derivatives to eliminate all derivatives
of \(y\) of order at least \(2\), gives
\[
\begin{aligned}
&\bigl[(\theta+1)(\theta+2)(\theta+5)^2(\theta+6)\partial_t^2\\
&\quad -2(10\theta^3+91\theta^2+255\theta+214)(\theta+1)\partial_t\\
&\quad +64(\theta+1)^2(\theta+2)\bigr]W=0,
\qquad \theta=t\partial_t .
\end{aligned}
\]
Taking the coefficient of \(t^k/(k!)^2\) gives, for \(k\ge0\),
\[
(k+5)^2(k+6)E_{k+2}
-2(10k^3+91k^2+255k+214)E_{k+1}
+64(k+1)^2(k+2)E_k=0.
\]

Now
\[
m_k=\sum_{j=0}^k(-1)^{k-j}\binom kj E_j,
\]
because \(\chi_{\rm adj}=|\operatorname{tr}U|^2-1\). Equivalently, for
ordinary generating functions,
\[
\sum_{k\ge0}m_kx^k
=\frac{1}{1+x}\sum_{k\ge0}E_k\left(\frac{x}{1+x}\right)^k .
\]
Substituting this relation into the ordinary-generating-function form of
the displayed recurrence for \(E_k\) gives
\[
x^2(1+x)^2\mathcal L_{\rm adj}G=0,
\]
where \(G(x)=\sum_{k\ge0}m_kx^k\), and \(\mathcal L_{\rm adj}G=0\) is the
ordinary-generating-function equation whose coefficient of \(x^k\) is the
recurrence displayed in the theorem. Comparing coefficients gives the
claim. Verification. The symbolic residual of the identity
\(\mathcal L_{\rm trace}F=x^2(1+x)^2\mathcal L_{\rm adj}G\) is zero, and
the displayed recurrence agrees with exact Schur-Weyl partition arithmetic
for \(0\le k\le20\). \(\square\)

## Assumption SU4-RatioMono

Let
\[
\rho_k:=m_{k+1}/m_k .
\]
Assume
\[
\rho_3=\frac92,\qquad \frac92\le\rho_k\le\rho_{k+1}
\quad\text{for every }k\ge3.
\]

## Proposition 14 (SU(4) tail conditional on ratio monotonicity)

Under Assumption SU4-RatioMono,
\[
Q_3^{SU(4),{\rm adj}}(n)>0
\]
for every odd \(n\ge5\).

Proof.
For \(0\le j\le n\), put
\[
B_j=m_{j+2}m_{n-j}+m_jm_{n+2-j}-2m_{j+1}m_{n+1-j}.
\]
Then
\[
\frac12Q_3(n)=\frac12\sum_{j=0}^{n}\binom njB_j.
\]
For \(3\le j\le n-3\), Assumption SU4-RatioMono gives
\[
\frac{B_j}{m_{j+1}m_{n+1-j}}
=\frac{\rho_{j+1}}{\rho_{n-j}}
 +\frac{\rho_{n+1-j}}{\rho_j}-2\ge0
\]
by AM-GM, since
\[
\frac{\rho_{j+1}\rho_{n+1-j}}{\rho_{n-j}\rho_j}\ge1.
\]
It remains to check the depth-2 boundary contribution
\[
L_2:=B_0+nB_1+\binom n2 B_2.
\]
For \(SU(4)\), the initial adjoint moments are
\[
m_0=1,\qquad m_1=0,\qquad m_2=1,\qquad m_3=2,\qquad m_4=9,
\]
so
\[
B_0=m_n+m_{n+2},\quad
B_1=2m_{n-1}-2m_n,\quad
B_2=9m_{n-2}+m_n-4m_{n-1}.
\]
Set
\[
\rho_{n-2}=\frac92+a,\qquad
\rho_{n-1}=\frac92+a+b,\qquad
\rho_n=\frac92+a+b+c,
\]
where \(a,b,c\ge0\) by Assumption SU4-RatioMono. Applying
Theorem 13 at \(k=n-2\) to eliminate \(m_{n+2}\), and using
\[
m_{n-1}=\rho_{n-2}m_{n-2},\quad
m_n=\rho_{n-1}\rho_{n-2}m_{n-2},\quad
m_{n+1}=\rho_n\rho_{n-1}\rho_{n-2}m_{n-2},
\]
gives
\[
L_2=\frac{m_{n-2}}{8(n+5)^2(n+6)}P_4(n,a,b,c),
\]
where
\[
\begin{aligned}
P_4={}&
128a^3n^3+1072a^3n^2+2608a^3n+1664a^3\\
&+256a^2bn^3+2144a^2bn^2+5216a^2bn+3328a^2b\\
&+128a^2cn^3+1072a^2cn^2+2608a^2cn+1664a^2c\\
&+4a^2n^5+44a^2n^4+1676a^2n^3+13980a^2n^2
  +34648a^2n+24864a^2\\
&+128ab^2n^3+1072ab^2n^2+2608ab^2n+1664ab^2\\
&+128abcn^3+1072abcn^2+2608abcn+1664abc\\
&+4abn^5+44abn^4+2252abn^3+18804abn^2
  +46384abn+32352ab\\
&+1152acn^3+9648acn^2+23472acn+14976ac\\
&+20an^5+172an^4+5884an^3+60296an^2+158052an+122688a\\
&+576b^2n^3+4824b^2n^2+11736b^2n+7488b^2\\
&+576bcn^3+4824bcn^2+11736bcn+7488bc\\
&+18bn^5+198bn^4+4950bn^3+41202bn^2+103104bn+78192b\\
&+2592cn^3+21708cn^2+52812cn+33696c\\
&+45n^5+423n^4+6327n^3+88263n^2+242226n+200232.
\end{aligned}
\]
Verification. Expanding \(L_2\), substituting the three displayed ratio
relations, and substituting the \(k=n-2\) instance of Theorem 13 produces
exactly the displayed numerator \(P_4\). Every coefficient in \(P_4\) is
nonnegative, and its constant term is positive for \(n\ge1\). Hence
\(L_2>0\). The interior terms are nonnegative, so \(Q_3(n)/2\ge L_2>0\).
\(\square\)

## Theorem 15 (SU(4) ratio monotonicity)

Assumption SU4-RatioMono holds.

Proof.
The initial adjoint moments used below are
\[
m_3=2,\quad m_4=9,\quad m_5=43,\quad m_6=245,\quad m_7=1557,\quad
m_8=10829,\quad m_9=80958,\quad m_{10}=642501.
\]
Hence
\[
\rho_3=\frac92,\quad \rho_4=\frac{43}{9},\quad
\rho_5=\frac{245}{43},\quad \rho_6=\frac{1557}{245},\quad
\rho_7=\frac{10829}{1557},\quad
\rho_8=\frac{80958}{10829},\quad \rho_9=\frac{214167}{26986},
\]
and
\[
\rho_4-\rho_3=\frac5{18},\quad
\rho_5-\rho_4=\frac{356}{387},\quad
\rho_6-\rho_5=\frac{6926}{10535},\quad
\rho_7-\rho_6=\frac{228856}{381465}.
\]
It remains to prove \(\rho_{k+1}\ge\rho_k\) for \(k\ge7\).

Set
\[
U_k:=15-\rho_k .
\]
For \(k\ge3\), Theorem 13 gives
\[
\rho_{k+3}
=\frac{B_k}{A_k}
-\frac{C_k\rho_k\rho_{k+1}+D_k\rho_k+E_k}
{A_k\rho_k\rho_{k+1}\rho_{k+2}},
\]
where
\[
\begin{aligned}
A_k&=k^3+22k^2+161k+392,\\
B_k&=16k^3+230k^2+1054k+1524,\\
C_k&=10k^3-340k-750,\\
D_k&=72k^3+522k^2+1242k+972,\\
E_k&=45k^3+270k^2+495k+270.
\end{aligned}
\]
Equivalently,
\[
U_{k+3}=\Phi_k(U_k,U_{k+1},U_{k+2}),
\]
where
\[
\begin{aligned}
\Phi_k(u,v,w):={}&15-\frac{B_k}{A_k}
 +\frac{C_k}{A_k(15-w)}
 +\frac{D_k}{A_k(15-v)(15-w)}\\
&+\frac{E_k}{A_k(15-u)(15-v)(15-w)} .
\end{aligned}
\]
Define
\[
G_k:=\frac{1125}{10k+71},\qquad
H_k:=\frac{1125}{10k+66}.
\]
For \(j\ge7\), \(0<G_j<H_j\le H_7=1125/136<15\). Since
\(A_k,C_k,D_k,E_k\) are positive for \(k\ge7\), \(\Phi_k\) is increasing
in \(u,v,w\) on the box \(0\le u\le H_k\), \(0\le v\le H_{k+1}\),
\(0\le w\le H_{k+2}\).
We prove
\[
G_k\le U_k\le H_k\qquad(k\ge7).
\]
The base cases are exact:
\[
\begin{array}{c|cc}
k&U_k-G_k&H_k-U_k\\ \hline
7&\dfrac{4847}{73179}&\dfrac{48089}{211752}\\[4pt]
8&\dfrac{120402}{1635179}&\dfrac{286983}{1581034}\\[4pt]
9&\dfrac{331053}{4344746}&\dfrac{103677}{701636}.
\end{array}
\]
Assume \(G_j\le U_j\le H_j\) for \(j=k,k+1,k+2\), with \(k\ge7\).
Since \(\Phi_k\) is increasing,
\[
\Phi_k(G_k,G_{k+1},G_{k+2})
\le U_{k+3}\le
\Phi_k(H_k,H_{k+1},H_{k+2}).
\]
Verification. Direct expansion gives
\[
\Phi_k(G_k,G_{k+1},G_{k+2})-G_{k+3}
=\frac{P_-(k-7)}
{40(k+7)^2(k+8)(5k-2)(5k+3)(5k+8)(10k+101)},
\]
where
\[
\begin{aligned}
P_-(s)={}&2432000s^5+100233600s^4+1609345435s^3\\
&+12643460157s^2+48855354258s+74691137880,
\end{aligned}
\]
and
\[
H_{k+3}-\Phi_k(H_k,H_{k+1},H_{k+2})
=\frac{P_+(k-7)}
{10(k+7)^2(k+8)(5k+48)(10k-9)(10k+1)(10k+11)},
\]
where
\[
\begin{aligned}
P_+(s)={}&128000s^5+24662400s^4+750722465s^3\\
&+8917574583s^2+46522820852s+89208809220.
\end{aligned}
\]
Every coefficient in \(P_-\) and \(P_+\) is positive, and the denominators
are positive for \(k\ge7\). This proves the induction.

Finally, for \(k\ge7\),
\[
U_{k+1}\le H_{k+1}=\frac{1125}{10k+76}
<\frac{1125}{10k+71}=G_k\le U_k.
\]
Thus \(\rho_{k+1}>\rho_k\) for \(k\ge7\). Together with the exact initial
differences above, this proves
\[
\rho_3=\frac92,\qquad \frac92\le\rho_k\le\rho_{k+1}
\quad(k\ge3).
\]
\(\square\)

## Corollary 16 (SU(4) adjoint Q_3 closure)

For every \(n\ge0\),
\[
Q_3^{SU(4),{\rm adj}}(n)\ge0.
\]

Proof.
Theorem 15 discharges Assumption SU4-RatioMono, so Proposition 14 gives
strict positivity for every odd \(n\ge5\). The exact cases \(n=1,3\) are
listed in Estimate SU4-tail. Even \(n\) are pointwise nonnegative.
\(\square\)

Thus the SU(4) tail branch is closed. The next finite-rank SU(N) branch is
\[
5\le N<n+2.
\]

## Estimate SU5-tail

The next branch is
\[
Q_3^{SU(5),{\rm adj}}(n)\ge0
\qquad\text{for every odd }n\ge5.
\]
The cases \(n=1,3\) are exact:
\[
Q_3^{SU(5),{\rm adj}}(1)/2=4,\qquad
Q_3^{SU(5),{\rm adj}}(3)/2=80.
\]

## Theorem 17 (SU(5) adjoint moment recurrence)

Let
\[
m_k:=E_{SU(5)}[\chi_{\rm adj}(U)^k].
\]
For every \(k\ge0\),
\[
\begin{aligned}
&192(k+1)(k+2)(k+3)(k+6)m_k\\
&\quad+16(k+2)(k+3)(22k^2+184k+309)m_{k+1}\\
&\quad+(k+3)(129k^3+1245k^2+2570k-1112)m_{k+2}\\
&\quad-(k+3)(30k^3+642k^2+4487k+10196)m_{k+3}\\
&\quad+(k+8)^2(k+10)^2m_{k+4}=0.
\end{aligned}
\]

Proof.
Let
\[
E_k:=E_{SU(5)}[|\operatorname{tr}U|^{2k}]
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le5}
(\dim\operatorname{Sp}_\lambda)^2 .
\]
As in Theorem 13, Gessel's determinant gives
\[
U_5(X):=\sum_{k\ge0}E_k\frac{X^{2k}}{(k!)^2}
=\det(b_{|i-j|}(X))_{1\le i,j\le5},
\qquad b_j(X)=I_j(2X).
\]
Put \(t=X^2\), \(W(t)=U_5(\sqrt t)\), and
\[
y(t):=I_0(2\sqrt t).
\]
The Bessel recurrence gives
\[
b_1=Xy',\quad b_2=y-y',\quad
b_3=Xy'-\frac{2(y-y')}{X},\quad
b_4=y-4y'+\frac{6(y-y')}{t}.
\]
The determinant simplifies to
\[
W=-\frac{4}{t^2}
\bigl(-2t(y')^2-(y')^2-yy'+2y^2\bigr)
\bigl(3t(y')^3-4ty(y')^2+y(y')^2-5y^2y'+4y^3\bigr).
\]
Since \(ty''+y'-y=0\), differentiating this expression and eliminating
higher derivatives gives
\[
\begin{aligned}
&(k+7)^2(k+9)^2E_{k+3}
-(35k^4+742k^3+5631k^2+17932k+19941)E_{k+2}\\
&\quad +(k+2)^2(259k^2+2176k+4242)E_{k+1}
-225(k+1)^2(k+2)^2E_k=0 .
\end{aligned}
\]
Verification. The symbolic residual of this differential identity is zero.

Now
\[
m_k=\sum_{j=0}^k(-1)^{k-j}\binom kj E_j,
\]
because \(\chi_{\rm adj}=|\operatorname{tr}U|^2-1\). Substituting this
binomial transform into the displayed trace recurrence and comparing
ordinary-generating-function coefficients gives the recurrence displayed in
the theorem. Verification. The displayed recurrence agrees with exact
Gessel-determinant coefficients for \(0\le k\le65\). \(\square\)

## Assumption SU5-RatioMono

Let \(m_k\) be the SU(5) adjoint moments and
\[
\rho_k:=m_{k+1}/m_k .
\]
Assume
\[
\rho_3=\frac92,\qquad \frac92\le\rho_k\le\rho_{k+1}
\quad\text{for every }k\ge3.
\]

## Proposition 18 (SU(5) tail conditional on ratio monotonicity)

Under Assumption SU5-RatioMono,
\[
Q_3^{SU(5),{\rm adj}}(n)>0
\]
for every odd \(n\ge5\).

Proof.
For \(0\le j\le n\), put
\[
B_j=m_{j+2}m_{n-j}+m_jm_{n+2-j}-2m_{j+1}m_{n+1-j}.
\]
Then
\[
\frac12Q_3(n)=\frac12\sum_{j=0}^{n}\binom njB_j .
\]
For \(3\le j\le n-3\), Assumption SU5-RatioMono gives
\[
\frac{B_j}{m_{j+1}m_{n+1-j}}
=\frac{\rho_{j+1}}{\rho_{n-j}}
 +\frac{\rho_{n+1-j}}{\rho_j}-2\ge0
\]
by AM-GM. The depth-2 boundary contribution is
\[
L_2:=B_0+nB_1+\binom n2B_2 .
\]
Since \(m_0=1\), \(m_1=0\), \(m_2=1\), \(m_3=2\), and \(m_4=9\),
\[
B_0=m_n+m_{n+2},\quad
B_1=2m_{n-1}-2m_n,\quad
B_2=9m_{n-2}+m_n-4m_{n-1}.
\]
Set
\[
\rho_{n-2}=\frac92+a,\qquad
\rho_{n-1}=\frac92+a+b,\qquad
\rho_n=\frac92+a+b+c,
\]
where \(a,b,c\ge0\) by Assumption SU5-RatioMono. Applying Theorem 17 at
\(k=n-2\) to eliminate \(m_{n+2}\), and using the three ratio relations,
gives
\[
L_2=\frac{m_{n-2}}{8(n+6)^2(n+8)^2}P_5(n,a,b,c).
\]
Verification. Expanding \(L_2\), substituting the ratio relations, and
substituting the \(k=n-2\) instance of Theorem 17 gives the numerator
\(P_5\). In the polynomial ring \(\mathbb Q[a,b,c,s]\), with \(s=n-5\),
the shifted polynomial \(P_5(s+5,a,b,c)\) has \(80\) monomials and every
coefficient is positive; the coefficient minimum is \(4\). Its constant
part is
\[
45s^6+2313s^5+32700s^4+353898s^3+4615728s^2
+33865668s+87594336.
\]
Thus \(P_5(n,a,b,c)>0\) for \(n\ge5\). Hence \(L_2>0\), and the interior
terms are nonnegative, so \(Q_3(n)/2\ge L_2>0\). \(\square\)

## Theorem 19 (SU(5) ratio monotonicity)

Assumption SU5-RatioMono holds.

Proof.
The initial ratios are
\[
\rho_3=\frac92,\qquad \rho_4=\frac{44}{9},\qquad
\rho_5=6,\qquad \rho_6=\frac{76}{11},
\]
with
\[
\rho_4-\rho_3=\frac7{18},\qquad
\rho_5-\rho_4=\frac{10}{9},\qquad
\rho_6-\rho_5=\frac{10}{11}.
\]
It remains to prove \(\rho_{k+1}\ge\rho_k\) for \(k\ge6\).

Set
\[
U_k:=24-\rho_k .
\]
Theorem 17 gives \(U_{k+3}=\Psi_k(U_k,U_{k+1},U_{k+2})\), where
\[
\begin{aligned}
\Psi_k(u,v,w):={}&24+\frac{q_{3,k}}{q_{4,k}}
+\frac{q_{2,k}}{q_{4,k}(24-w)}
+\frac{q_{1,k}}{q_{4,k}(24-v)(24-w)}\\
&+\frac{q_{0,k}}{q_{4,k}(24-u)(24-v)(24-w)}
\end{aligned}
\]
and \(q_{j,k}\) are the five polynomial coefficients in Theorem 17.

Define
\[
G_k:=\frac{2880}{10k+109},\qquad H_k:=\frac{288}{k+10}.
\]
For \(j\ge6\), \(0<G_j<H_j\le H_6=18<24\). Since
\(q_{0,k},q_{1,k},q_{2,k},q_{4,k}\) are positive for \(k\ge6\), \(\Psi_k\)
is increasing in \(u,v,w\) on the box \(0\le u\le H_k\),
\(0\le v\le H_{k+1}\), \(0\le w\le H_{k+2}\).
We prove
\[
G_k\le U_k\le H_k\qquad(k\ge6).
\]
The base cases are exact:
\[
\begin{array}{c|cc}
k&U_k-G_k&H_k-U_k\\ \hline
6&\dfrac{92}{1859}&\dfrac{10}{11}\\[4pt]
7&\dfrac{19597}{163248}&\dfrac{11345}{15504}\\[4pt]
8&\dfrac{3424}{21315}&\dfrac{4272}{7105}.
\end{array}
\]
Assume \(G_j\le U_j\le H_j\) for \(j=k,k+1,k+2\), with \(k\ge6\).
Monotonicity of \(\Psi_k\) gives
\[
\Psi_k(G_k,G_{k+1},G_{k+2})
\le U_{k+3}\le
\Psi_k(H_k,H_{k+1},H_{k+2}).
\]
Verification. Direct expansion gives
\[
\Psi_k(G_k,G_{k+1},G_{k+2})-G_{k+3}
=\frac{R_-(k-6)}
{2(k+8)^2(k+10)^2(10k-11)(10k-1)(10k+9)(10k+139)},
\]
where
\[
\begin{aligned}
R_-(s)={}&3750000s^6+174475000s^5+3214281000s^4\\
&+30095015750s^3+153613752915s^2\\
&+416461589985s+482743428300,
\end{aligned}
\]
and
\[
H_{k+3}-\Psi_k(H_k,H_{k+1},H_{k+2})
=\frac{R_+(k-6)}
{k(k-2)(k-1)(k+8)^2(k+10)^2(k+13)},
\]
where
\[
R_+(s)=1900s^5+76662s^4+1135088s^3+7550562s^2
+22473804s+23833584.
\]
Every coefficient in \(R_-\) and \(R_+\) is positive, and the denominators
are positive for \(k\ge6\). This proves the induction.

Finally, for \(k\ge6\),
\[
U_{k+1}\le H_{k+1}=\frac{2880}{10k+110}
<\frac{2880}{10k+109}=G_k\le U_k.
\]
Thus \(\rho_{k+1}>\rho_k\) for \(k\ge6\), and the exact initial differences
prove Assumption SU5-RatioMono. \(\square\)

## Corollary 20 (SU(5) adjoint Q_3 closure)

For every \(n\ge0\),
\[
Q_3^{SU(5),{\rm adj}}(n)\ge0.
\]

Proof.
Theorem 19 discharges Assumption SU5-RatioMono, so Proposition 18 gives
strict positivity for every odd \(n\ge5\). The exact cases \(n=1,3\) are
listed in Estimate SU5-tail. Even \(n\) are pointwise nonnegative.
\(\square\)

Thus the SU(5) tail branch is closed. The remaining finite-rank SU(N)
branch is
\[
6\le N<n+2.
\]

## Proposition 21 (general SU(N) adjoint moment formula)

For \(N\ge2\), define
\[
E_k^{(N)}:=E_{SU(N)}[|\operatorname{tr}U|^{2k}],\qquad
m_k^{(N)}:=E_{SU(N)}[\chi_{\rm adj}(U)^k].
\]
Then
\[
E_k^{(N)}
=\sum_{\lambda\vdash k,\ \ell(\lambda)\le N}
(\dim\operatorname{Sp}_\lambda)^2,
\]
and
\[
m_k^{(N)}
=\sum_{a=0}^k(-1)^{k-a}\binom ka E_a^{(N)}.
\]
Equivalently, Gessel's determinant gives
\[
\sum_{k\ge0}E_k^{(N)}\frac{X^{2k}}{(k!)^2}
=\det(I_{|i-j|}(2X))_{1\le i,j\le N}.
\]
For every \(n\ge0\),
\[
\frac12Q_3^{SU(N),{\rm adj}}(n)
=\frac12\sum_{j=0}^{n}\binom nj
\left(m_{j+2}^{(N)}m_{n-j}^{(N)}
+m_j^{(N)}m_{n+2-j}^{(N)}
-2m_{j+1}^{(N)}m_{n+1-j}^{(N)}\right).
\]

Proof.
The function \(U\mapsto|\operatorname{tr}U|^{2k}\) is invariant under
scalar phases, so the \(SU(N)\) and \(U(N)\) Haar integrals agree.
Schur-Weyl duality gives
\[
(\mathbb C^N)^{\otimes k}
\cong\bigoplus_{\lambda\vdash k,\ \ell(\lambda)\le N}
S_\lambda(\mathbb C^N)\otimes\operatorname{Sp}_\lambda ,
\]
and Schur orthogonality gives the displayed finite partition sum for
\(E_k^{(N)}\). Since
\[
\chi_{\rm adj}(U)=|\operatorname{tr}U|^2-1,
\]
the displayed binomial transform gives \(m_k^{(N)}\). Gessel, §7,
Theorem 16, identifies the Bessel determinant with the same partition
sum by RSK. Finally, expanding
\[
\frac12E[(X-Y)^2(X+Y)^n]
\]
for independent copies \(X,Y\) with moments \(m_k^{(N)}\) gives the
displayed \(Q_3/2\) formula. \(\square\)

## Assumption SUN-RatioMono

Fix \(N\ge4\), and let
\[
m_k^{(N)}:=E_{SU(N)}[\chi_{\rm adj}(U)^k],\qquad
\rho_k^{(N)}:=m_{k+1}^{(N)}/m_k^{(N)} .
\]
Assume
\[
\rho_3^{(N)}=\frac92,\qquad
\frac92\le \rho_k^{(N)}\le \rho_{k+1}^{(N)}
\quad\text{for every }k\ge3.
\]

## Proposition 22 (universal SU(N) boundary reduction)

Under Assumption SUN-RatioMono, for the fixed \(N\ge4\),
\[
Q_3^{SU(N),{\rm adj}}(n)>0
\]
for every odd \(n\ge5\).

Proof.
Write \(m_k=m_k^{(N)}\) and \(\rho_k=\rho_k^{(N)}\). For \(0\le j\le n\),
put
\[
B_j=m_{j+2}m_{n-j}+m_jm_{n+2-j}-2m_{j+1}m_{n+1-j}.
\]
Then
\[
\frac12Q_3(n)=\frac12\sum_{j=0}^{n}\binom njB_j .
\]
For \(3\le j\le n-3\), Assumption SUN-RatioMono gives
\[
\frac{B_j}{m_{j+1}m_{n+1-j}}
=\frac{\rho_{j+1}}{\rho_{n-j}}
 +\frac{\rho_{n+1-j}}{\rho_j}-2\ge0
\]
by AM-GM, since
\[
\frac{\rho_{j+1}\rho_{n+1-j}}{\rho_{n-j}\rho_j}\ge1.
\]

It remains to check the depth-2 boundary contribution
\[
L_2:=B_0+nB_1+\binom n2B_2 .
\]
For \(N\ge4\), the stable initial moments give
\[
m_0=1,\qquad m_1=0,\qquad m_2=1,\qquad m_3=2,\qquad m_4=9.
\]
Therefore
\[
B_0=m_n+m_{n+2},\quad
B_1=2m_{n-1}-2m_n,\quad
B_2=9m_{n-2}+m_n-4m_{n-1}.
\]
Set
\[
M=m_{n-2},\qquad R=\rho_{n-2},\qquad S=\rho_{n-1},
\qquad T=\rho_n,\qquad U=\rho_{n+1}.
\]
Then
\[
m_{n-1}=RM,\quad m_n=SRM,\quad m_{n+2}=UTSRM.
\]
By Assumption SUN-RatioMono, \(U\ge T\), so
\[
\frac{L_2}{M}\ge
RS\left(T^2+\frac{n^2-5n+2}{2}\right)-2n(n-2)R
\,+\frac92n(n-1).
\]
Verification. Expanding \(L_2/M\), substituting the displayed ratio
relations, and replacing \(U\) by \(T\) gives exactly the displayed lower
bound.

For \(n\ge5\), \(R,S,T\ge9/2\), and
\[
\frac92\left(\frac{81}{4}+\frac{n^2-5n+2}{2}\right)-2n(n-2)
=\frac{2n^2-58n+765}{8}>0,
\]
because the quadratic has negative discriminant and positive leading
coefficient. Hence
\[
\frac{L_2}{M}
\ge
\frac92\cdot\frac{2n^2-58n+765}{8}
\,+\frac92n(n-1)
=\frac{9(10n^2-66n+765)}{16}>0.
\]
Thus \(L_2>0\). The interior terms are nonnegative, so
\(Q_3(n)/2\ge L_2>0\). \(\square\)

The next proposition discharges the stable-index part of
Assumption SUN-RatioMono.

## Proposition 23 (stable part of SUN-RatioMono)

For \(N\ge4\), the inequalities
\[
\rho_3^{(N)}=\frac92,\qquad
\frac92\le\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\]
hold for every \(3\le k\le N-2\).

Proof.
For \(r\le N\), Lemma 3 gives
\[
m_r^{(N)}=!r.
\]
Thus \(\rho_3^{(N)}=!4/!3=9/2\). It remains to check that the derangement
ratios are increasing. Let \(D_r:=!r\). The recurrence
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1}
\]
gives
\[
\frac{D_{r+2}}{D_{r+1}}-\frac{D_{r+1}}{D_r}
=1+\frac{(-1)^{r+2}}{D_{r+1}}-\frac{(-1)^{r+1}}{D_r}.
\]
If \(r\) is even, the right side is positive. If \(r\) is odd and
\(r\ge3\), the recurrence \(D_{r+1}=r(D_r+D_{r-1})\) shows
\(D_r\ge D_3\) and \(D_{r+1}\ge D_4\), so the right side is at least
\[
1-\frac1{D_4}-\frac1{D_3}=\frac7{18}>0.
\]
Therefore \(D_{r+2}/D_{r+1}\ge D_{r+1}/D_r\) for every \(r\ge3\).
For \(3\le k\le N-2\), both ratios \(\rho_k^{(N)}\) and
\(\rho_{k+1}^{(N)}\) are derangement ratios, so the claim follows.
\(\square\)

Consequently the remaining finite-rank SU(N) branch is reduced to the
transition range \(k\ge N-1\).

## Proposition 24 (first transition ratios)

**Superseded proof route.** Propositions 24--62 are no longer invoked by the
paper's theorem assembly.  The uniform transition-strip lemma proves all of
their conclusions from one derangement-error estimate, with exact GMP replay
only for ranks 6 through 18.

For every \(N\ge6\),
\[
\rho_{N-1}^{(N)}\le\rho_N^{(N)}
\qquad\text{and}\qquad
\rho_N^{(N)}\le\rho_{N+1}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). Proposition 21 and the first two nonstable RSK omissions
give
\[
m_{N-1}^{(N)}=D_{N-1},\qquad
m_N^{(N)}=D_N,\qquad
m_{N+1}^{(N)}=D_{N+1}-1,
\]
because the only partition of \(N+1\) with length \(>N\) is
\((1^{N+1})\), whose Specht dimension is \(1\), so the binomial transform
subtracts \(1\) from \(D_{N+1}\). Similarly,
\[
m_{N+2}^{(N)}=D_{N+2}-N(N+1),
\]
because the omitted partitions of \(N+2\) are \((1^{N+2})\) and
\((2,1^N)\), with Specht dimensions \(1\) and \(N+1\), and the binomial
correction is
\[
(N+2)-\bigl((N+1)^2+1\bigr)=-N(N+1).
\]

For the first inequality, the numerator of
\[
\frac{m_{N+1}^{(N)}}{m_N^{(N)}}-\frac{m_N^{(N)}}{m_{N-1}^{(N)}}
\]
is
\[
D_{N-1}(D_{N+1}-1)-D_N^2.
\]
Using
\[
\frac{D_{r+1}}{D_r}=r+1+\frac{(-1)^{r+1}}{D_r},
\]
this numerator equals
\[
\begin{cases}
D_N(D_{N-1}+1),&N\ \text{odd},\\
D_ND_{N-1}-2D_{N-1}-D_N,&N\ \text{even}.
\end{cases}
\]
The odd case is positive. In the even case, \(N\ge6\) gives
\(D_N\ge D_6=265\) and \(D_{N-1}\ge D_5=44\), so the numerator is positive.

For the second inequality, the numerator of
\[
\frac{m_{N+2}^{(N)}}{m_{N+1}^{(N)}}-\frac{m_{N+1}^{(N)}}{m_N^{(N)}}
\]
is
\[
D_N(D_{N+2}-N(N+1))-(D_{N+1}-1)^2.
\]
The same ratio identity gives this numerator as
\[
\begin{cases}
D_N(D_{N+1}-N(N+1)-1)+D_{N+1}-1,&N\ \text{odd},\\
D_N(D_{N+1}-N(N+1)+1)+3D_{N+1}-1,&N\ \text{even}.
\end{cases}
\]
For \(N\ge6\), \(D_{N+1}>N(N+1)+1\): this follows by induction from
\(D_7=1854\) and \(D_{r+1}=r(D_r+D_{r-1})\). Hence both displayed
expressions are positive. \(\square\)

## Proposition 25 (third transition ratio)

For every \(N\ge6\),
\[
\rho_{N+1}^{(N)}\le\rho_{N+2}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). Proposition 21 and the omitted partitions of
\(N+3\) give
\[
m_{N+3}^{(N)}
=D_{N+3}-\frac{(N+2)(N^3+2N^2+3)}{2}.
\]
Indeed, the omitted partitions of \(N+3\) are
\[
(1^{N+3}),\quad (2,1^{N+1}),\quad (3,1^N),\quad (2,2,1^{N-1}),
\]
with Specht dimensions
\[
1,\qquad N+2,\qquad \binom{N+2}{2},\qquad \frac{N(N+3)}2.
\]
Combining these omissions with the binomial corrections from
\(m_{N+1}^{(N)}\) and \(m_{N+2}^{(N)}\) gives the displayed formula.

Set \(A:=D_N\). The inequality is equivalent to positivity of
\[
\Delta:=m_{N+3}^{(N)}m_{N+1}^{(N)}-(m_{N+2}^{(N)})^2.
\]
Verification. Substituting the displayed formulas for
\(m_{N+1}^{(N)},m_{N+2}^{(N)},m_{N+3}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives the following two cases. If \(N\) is even, then
\[
\begin{aligned}
2\Delta={}&2A^2(N+1)^2(N+2)\\
&-A(N+1)(N+2)(N^3-2N^2-2N+15)+14N+26.
\end{aligned}
\]
If \(N\) is odd, then
\[
\begin{aligned}
\frac{2\Delta}{N+1}={}&2A^2(N+1)(N+2)\\
&-A(N-1)(N+2)(N^2-N-3)-2N^3+2N^2+2N-2.
\end{aligned}
\]
For \(N\ge6\), \(A=D_N\ge N^2\): the base case is \(D_6=265\ge36\), and
\(D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^2\) for \(r\ge6\) once
\(D_r\ge r^2\). In both displayed quadratics in \(A\), the derivative with
respect to \(A\) is positive for \(A\ge N^2\). At \(A=N^2\), the even
derivative is
\[
(N+1)(N+2)(3N^3+6N^2+2N-15)>0,
\]
and the odd derivative is
\[
(N+2)(3N^3+6N^2+2N-3)>0.
\]
Thus it suffices to substitute \(A=N^2\). In the even case this gives
\[
N^7+7N^6+16N^5-N^4-41N^3-30N^2+14N+26,
\]
and in the odd case this gives
\[
N^6+6N^5+10N^4-N^3-4N^2+2N-2.
\]
After substituting \(N=s+6\), both polynomials have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 26 (fourth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+2}^{(N)}\le\rho_{N+3}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). The omitted-partition hook calculation gives
\[
m_{N+4}^{(N)}
=D_{N+4}
-\frac{(N+3)(N^5+6N^4+10N^3+9N^2+16N+24)}{6}.
\]
Set
\[
\Delta:=m_{N+4}^{(N)}m_{N+2}^{(N)}-(m_{N+3}^{(N)})^2.
\]
Verification. Substituting the formulas for
\(m_{N+2}^{(N)},m_{N+3}^{(N)},m_{N+4}^{(N)}\), writing \(A:=D_N\), and
using \(D_{r+1}=(r+1)D_r+(-1)^{r+1}\), gives a quadratic in \(A\) in each
parity of \(N\). For \(N\ge6\), \(A=D_N\ge N^3\): the base case is
\(D_6=265\ge216\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^3
\]
for \(r\ge6\) once \(D_r\ge r^3\).

At \(A=N^3\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)
(11N^5+36N^4+32N^3-15N^2-34N+6)}{6}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)
(11N^5+36N^4+32N^3-15N^2-10N+18)}{6}>0.
\]
The two quintic factors are positive for \(N\ge6\) after substituting
\(N=s+6\), where their coefficient lists are positive.
Thus it suffices to substitute \(A=N^3\). In the even case this gives
\[
\frac{
10N^{11}+96N^{10}+366N^9+665N^8+406N^7-476N^6
-846N^5-314N^4+178N^3+407N^2+108N-264}{12},
\]
and in the odd case this gives
\[
\frac{
10N^{11}+96N^{10}+366N^9+665N^8+450N^7-180N^6
-154N^5+330N^4+162N^3-77N^2+36N+24}{12}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 27 (fifth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+3}^{(N)}\le\rho_{N+4}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). The omitted-partition hook calculation gives
\[
m_{N+5}^{(N)}
=D_{N+5}
-\frac{(N+4)(N^7+12N^6+52N^5+102N^4+112N^3+204N^2+459N+450)}{24}.
\]
Set
\[
\Delta:=m_{N+5}^{(N)}m_{N+3}^{(N)}-(m_{N+4}^{(N)})^2.
\]
For \(N=6\) and \(N=7\), direct substitution of the hook-length formulas gives
\[
\Delta=155283349651,\qquad \Delta=18890757347679,
\]
respectively.

Now assume \(N\ge8\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+3}^{(N)},m_{N+4}^{(N)},m_{N+5}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^4\):
the base case is \(D_8=14833\ge8^4\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^4
\]
for \(r\ge8\) once \(D_r\ge r^4\).

At \(A=N^4\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)
(47N^7+284N^6+536N^5+302N^4-88N^3-184N^2-327N-450)}{24}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)
(47N^7+284N^6+536N^5+302N^4-88N^3-88N^2+57N-18)}{24}>0.
\]
The two degree-7 factors are positive for \(N\ge8\), after substituting
\(N=s+8\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^4\). In the even case this gives
\[
\frac{
138N^{15}+2220N^{14}+14862N^{13}+53567N^{12}+111372N^{11}
+128518N^{10}+55161N^9-63990N^8-162738N^7-218765N^6
-173478N^5-28378N^4+66303N^3+80388N^2+82692N+47232}{144},
\]
and in the odd case this gives
\[
\frac{
138N^{15}+2220N^{14}+14862N^{13}+53567N^{12}+111372N^{11}
+129082N^{10}+63081N^9-18726N^8-27594N^7+2323N^6
+5514N^5-814N^4+135N^3-1980N^2-3708N-1152}{144}.
\]
After substituting \(N=s+8\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 28 (sixth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+4}^{(N)}\le\rho_{N+5}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). The omitted-partition hook calculation gives
\[
m_{N+6}^{(N)}
=D_{N+6}
-\frac{(N+5)(N^9+20N^8+160N^7+650N^6+1408N^5+1832N^4
+3495N^3+10170N^2+16776N+11808)}{120}.
\]
Set
\[
\Delta:=m_{N+6}^{(N)}m_{N+4}^{(N)}-(m_{N+5}^{(N)})^2.
\]
For \(N=6,7,8\), direct substitution of the hook-length formulas gives
\[
\Delta=15741914871746,\quad
\Delta=2425014101600951,\quad
\Delta=397297400564084400,
\]
respectively.

Now assume \(N\ge9\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+4}^{(N)},m_{N+5}^{(N)},m_{N+6}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^5\):
the base case is \(D_9=133496\ge9^5\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^5
\]
for \(r\ge9\) once \(D_r\ge r^5\).

At \(A=N^5\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)
(239N^9+2390N^8+8380N^7+12050N^6+5812N^5-772N^4
-2755N^3-5940N^2-8796N-5928)}{120}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)
(239N^9+2390N^8+8380N^7+12050N^6+5812N^5-772N^4
-2275N^3-2100N^2+804N+1032)}{120}>0.
\]
The two degree-9 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^5\). In the even case this gives
\[
\frac{
2856N^{19}+71400N^{18}+771480N^{17}+4720199N^{16}
+17984088N^{15}+44118924N^{14}+69221160N^{13}+64946536N^{12}
+24582372N^{11}-26861198N^{10}-72325716N^9-106255284N^8
-105581136N^7-58513648N^6-1006548N^5+30274179N^4
+41167860N^3+40503372N^2+25978464N+7971840}{2880},
\]
and in the odd case this gives
\[
\frac{
2856N^{19}+71400N^{18}+771480N^{17}+4720199N^{16}
+17984088N^{15}+44118924N^{14}+69232632N^{13}+65210392N^{12}
+27164292N^{11}-12856958N^{10}-26426340N^9-13589172N^8
+5949696N^7+11419232N^6+5099244N^5-498669N^4
-3028092N^3-2740308N^2-1237536N-322560}{2880}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 29 (seventh transition ratio)

For every \(N\ge6\),
\[
\rho_{N+5}^{(N)}\le\rho_{N+6}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). The omitted-partition hook calculation gives
\[
m_{N+7}^{(N)}
=D_{N+7}
-\frac{(N+6)(N^{11}+30N^{10}+380N^9+2625N^8+10608N^7
+25152N^6+37745N^5+74820N^4+255976N^3+595893N^2
+746490N+432600)}{720}.
\]
Set
\[
\Delta:=m_{N+7}^{(N)}m_{N+5}^{(N)}-(m_{N+6}^{(N)})^2.
\]
For \(N=6,7,8,9\), direct substitution of the hook-length formulas gives
\[
\Delta=1867810044306561,\quad
\Delta=363594262731677942,\quad
\Delta=71436857659292750544,\quad
\Delta=15332127989001241084746,
\]
respectively.

Now assume \(N\ge10\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+5}^{(N)},m_{N+6}^{(N)},m_{N+7}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^6\):
the base case is \(D_{10}=1334961\ge10^6\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^6
\]
for \(r\ge10\) once \(D_r\ge r^6\).

At \(A=N^6\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)
(1439N^{11}+21582N^{10}+122290N^9+323805N^8+394992N^7
+173244N^6-11981N^5-53160N^4-118126N^3-207111N^2
-297894N-212280)}{720}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)
(1439N^{11}+21582N^{10}+122290N^9+323805N^8+394992N^7
+173244N^6-11981N^5-50280N^4-80686N^3-34311N^2
+33306N+8040)}{720}>0.
\]
The two degree-11 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^6\). In the even case this gives
\[
\frac{
86280N^{23}+3105720N^{22}+49600440N^{21}+463201199N^{20}
+2805455120N^{19}+11548091910N^{18}+32890230195N^{17}
+64620971694N^{16}+85319021016N^{15}+70260689990N^{14}
+25408017080N^{13}-20918551116N^{12}-60315081378N^{11}
-100556714004N^{10}-130345248960N^9-120437174881N^8
-64612274512N^7-3681013662N^6+26893258635N^5
+39008961654N^4+47028359124N^3+46516837896N^2
+32413325760N+11113804800}{86400},
\]
and in the odd case this gives
\[
\frac{
86280N^{23}+3105720N^{22}+49600440N^{21}+463201199N^{20}
+2805455120N^{19}+11548091910N^{18}+32890230195N^{17}
+64621317054N^{16}+85330762536N^{15}+70436113190N^{14}
+26922214520N^{13}-12571890876N^{12}-29678109858N^{11}
-24899726964N^{10}-6606693840N^9+6988715039N^8
+7205818928N^7+3423176178N^6+107436075N^5-3543660666N^4
-4548145356N^3-2779606584N^2-1002536640N-154483200}{86400}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 30 (general transition correction formula)

For \(N\ge1\) and \(j\ge0\), define
\[
\Omega_t(N):=
\sum_{r=1}^t\ \sum_{\substack{\nu\vdash t-r\\ \nu_1\le N+r}}
\bigl(f^{(N+r,\nu)}\bigr)^2
\qquad(t\ge1),\qquad \Omega_0(N):=0,
\]
where \(f^\mu\) denotes the Specht dimension of the partition \(\mu\). Then
\[
m_{N+j}^{(N)}
=D_{N+j}
-\sum_{t=1}^j(-1)^{j-t}\binom{N+j}{N+t}\Omega_t(N).
\]
For fixed \(j\), if \(N\ge j\), the constraint \(\nu_1\le N+r\) is automatic.

Proof.
By Proposition 21,
\[
m_{N+j}^{(N)}
=\sum_{a=0}^{N+j}(-1)^{N+j-a}\binom{N+j}{a}E_a^{(N)}.
\]
For \(a\le N\), \(E_a^{(N)}=a!\). Hence the terms with \(a\le N\) and the
full symmetric-group sums for \(a>N\) combine to \(D_{N+j}\), and the only
correction comes from omitted partitions of \(a=N+t\), \(1\le t\le j\).
By RSK,
\[
(N+t)!-E_{N+t}^{(N)}
=\sum_{\substack{\lambda\vdash N+t\\ \ell(\lambda)>N}}(f^\lambda)^2.
\]
Conjugating \(\lambda\) preserves Specht dimension and changes the condition
\(\ell(\lambda)>N\) into \(\lambda'_1>N\). Writing
\(\lambda'=(N+r,\nu)\) gives \(1\le r\le t\), \(\nu\vdash t-r\), and
\(\nu_1\le N+r\), and this parametrization is bijective. Thus
\[
(N+t)!-E_{N+t}^{(N)}=\Omega_t(N).
\]
Substituting this identity into the binomial transform gives the displayed
formula. If \(N\ge j\), then \(\nu_1\le|\nu|=t-r\le j-1\le N\le N+r\), so the
constraint is automatic. \(\square\)

## Proposition 31 (eighth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+6}^{(N)}\le\rho_{N+7}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). Proposition 30 with \(j=8\) gives
\[
m_{N+8}^{(N)}
=D_{N+8}
-\frac{(N+7)(N^{13}+42N^{12}+770N^{11}+8043N^{10}+52248N^9
+215628N^8+557429N^7+945144N^6+1926442N^5+7350735N^4
+22071126N^3+39473448N^2+41032944N+21150720)}{5040}.
\]
Set
\[
\Delta:=m_{N+8}^{(N)}m_{N+6}^{(N)}-(m_{N+7}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11\), direct substitution of the hook-length formulas gives
\[
\begin{gathered}
\Delta=254881456479831585,\quad
\Delta=62735747654515502252,\quad
\Delta=14737022237545368612140,\\
\Delta=3661918679655379068792732,\quad
\Delta=1004891295943303392780074729,\\
\Delta=308070944933634724735782382320,
\end{gathered}
\]
respectively.

Now assume \(N\ge12\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+6}^{(N)},m_{N+7}^{(N)},m_{N+8}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^7\):
the base case is \(D_{12}=176214841\ge12^7\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^7
\]
for \(r\ge12\) once \(D_r\ge r^7\).

At \(A=N^7\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)
(10079N^{13}+211652N^{12}+1763692N^{11}+7407211N^{10}
+16367022N^9+17786244N^8+7260535N^7-241658N^6
-1230964N^5-3086881N^4-5746692N^3-10496808N^2
-14623512N-9231120)}{5040}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)
(10079N^{13}+211652N^{12}+1763692N^{11}+7407211N^{10}
+16367022N^9+17786244N^8+7260535N^7-241658N^6
-1210804N^5-2703841N^4-2964612N^3-920808N^2
+819048N-68400)}{5040}>0.
\]
The two degree-13 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^7\). In the even case this gives
\[
\frac{
3628080N^{27}+177770880N^{26}+3936229920N^{25}+52083425519N^{24}
+458831595540N^{23}+2837396488142N^{22}+12649241058828N^{21}
+41078558241680N^{20}+96956564398746N^{19}+163797560671127N^{18}
+191537655117456N^{17}+143974737619232N^{16}+52258478461056N^{15}
-26615657296432N^{14}-88019858899116N^{13}-161710593010423N^{12}
-248554381128132N^{11}-301845709496446N^{10}-259447058738268N^9
-130097528943160N^8-11300305552998N^7+41139115210065N^6
+64686077172600N^5+90021632989752N^4+110361100043808N^3
+107698264845264N^2+69477664913280N+21183410304000}{3628800},
\]
and in the odd case this gives
\[
\frac{
3628080N^{27}+177770880N^{26}+3936229920N^{25}+52083425519N^{24}
+458831595540N^{23}+2837396488142N^{22}+12649241058828N^{21}
+41078558241680N^{20}+96956578912506N^{19}+163798242807767N^{18}
+191552052313776N^{17}+144154946036192N^{16}+53746198276896N^{15}
-18095672650672N^{14}-53377258997196N^{13}-61157886518743N^{12}
-42284759110212N^{11}-11123170298206N^{10}+6189073739652N^9
+6942335995400N^8+4815792092442N^7+1217599240305N^6
-5537705272680N^5-9503247564648N^4-7717958864352N^3
-3640231517616N^2-1072418140800N-114916838400}{3628800}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 32 (ninth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+7}^{(N)}\le\rho_{N+8}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). Proposition 30 with \(j=9\) gives
\[
m_{N+9}^{(N)}
=D_{N+9}
-\frac{(N+8)(N^{15}+56N^{14}+1400N^{13}+20524N^{12}
+194544N^{11}+1236984N^{10}+5284802N^9+14842132N^8
+28053592N^7+58395428N^6+238401576N^5+876264984N^4
+2046317673N^3+3017236932N^2+2792508012N+1350412560)}{40320}.
\]
Set
\[
\Delta:=m_{N+9}^{(N)}m_{N+7}^{(N)}-(m_{N+8}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12\), direct substitution of the hook-length formulas
gives
\[
\begin{gathered}
\Delta=39406224712279244637,\quad
\Delta=12298164094652542220156,\\
\Delta=3450352965360015965435625,\quad
\Delta=988760159115858016357960584,\\
\Delta=306949458230346825740138150231,\quad
\Delta=105322353384151727669764906060824,\\
\Delta=40047731347388887220476970687998229,
\end{gathered}
\]
respectively.

Now assume \(N\ge13\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+7}^{(N)},m_{N+8}^{(N)},m_{N+9}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^8\):
the base case is \(D_{13}=2290792932\ge13^8\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^8
\]
for \(r\ge13\) once \(D_r\ge r^8\).

At \(A=N^8\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)
(80639N^{15}+2257880N^{14}+25965408N^{13}+158048380N^{12}
+545823040N^{11}+1058907360N^{10}+1053872374N^9+406395436N^8
-5819232N^7-33752236N^6-98449288N^5-196343616N^4
-375403617N^3-722839524N^2-948259404N-546673680)}{40320}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)
(80639N^{15}+2257880N^{14}+25965408N^{13}+158048380N^{12}
+545823040N^{11}+1058907360N^{10}+1053872374N^9+406395436N^8
-5819232N^7-33590956N^6-94256008N^5-152636736N^4
-142999137N^3-63043044N^2-10093644N-32754960)}{40320}>0.
\]
The two degree-15 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^8\). In the even case this gives
\[
\frac{
203207760N^{31}+13005236160N^{30}+381213816480N^{29}
+6782145209279N^{28}+81808950983676N^{27}+707833979412112N^{26}
+4531918414599555N^{25}+21829872114747802N^{24}
+79624432197181536N^{23}+219492869261645473N^{22}
+452529177387326514N^{21}+683416238546542144N^{20}
+728751829374595089N^{19}+510942615709155784N^{18}
+187599979046304804N^{17}-48360938835968219N^{16}
-206711163460260120N^{15}-410345447951221652N^{14}
-716651012383157271N^{13}-1060168478665942646N^{12}
-1222432699215532140N^{11}-986837994923591317N^{10}
-470689353671021166N^9-58143601846967556N^8
+101547737679478083N^7+173228684216384412N^6
+271504706578245504N^5+391246206296574096N^4
+482897364810579216N^3+449865637812246528N^2
+267845260488986880N+74654636674867200}{203212800},
\]
and in the odd case this gives
\[
\frac{
203207760N^{31}+13005236160N^{30}+381213816480N^{29}
+6782145209279N^{28}+81808950983676N^{27}+707833979412112N^{26}
+4531918414599555N^{25}+21829872114747802N^{24}
+79624432197181536N^{23}+219492870074486593N^{22}
+452529227783354994N^{21}+683417663449384864N^{20}
+728776156677658929N^{19}+511222487794714024N^{18}
+189892224270589764N^{17}-34596740899566299N^{16}
-145252925474280600N^{15}-205660757273597012N^{14}
-211681808498615511N^{13}-153902645048552246N^{12}
-78412551144703020N^{11}-31308523283764597N^{10}
-8203193677245486N^9+7375763836463964N^8
+6774568413863523N^7-9025723017521028N^6
-22109403840574176N^5-22153913264092464N^4
-12445576332334704N^3-3843300340829952N^2
-494571119051520N+209354703667200}{203212800}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 33 (tenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+8}^{(N)}\le\rho_{N+9}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge10\), Proposition 30 with \(j=10\) gives
\[
m_{N+10}^{(N)}
=D_{N+10}
-\frac{(N+9)(N^{17}+72N^{16}+2352N^{15}+45972N^{14}
+596400N^{13}+5376672N^{12}+34204714N^{11}+152398260N^{10}
+463471296N^9+968105916N^8+2052382008N^7+8648198496N^6
+37534002913N^5+109023460740N^4+208533189276N^3
+270696395952N^2+235518865920N+110402611200)}{362880}.
\]
Set
\[
\Delta:=m_{N+10}^{(N)}m_{N+8}^{(N)}-(m_{N+9}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13\), direct substitution of the hook-length
formulas gives
\[
\begin{gathered}
\Delta=6813854089176842631054,\quad
\Delta=2708495048576903530387324,\\
\Delta=908027814739482279470762685,\quad
\Delta=299384526399163870755924652400,\\
\Delta=104686217257849609761817491133909,\quad
\Delta=39998354664073776381274210861634460,\\
\Delta=16818179640707502071612872009415934596,\quad
\Delta=7771529846902770914762363686452853966096,
\end{gathered}
\]
respectively.

Now assume \(N\ge14\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+8}^{(N)},m_{N+9}^{(N)},m_{N+10}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^9\):
the base case is \(D_{14}=32071101049\ge14^9\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^9
\]
for \(r\ge14\) once \(D_r\ge r^9\).

At \(A=N^9\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)
(725759N^{17}+26127306N^{16}+396263688N^{15}
+3292030404N^{14}+16292448984N^{13}+48831383160N^{12}
+85728262286N^{11}+79532524224N^{10}+29259948288N^9
-162982836N^8-1067275056N^7-3623448024N^6-8074733689N^5
-15531703650N^4-33397172388N^3-64616905944N^2
-77871043872N-41039913600)}{362880}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)
(725759N^{17}+26127306N^{16}+396263688N^{15}
+3292030404N^{14}+16292448984N^{13}+48831383160N^{12}
+85728262286N^{11}+79532524224N^{10}+29259948288N^9
-162982836N^8-1065823536N^7-3574096344N^6-7379455609N^5
-10293167970N^4-10725881508N^3-8668067544N^2
-5698566432N-4045749120)}{362880}>0.
\]
The two degree-17 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^9\). In the even case this gives
\[
\frac{
14631281280N^{35}+1185133057920N^{34}+44420508034560N^{33}
+1022372961281279N^{32}+16168413464622608N^{31}
+186264906475449496N^{30}+1616864226783938128N^{29}
+10788828244891243728N^{28}+55949597276961728760N^{27}
+226444556605714785636N^{26}+714063120728743907736N^{25}
+1741096373096470982472N^{24}+3235453092149751612528N^{23}
+4474734641068840143408N^{22}+4432370376224126750280N^{21}
+2935568247887338081698N^{20}+1081335774148378761960N^{19}
-91342115186841387264N^{18}-734777294170600610736N^{17}
-1559168806445732162160N^{16}-3011902406495137547544N^{15}
-5175015858974751915116N^{14}-7379221267437810603032N^{13}
-8030880144631698328216N^{12}-6075151132794019946704N^{11}
-2752994147208653644368N^{10}-410619060599688867816N^9
+383286964954326542679N^8+710833705851923145864N^7
+1235972266511992376472N^6+2017990649694284772480N^5
+2943169600358170980720N^4+3492876402508716978048N^3
+3018460839226125878016N^2+1645476259521404620800N
+420509959346818252800}{14631321600},
\]
and in the odd case this gives
\[
\frac{
14631281280N^{35}+1185133057920N^{34}+44420508034560N^{33}
+1022372961281279N^{32}+16168413464622608N^{31}
+186264906475449496N^{30}+1616864226783938128N^{29}
+10788828244891243728N^{28}+55949597276961728760N^{27}
+226444556605714785636N^{26}+714063120787269113496N^{25}
+1741096377719960785992N^{24}+3235453260643698035568N^{23}
+4474738398031240289328N^{22}+4432427691428529858120N^{21}
+2936201979699552619938N^{20}+1086583493800202025960N^{19}
-58185803740848374784N^{18}-573456516535312353456N^{17}
-953658007305010078320N^{16}-1267427857864125984024N^{15}
-1365467663704600552556N^{14}-1211355638946195669272N^{13}
-893012536028196188056N^{12}-529415940950777408464N^{11}
-206265796375704736848N^{10}-1450110262309555176N^9
+46239024958444292439N^8+8495692177623849864N^7
-27206553304365113448N^6-28786618579761429120N^5
-77513304330579600N^4+25158908227716402048N^3
+25943640928385446656N^2+12763895810826424320N
+3533118047856230400}{14631321600}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 34 (eleventh transition ratio)

For every \(N\ge6\),
\[
\rho_{N+9}^{(N)}\le\rho_{N+10}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge11\), Proposition 30 with \(j=11\) gives
\[
m_{N+11}^{(N)}
=D_{N+11}
-\frac{(N+10)(N^{19}+90N^{18}+3720N^{17}+93375N^{16}
+1584816N^{15}+19147824N^{14}+168617230N^{13}+1087375680N^{12}
+5073075816N^{11}+16658275698N^{10}+38231053284N^9
+82688229360N^8+348168913021N^7+1733689473462N^6
+6072943640652N^5+14383340130555N^4+23914938342996N^3
+28877873970276N^2+24560356173264N+11321406391680)}{3628800}.
\]
Set
\[
\Delta:=m_{N+11}^{(N)}m_{N+9}^{(N)}-(m_{N+10}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15\), direct substitution of the hook-length
formulas gives
\[
\begin{gathered}
\Delta=1302995252460920344262919,\quad
\Delta=663571426044273782881718696,\\
\Delta=266305274891498475082793422359,\quad
\Delta=100911166262323382289653695410112,\\
\Delta=39617618166789524613364731947471087,\quad
\Delta=16783173678848999079599868332123462790,\\
\Delta=7768530068175868341187849118317526816914,\quad
\Delta=3932268596380642080901282525338082573116284,\\
\Delta=2170736484858344363193806147949629795366071823,\quad
\Delta=1302452468793726309607411028670900202696991799844,
\end{gathered}
\]
respectively.

Now assume \(N\ge16\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+9}^{(N)},m_{N+10}^{(N)},m_{N+11}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{10}\):
the base case is \(D_{16}=7697064251745\ge16^{10}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{10}
\]
for \(r\ge16\) once \(D_r\ge r^{10}\).

At \(A=N^{10}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)
(7257599N^{19}+326591930N^{18}+6314109810N^{17}
+68584279875N^{16}+459209653104N^{15}+1954649490456N^{14}
+5252162543210N^{13}+8510945490680N^{12}+7450478052204N^{11}
+2633511008442N^{10}-5204602524N^9-38165698200N^8
-149585030821N^7-381053235202N^6-764556288642N^5
-1648317189825N^4-3882625480836N^3-7113947828076N^2
-7817512952304N-3808984798080)}{3628800}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)
(7257599N^{19}+326591930N^{18}+6314109810N^{17}
+68584279875N^{16}+459209653104N^{15}+1954649490456N^{14}
+5252162543210N^{13}+8510945490680N^{12}+7450478052204N^{11}
+2633511008442N^{10}-5204602524N^9-38151183000N^8
-148960877221N^7-369658803202N^6-649595904642N^5
-950136069825N^4-1282720917636N^3-1356827597676N^2
-952041080304N-479422903680)}{3628800}>0.
\]
The two degree-19 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^{10}\). In the even case this gives
\[
\frac{
1316818581120N^{39}+131681849040000N^{38}+6142957702963200N^{37}
+177572936212559999N^{36}+3563238414935969136N^{35}
+52693356770756301882N^{34}+594999088603658711811N^{33}
+5244622260023880344598N^{32}+36581328421669881488016N^{31}
+203445751280724833675692N^{30}+904690127051714135266416N^{29}
+3211998509164533260121144N^{28}+9052549741086548823216684N^{27}
+20038843038951799268232168N^{26}+34260849541749643515490320N^{25}
+44105444323472123596110234N^{24}+41099123431247300579571552N^{23}
+25925394196131242940119412N^{22}+9476445850738412520303810N^{21}
+234407091245411712671028N^{20}-3769682222055109172382792N^{19}
-8499148763869587087601700N^{18}-17850380992001169060094992N^{17}
-34524051409529375227722264N^{16}-57976936513196208273376308N^{15}
-78628304211796981334838264N^{14}-80291411192689171659636480N^{13}
-56974988015532638264881441N^{12}-24555802553568406005486768N^{11}
-4026448477128554472975246N^{10}+2127121135794357050523123N^9
+4284348928552733470564950N^8+8162272272575944895997384N^7
+14735586293829020959544016N^6+24363866045547124187103792N^5
+34394622649121195963590752N^4+38148137360984717765438976N^3
+30357071071567664155799040N^2+15237782678225940166656000N
+3612581898146881486848000}{1316818944000},
\]
and in the odd case this gives
\[
\frac{
1316818581120N^{39}+131681849040000N^{38}+6142957702963200N^{37}
+177572936212559999N^{36}+3563238414935969136N^{35}
+52693356770756301882N^{34}+594999088603658711811N^{33}
+5244622260023880344598N^{32}+36581328421669881488016N^{31}
+203445751280724833675692N^{30}+904690127051714135266416N^{29}
+3211998509169800535171384N^{28}+9052549741602741759996204N^{27}
+20038843062496516823169768N^{26}+34260850205452541139163920N^{25}
+44105457271079145472901274N^{24}+41099308908057065758527072N^{23}
+25927414432717452619257972N^{22}+9493536821293321168643010N^{21}
+348130376849852572319028N^{20}-3170757880427662853418312N^{19}
-6000006159998510332123940N^{18}-9619763312238077789221392N^{17}
-13324144840214184478544664N^{16}-15929279981260132448965428N^{15}
-15959188414377164245713144N^{14}-12737371153893825642416640N^{13}
-7550836119681809728165921N^{12}-2866709665686902165370288N^{11}
-283386411904110788417166N^{10}+411941229139800384883443N^9
+438046134346735474883670N^8+541073075929732966649544N^7
+820339818922365544763856N^6+1165830436343236747882032N^5
+1300939053569522772289632N^4+1072872451521935912295936N^3
+616056416143951805038080N^2+220641358466597393203200N
+42952595075103522816000}{1316818944000}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 35 (twelfth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+10}^{(N)}\le\rho_{N+11}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge12\), Proposition 30 with \(j=12\) gives
\[
m_{N+12}^{(N)}
=D_{N+12}
-\frac{(N+11)(N^{21}+110N^{20}+5610N^{19}+175725N^{18}
+3773616N^{17}+58691424N^{16}+680317990N^{15}
+5949406760N^{14}+39211703916N^{13}+191979965958N^{12}
+678869564844N^{11}+1704866547480N^{10}+3783389804821N^9
+15470929796762N^8+86095094847102N^7+355867186762065N^6
+1014225525106596N^5+2064813076161396N^4
+3146510994392304N^3+3697149122140320N^2
+3140744998896000N+1423735367040000)}{39916800}.
\]
Set
\[
\Delta:=m_{N+12}^{(N)}m_{N+10}^{(N)}-(m_{N+11}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16\), direct substitution of the
hook-length formulas gives
\[
\begin{gathered}
\Delta=272877302308700356233769809,\quad
\Delta=179268984244846216713178911600,\\
\Delta=86367303292900676831439055272561,\quad
\Delta=37611774576211473420036594050576384,\\
\Delta=16542336378982215305089462870750721417,\\
\Delta=7742480765706091807652721778761118783311,\\
\Delta=3929668632220844281683139813639263059975000,\\
\Delta=2170492735985555040408851825992397140907456975,\\
\Delta=1302430714340459432830328053390292912495960817361,\\
\Delta=846593105592369170486931923524422038900376893954096,\\
\Delta=594309584979454468445716960726247172905985251300624100,
\end{gathered}
\]
respectively.

Now assume \(N\ge17\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+10}^{(N)},m_{N+11}^{(N)},m_{N+12}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{11}\):
the base case is \(D_{17}=130850092279664\ge17^{11}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{11}
\]
for \(r\ge17\) once \(D_r\ge r^{11}\).

At \(A=N^{11}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)(N+11)
(79833599N^{21}+4390847912N^{20}+105380348480N^{19}
+1448979755685N^{18}+12595585234914N^{17}+72014283419628N^{16}
+272785712700578N^{15}+671360116628300N^{14}
+1018162423442384N^{13}+848521793094294N^{12}
+289694483034552N^{11}-186564394992N^{10}-1521622553341N^9
-6817436987020N^8-19942140491548N^7-43546443920991N^6
-92594058307146N^5-233736024964464N^4-553160282329512N^3
-930156772252752N^2-937303362624960N-427734367200000)}{39916800}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)(N+11)
(79833599N^{21}+4390847912N^{20}+105380348480N^{19}
+1448979755685N^{18}+12595585234914N^{17}+72014283419628N^{16}
+272785712700578N^{15}+671360116628300N^{14}
+1018162423442384N^{13}+848521793094294N^{12}
+289694483034552N^{11}-186564394992N^{10}-1521462886141N^9
-6808974625420N^8-19748144843548N^7-41028492176991N^6
-72268423747146N^5-128337151569264N^4-203842457843112N^3
-221353356316752N^2-145476355700160N-61483596652800)}{39916800}>0.
\]
The two degree-21 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^{11}\). In the even case this gives
\[
\frac{
144850080211200N^{43}+17526859585804800N^{42}
+995844285565113600N^{41}+35316621309172895999N^{40}
+876479974044461740620N^{39}+16175993422351615760570N^{38}
+230286932289932088180900N^{37}+2588948687102732050241688N^{36}
+23338543426025790522694662N^{35}+170342114406930109024629625N^{34}
+1012111786918885016790744300N^{33}+4904860360782145744127823292N^{32}
+19362309629769577232072295936N^{31}+61981093508742751734342338016N^{30}
+159595181612592099375364290360N^{29}+326468006464724468574925015326N^{28}
+520981647609010190816244294792N^{27}+631582760523834845156800348644N^{26}
+558779868593183224933445177928N^{25}+337733266025304868874645321280N^{24}
+121458099190594581337056423300N^{23}+10763655061076407001864969050N^{22}
-26865260107757995530841353000N^{21}-64140900325567823420149691720N^{20}
-144884248379844569778161010720N^{19}
-308801202516785504896833380640N^{18}
-590809357958582684104457597640N^{17}
-950985641407182559837049472829N^{16}
-1214488399507928564322144771252N^{15}
-1163260160067421476281438937070N^{14}
-777799106420845184761434411660N^{13}
-319856412997763340268734348168N^{12}
-54663084987372436600829621466N^{11}
+16795141958845983269969720109N^{10}
+36560451254695638645613795740N^9+75462744797709934688760517308N^8
+148124335221662483614214145888N^7
+270325552124589192646930131696N^6
+434919142129449838284431681472N^5
+577842484413482400704730997824N^4
+593424463022873951093633687040N^3
+436455061465925598018488755200N^2
+203796016972510718954354688000N
+45401378120881785689702400000}{144850083840000},
\]
and in the odd case this gives
\[
\frac{
144850080211200N^{43}+17526859585804800N^{42}
+995844285565113600N^{41}+35316621309172895999N^{40}
+876479974044461740620N^{39}+16175993422351615760570N^{38}
+230286932289932088180900N^{37}+2588948687102732050241688N^{36}
+23338543426025790522694662N^{35}+170342114406930109024629625N^{34}
+1012111786918885016790744300N^{33}+4904860360782145744127823292N^{32}
+19362309629770156632400398336N^{31}+61981093508811700373147022816N^{30}
+159595181616438158722013250360N^{29}+326468006598366049320692464926N^{28}
+520981650851292635070102189192N^{27}+631582818866081258973600681444N^{26}
+558780675937031092682785107528N^{25}+337742057052753954434349641280N^{24}
+121534533809244320522629940100N^{23}+11299038738424123164646639450N^{22}
-23830784101376408890348309800N^{21}-50215509290406858678153774920N^{20}
-93299818374160698034251685920N^{19}
-155613370234443032663567550240N^{18}
-230265610108211461066833610440N^{17}
-290006990781385261403058893629N^{16}
-294674335670920091800625436852N^{15}
-229129025806408248057889650670N^{14}
-128483797910346747829903582860N^{13}
-46217687294159923243066533768N^{12}
-5909849570580122280232181466N^{11}
+4752306296289980491077713709N^{10}
+7659983809461214892601763740N^9+13016327620033616290065234108N^8
+21611168782603114824011304288N^7
+31481242946330084702987437296N^6
+37475002514942935916718324672N^5
+34987935865842496849840041024N^4
+24617718402167312617963875840N^3
+12169366880448632893505971200N^2
+3737597377565470139578368000N
+607263651091980546048000000}{144850083840000}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 36 (thirteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+11}^{(N)}\le\rho_{N+12}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge13\), Proposition 30 with \(j=13\) gives
\[
m_{N+13}^{(N)}
=D_{N+13}
-\frac{(N+12)(N^{23}+132N^{22}+8140N^{21}+311058N^{20}
+8236096N^{19}+159993636N^{18}+2353952293N^{17}
+26674834758N^{16}+234171846568N^{15}+1585762966260N^{14}
+8156090263000N^{13}+30982990092216N^{12}
+84880263053575N^{11}+194924932705536N^{10}
+756476882723056N^9+4579908770934882N^8
+21968738042464312N^7+73891428189425484N^6
+179394061557788451N^5+329706813191254686N^4
+482843373351470748N^3+569544744919538952N^2
+485178256571326560N+214656221754374400)}{479001600}.
\]
Set
\[
\Delta:=m_{N+13}^{(N)}m_{N+11}^{(N)}-(m_{N+12}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16,17\), direct substitution of the
hook-length formulas gives
\[
\begin{gathered}
\Delta=62053268258617514038409289215,\quad
\Delta=52987197630443447008549718041680,\\
\Delta=30759385305870137640156148668462791,\quad
\Delta=15407670359929730333204554782925602048,\\
\Delta=7581355526868485858726286698650683204175,\\
\Delta=3909293271892088334628716680283816520060365,\\
\Delta=2168137627574512972871926918228698223754463585,\\
\Delta=1302177092110927967359290836350877822641274743071,\\
\Delta=846567287756840444021060413294524282963701057160479,\\
\Delta=594307072710621917346273218221073015661936271940019944,\\
\Delta=449297920809269172910381633775090196454278525118772877750,\\
\Delta=364830091306238796022944398801054022000244264883978913595172,
\end{gathered}
\]
respectively.

Now assume \(N\ge18\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+11}^{(N)},m_{N+12}^{(N)},m_{N+13}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{12}\):
the base case is \(D_{18}=2355301661033953\ge18^{12}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{12}
\]
for \(r\ge18\) once \(D_r\ge r^{12}\).

At \(A=N^{12}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)(N+11)(N+12)
(958003199N^{23}+63228211092N^{22}+1844156154632N^{21}
+31297964381706N^{20}+342412374434504N^{19}
+2526788955986328N^{18}+12779316710871095N^{17}
+44064052732339794N^{16}+100837554265997648N^{15}
+144579845874675684N^{14}+115481323388802824N^{13}
+38240154606998640N^{12}-7416905850415N^{11}
-66880101881736N^{10}-339465606541996N^9-1135841364259206N^8
-2766420293208880N^7-5975386059006672N^6-15093199594333575N^5
-41367878623968534N^4-92036223166693356N^3
-141255468509982696N^2-132137517351494880N
-57237976099180800)}{479001600}>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)(N+2)(N+3)(N+4)(N+5)(N+6)(N+7)(N+8)(N+9)(N+10)(N+11)(N+12)
(958003199N^{23}+63228211092N^{22}+1844156154632N^{21}
+31297964381706N^{20}+342412374434504N^{19}
+2526788955986328N^{18}+12779316710871095N^{17}
+44064052732339794N^{16}+100837554265997648N^{15}
+144579845874675684N^{14}+115481323388802824N^{13}
+38240154606998640N^{12}-7416905850415N^{11}
-66878185875336N^{10}-339342982132396N^9-1132396384752006N^8
-2710597446744880N^7-5399108814078672N^6-11145429351671175N^5
-23263413441965334N^4-37420629343256556N^3
-38191493523333096N^2-23221342232851680N
-8892860974905600)}{479001600}>0.
\]
The two degree-23 factors are positive for \(N\ge6\), after substituting
\(N=s+6\) and checking that all coefficients are positive. Thus it suffices
to substitute \(A=N^{12}\). In the even case this gives
\[
\frac{
19120211026963200N^{47}+2753310386206195200N^{46}
+187186865685768748800N^{45}+7990565610251187532799N^{44}
+240268222736413888262180N^{43}+5411581514892142553756012N^{42}
+94779986498920996122987203N^{41}+1322743012604026441620949394N^{40}
+14953257659111806051799675432N^{39}+138450845556225660099864408439N^{38}
+1057408463526914093008661740042N^{37}+6688256180845698271701481430068N^{36}
+35082032042042365280929191015919N^{35}+152427850600119942862778339290972N^{34}
+546610459201586770341186854673892N^{33}+1607585005222803517728718605316854N^{32}
+3840719702871919184253441435037008N^{31}
+7353481070419569452016633706265136N^{30}
+11069020644132581891435572020258366N^{29}
+12746798832837278520874165499238540N^{28}
+10782588401559820926950746017888120N^{27}
+6274375519252645588079212732431022N^{26}
+2208707099230164544214042358427612N^{25}
+271905160499617001847074568301840N^{24}
-257743753535544124243777792680530N^{23}
-650312394283471176664985373884240N^{22}
-1569095665382405838597917263048640N^{21}
-3632682396665117919305459668733749N^{20}
-7735538165343960924474535326781732N^{19}
-14289363334315139001852459432006652N^{18}
-21748989500265634747608123428279089N^{17}
-26101322024773668352504772881229758N^{16}
-23555724231082670118215636063967904N^{15}
-14943259751926552335501607920590325N^{14}
-5893695902623408498778989130263494N^{13}
-1020098222785100192290690706862084N^{12}
+183082669467677016475161744927027N^{11}
+427423675500660341803353116845716N^{10}
+946242703607057100217759275086652N^9
+1995815786948078822251483055779440N^8
+3951047050850409073505076593232720N^7
+7045082579438352208508276532438720N^6
+10735421567475784671063090876430656N^5
+13279664562514136411965429347046656N^4
+12655604191737877099340775123778560N^3
+8684692127219177055479769822720000N^2
+3820059519788590965342664040448000N
+810357106786201371079994572800000}{19120211066880000},
\]
and in the odd case this gives
\[
\frac{
19120211026963200N^{47}+2753310386206195200N^{46}
+187186865685768748800N^{45}+7990565610251187532799N^{44}
+240268222736413888262180N^{43}+5411581514892142553756012N^{42}
+94779986498920996122987203N^{41}+1322743012604026441620949394N^{40}
+14953257659111806051799675432N^{39}+138450845556225660099864408439N^{38}
+1057408463526914093008661740042N^{37}+6688256180845698271701481430068N^{36}
+35082032042042365280929191015919N^{35}+152427850600120019343622526977372N^{34}
+546610459201597630621058153131492N^{33}+1607585005223530621113881243985654N^{32}
+3840719702902437870835721177283408N^{31}
+7353481071320300336613192514303536N^{30}
+11069020664005899517350497203592766N^{29}
+12746799173027271169372031766246540N^{28}
+10782593029626938023683270080736120N^{27}
+6274426365344730363114624831797422N^{26}
+2209163000285912190865319697077212N^{25}
+275262933419085734928011425722640N^{24}
-237364403972958666757924909608530N^{23}
-548353028571137020279118300041040N^{22}
-1149628399015253634496734141563840N^{21}
-2221204672031764431724443523133749N^{20}
-3884531082863188313691266938813732N^{19}
-5879092997314191578248334117113852N^{18}
-7316187407263874095457154211664689N^{17}
-7152458023866249803115462729524158N^{16}
-5270406590426180819721980717612704N^{15}
-2785347601597026220247215741991925N^{14}
-953764311411078640812474307889094N^{13}
-134349072420915108371801602042884N^{12}
+69820848278418242317791928786227N^{11}
+137922771831406432545227386708116N^{10}
+269657430839421123109886036481852N^9
+492826685856772211099724519062640N^8
+791932276977010728209213779267920N^7
+1065956258392594488799621337661120N^6
+1158409794067733799412749995432256N^5
+982932188230926113459633226061056N^4
+623610428713325490852551388994560N^3
+275113702952976811600638840729600N^2
+75063266080958654307515621376000N
+11081374031640033700550737920000}{19120211066880000}.
\]
After substituting \(N=s+6\), both numerators have only positive
coefficients. Hence \(\Delta>0\). \(\square\)

## Proposition 37 (fourteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+12}^{(N)}\le\rho_{N+13}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge14\), Proposition 30 with \(j=14\) gives
\[
m_{N+14}^{(N)}
=D_{N+14}
-\frac{(N+13)(N^{25}+156N^{24}+11440N^{23}+523614N^{22}
+16745872N^{21}+396866184N^{20}+7210326409N^{19}
+102425880414N^{18}+1149091574488N^{17}+10201023393132N^{16}
+71221302430408N^{15}+384943776507888N^{14}
+1567279683373519N^{13}+4671875619383736N^{12}
+11213209579641940N^{11}+40679629366987182N^{10}
+260076301677666760N^9+1427610032828721456N^8
+5587295295274851015N^7+15897860694126056646N^6
+34378050542667476580N^5+59744145385520856936N^4
+87430878993039573888N^3+104869944038560237056N^2
+88813505137905031680N+38018445536403763200)}{6227020800}.
\]
Set
\[
\Delta:=m_{N+14}^{(N)}m_{N+12}^{(N)}-(m_{N+13}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16,17,18\), direct substitution of the
hook-length formulas gives
\[
\begin{gathered}
\Delta=15209125007422457014538885532181,\quad
\Delta=17014764831940166467918891896443376,\\
\Delta=11954044471662473620117009488464470109,\quad
\Delta=6898468500189873697028158457457230556156,\\
\Delta=3795269633941222591232069486722596998316938,\\
\Delta=2151373239100398136354898202598028311909426671,\\
\Delta=1299944721670140108414492412514731043280703843509,\\
\Delta=846292474083657130704814885281988942745162196347679,\\
\Delta=594275308589286517790616682712447549552838159652345406,\\
\Delta=449294432061516979602060891762995126213468455859684506780,\\
\Delta=364829723736391211695384917061819279725169866349053208057200,\\
\Delta=317402159612922117541545195745155599331206068564469210223430544,\\
\Delta=295184041225304382678750569421169750621410953864780549842784015914,
\end{gathered}
\]
respectively.

Now assume \(N\ge19\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+12}^{(N)},m_{N+13}^{(N)},m_{N+14}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{13}\):
the base case is \(D_{19}=44750731559645106\ge19^{13}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{13}
\]
for \(r\ge19\) once \(D_r\ge r^{13}\).

At \(A=N^{13}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)\cdots(N+13)}{6227020800}
\bigl(12454041599N^{25}+971415244670N^{24}
+33837631019348N^{23}+694561899740150N^{22}
+9333843372213932N^{21}+86264587846362068N^{20}
+560310206577993443N^{19}+2566406201218551692N^{18}
+8184881128552111412N^{17}+17610199456814017196N^{16}
+24055720419964618112N^{15}+18512217533091647144N^{14}
+5965492183593410249N^{13}-323973057692866N^{12}
-3210882481160536N^{11}-18321997142978386N^{10}
-69669460969051372N^9-191082632572533988N^8
-434421646112608011N^7-1077043312773247248N^6
-3169712023591052160N^5-8542059089644997280N^4
-17406133220569715856N^3-24595612863291130752N^2
-21672358008305852160N-9020441693096064000\bigr)>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)\cdots(N+13)}{6227020800}
\bigl(12454041599N^{25}+971415244670N^{24}
+33837631019348N^{23}+694561899740150N^{22}
+9333843372213932N^{21}+86264587846362068N^{20}
+560310206577993443N^{19}+2566406201218551692N^{18}
+8184881128552111412N^{17}+17610199456814017196N^{16}
+24055720419964618112N^{15}+18512217533091647144N^{14}
+5965492183593410249N^{13}-323973057692866N^{12}
-3210857573077336N^{11}-18320104128655186N^{10}
-69605546827560172N^9-189819518765378788N^8
-418221677880160011N^7-935823049409868048N^6
-2318501818369676160N^5-5007759801547819680N^4
-7546268909287130256N^3-7101722488831431552N^2
-4052948205635976960N-1478603920519756800\bigr)>0.
\]
Verification. After substituting \(N=s+6\), the two derivative numerators
have positive coefficients; each shifted numerator has degree \(38\), and the
minimum coefficient is \(12454041599\). Thus the derivative is positive for
\(A\ge N^{13}\). The lower value \(\Delta|_{A=N^{13}}\) has denominator
\(2982752926433280000\) in both parities. After substituting \(N=s+6\), the
two shifted numerators have positive coefficients, degree \(51\), and minimum
coefficient \(2982752925954278400\). Hence \(\Delta>0\). \(\square\)

## Proposition 38 (fifteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+13}^{(N)}\le\rho_{N+14}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge15\), Proposition 30 with \(j=15\) gives
\[
\begin{aligned}
m_{N+15}^{(N)}
={}&D_{N+15}
-\frac{(N+14)}{87178291200}
\bigl(N^{27}+182N^{26}+15652N^{25}+845117N^{24}\\
&+32101888N^{23}+910822640N^{22}+19998674839N^{21}
+347260018820N^{20}\\
&+4829644115296N^{19}+54114708174527N^{18}
+488287987296382N^{17}\\
&+3522742320543992N^{16}+20005270556076427N^{15}
+87104044310471822N^{14}\\
&+281849275017488152N^{13}+714027223326268943N^{12}
+2409281399682671332N^{11}\\
&+15722868296890014260N^{10}+97477849664268679749N^9
+439151669759458120608N^8\\
&+1446381580164829942140N^7+3624936774990311937429N^6
+7288879718528816952726N^5\\
&+12505498175874167902620N^4+18729349570303453623816N^3\\
&+22736991370825767839040N^2
+18891420970769515958400N\\
&+7778222300950028928000\bigr).
\end{aligned}
\]
Set
\[
\Delta:=m_{N+15}^{(N)}m_{N+13}^{(N)}-(m_{N+14}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16,17,18,19,20\), direct substitution
of the hook-length formulas gives
\[
\begin{gathered}
\Delta=3991738155322404286534175980834919,\quad
\Delta=5898203365383335273495293855085752664,\\
\Delta=5040283082217466381135182525996998449231,\quad
\Delta=3358455318815892241759813550104154244234210,\\
\Delta=2066046376131019033988773352286445063854396123,\\
\Delta=1285431991085406200840453706920888567179729024575,\\
\Delta=844076626835699601411276046726329072426682329869391,\\
\Delta=593964868226080053898149040632146112754689049662865446,\\
\Delta=449253858958749376621469392789046501424543150810317751145,\\
\Delta=364824713378873099503291213872898768754673407151675754764242,\\
\Delta=317401569079244644355259058420717219424324102710527553226435400,\\
\Delta=295183974256472300604972113684239429123851448490317783207619554456,\\
\Delta=292822565015749493359453262085177869944067058245358305449162399864151,\\
\Delta=309220635996501608275559151967161054496212942045035334350824518831013440,\\
\Delta=346945554424780708882614572901870497464138125683994340510998555005278405767,
\end{gathered}
\]
respectively.

Now assume \(N\ge21\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+13}^{(N)},m_{N+14}^{(N)},m_{N+15}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{14}\):
the base case is \(D_{21}=18795307255050944540\ge21^{14}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{14}
\]
for \(r\ge21\) once \(D_r\ge r^{14}\), because \(r^{15}\ge(r+1)^{14}\)
on this range. The coefficient of \(A^2\) in
\(\Delta\) is
\[
(N+14)\bigl((N+1)\cdots(N+13)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

At \(A=N^{14}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)\cdots(N+14)}{87178291200}
\bigl(174356582399N^{27}+15866448998246N^{26}
+650524408923298N^{25}+15882315446901449N^{24}
+257084073105635776N^{23}\\
+2906463726521627636N^{22}+23544497925544776037N^{21}
+137906144949067441232N^{20}\\
+581674269185793765982N^{19}+1736191189127146867115N^{18}
+3541836530331088740478N^{17}\\
+4637312560373711455328N^{16}+3452741137697038073629N^{15}
+1085721297016662935030N^{14}\\
-15439091131562702N^{13}-167075509815868381N^{12}
-1064876337758914772N^{11}\\
-4570578017070465400N^{10}-14137817620029393433N^9
-34726122683157742980N^8\\
-85376787799436196570N^7-256571828517714770151N^6
-789639704500039800498N^5\\
-1978749609505276028148N^4-3692850298930072823544N^3\\
-4874352718337913341376N^2-4104208141068604049280N
-1657953189071661081600\bigr)>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)\cdots(N+14)}{87178291200}
\bigl(174356582399N^{27}+15866448998246N^{26}
+650524408923298N^{25}+15882315446901449N^{24}
+257084073105635776N^{23}\\
+2906463726521627636N^{22}+23544497925544776037N^{21}
+137906144949067441232N^{20}\\
+581674269185793765982N^{19}+1736191189127146867115N^{18}
+3541836530331088740478N^{17}\\
+4637312560373711455328N^{16}+3452741137697038073629N^{15}
+1085721297016662935030N^{14}\\
-15439091131562702N^{13}-167075161102703581N^{12}
-1064845302287247572N^{11}\\
-4569338690482766200N^{10}-14108501652977822233N^9
-34269436415001224580N^8\\
-80451309894043351770N^7-218952797712480491751N^6
-585239257116388882098N^5\\
-1197469038713633506548N^4-1653440529034206455544N^3\\
-1443792932948869834176N^2-791889827992261879680N
-285338711673067852800\bigr)>0.
\]
Verification. After substituting \(N=s+6\), the two derivative numerators
have positive coefficients; each shifted numerator has degree \(41\), and the
minimum coefficient is \(174356582399\). Thus the derivative is positive for
\(A\ge N^{14}\). The lower value \(\Delta|_{A=N^{14}}\) has denominator
\(542861032610856960000\) in both parities. After substituting \(N=s+6\), the
two shifted numerators have positive coefficients, degree \(55\), and minimum
coefficient \(542861032604629939200\). Hence \(\Delta>0\). \(\square\)

## Proposition 39 (sixteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+14}^{(N)}\le\rho_{N+15}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge16\), Proposition 30 with \(j=16\) gives
\[
\begin{aligned}
m_{N+16}^{(N)}
={}&D_{N+16}
-\frac{(N+15)}{1307674368000}
\bigl(N^{29}+210N^{28}+20930N^{27}+1316175N^{26}\\
&+58555224N^{25}+1958450676N^{24}+51072724755N^{23}
+1062699967680N^{22}\\
&+17901419606070N^{21}+246171056576445N^{20}
+2772207937708482N^{19}\\
&+25516567687709280N^{18}+190495008548308995N^{17}
+1136097197103512130N^{16}\\
&+5278602639316005162N^{15}+18497789166188188773N^{14}
+49892797028411708700N^{13}\\
&+157539556679896110720N^{12}+1010468875443173529617N^{11}\\
&+6977618320983353096196N^{10}+35875336937824636399654N^9\\
&+135493652883885672496575N^8+389572091422124368679922N^7\\
&+894072062610901708559844N^6
+1745545274863806479962536N^5\\
&+3058057165677662636967744N^4
+4706161889168267162845632N^3\\
&+5696450026050784372791552N^2
+4590790794846222727864320N\\
&+1814379810473622712320000\bigr).
\end{aligned}
\]
Set
\[
\Delta:=m_{N+16}^{(N)}m_{N+14}^{(N)}-(m_{N+15}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21\), direct
substitution of the hook-length formulas gives
\[
\begin{gathered}
\Delta=1115481831981990204496480262186834465,\quad
\Delta=2194708106876995361269892500092216198716,\\
\Delta=2293537578108960944158891686734298244450241,\quad
\Delta=1769453246942105629079160680139924582432525103,\\
\Delta=1217951811039678309054340345767022745700429476529,\\
\Delta=830859267100433405333926788099402392737189906592175,\\
\Delta=591660961924964874782476800661153456785285973766493209,\\
\Delta=448888014877920541630932449200488951682815982050930662700,\\
\Delta=364770855144474401166911135116342250325661386553820192293313,\\
\Delta=317394118049093606003111058599154015114821077572442778071375559,\\
\Delta=295182995153875771100678711714873976122974096215896863311208818900,\\
\Delta=292822441759660595011593011527308269859947542508219387767018110179324,\\
\Delta=309220621027785476070705423458016136016421351749947220532451716922439513,\\
\Delta=346945552660968984216984255900562603811274547652462445340345961099178935232,\\
\Delta=412865209655962161124927369627914126381502823515186640509684030744223806570941,\\
\Delta=520210164409184248207915873788176328958526226998982750879254684953948126966639580,
\end{gathered}
\]
respectively.

Now assume \(N\ge22\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+14}^{(N)},m_{N+15}^{(N)},m_{N+16}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{15}\):
the base case is \(D_{22}=413496759611120779881\ge22^{15}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{15}
\]
for \(r\ge22\) once \(D_r\ge r^{15}\), because \(r^{16}\ge(r+1)^{15}\)
on this range. The coefficient of \(A^2\) in \(\Delta\) is
\[
(N+15)\bigl((N+1)\cdots(N+14)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

At \(A=N^{15}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)\cdots(N+15)}{1307674368000}
\bigl(2615348735999N^{29}+274611617279820N^{28}
+13089820423664740N^{27}+374844857586390975N^{26}\\
+7191547340739768006N^{25}+97584611257707512484N^{24}
+963524851600498245825N^{23}\\
+7012936740800288467650N^{22}+37685404502613900103500N^{21}\\
+148194464595082452511965N^{20}
+417727699226732356522788N^{19}\\
+813345367342819161210360N^{18}
+1025626778497995734759805N^{17}\\
+741361505978465264370360N^{16}
+228001586046060871750488N^{15}\\
-798783614945242083N^{14}-9360560540154422130N^{13}\\
-66268645835071134000N^{12}-319045152979705839557N^{11}\\
-1110577140267104098446N^{10}-2994203278814752607824N^9\\
-7466582650364112120105N^8-22233463153085288041392N^7\\
-74742380604263634003924N^6-220208076751774694354136N^5\\
-505298489664009657405024N^4-873882142435287888367392N^3\\
-1094799461451596123980032N^2-892041587287113561822720N\\
-352296665209348763904000\bigr)>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)\cdots(N+15)}{1307674368000}
\bigl(2615348735999N^{29}+274611617279820N^{28}
+13089820423664740N^{27}+374844857586390975N^{26}\\
+7191547340739768006N^{25}+97584611257707512484N^{24}
+963524851600498245825N^{23}\\
+7012936740800288467650N^{22}+37685404502613900103500N^{21}\\
+148194464595082452511965N^{20}
+417727699226732356522788N^{19}\\
+813345367342819161210360N^{18}
+1025626778497995734759805N^{17}\\
+741361505978465264370360N^{16}
+228001586046060871750488N^{15}\\
-798783614945242083N^{14}-9360555309456950130N^{13}\\
-66268107073231518000N^{12}-319020045631840239557N^{11}\\
-1109877142177913698446N^{10}-2981196631711574879824N^9\\
-7296796365470350584105N^8-20634827330874276457392N^7\\
-63776377424409671763924N^6-165564774239333363666136N^5\\
-310638423249326732349024N^4-394147693976320298479392N^3\\
-324697131823751854924032N^2-175865524382718156510720N\\
-64047625000105114368000\bigr)>0.
\]
Verification. After substituting \(N=s+6\), the two derivative numerators
have positive coefficients; each shifted numerator has degree \(44\), and the
minimum coefficient is \(2615348735999\). Thus the derivative is positive for
\(A\ge N^{15}\). The lower value \(\Delta|_{A=N^{15}}\) has denominator
\(114000816848279961600000\) in both parities. After substituting \(N=s+6\),
the two shifted numerators have positive coefficients, degree \(59\), and
minimum coefficient \(114000816848192783308800\). Hence \(\Delta>0\).
\(\square\)

## Proposition 40 (seventeenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+15}^{(N)}\le\rho_{N+16}^{(N)}.
\]

Proof.
Let \(D_r:=!r\). For \(N\ge17\), Proposition 30 with \(j=17\) gives
\[
\begin{aligned}
m_{N+17}^{(N)}
={}&D_{N+17}
-\frac{(N+16)}{20922789888000}
\bigl(N^{31}+240N^{30}+27440N^{29}+1987800N^{28}\\
&+102357024N^{27}+3983660016N^{26}+121630116660N^{25}
+2984270026920N^{24}\\
&+59777115125520N^{23}+987372269877240N^{22}
+13522135059409392N^{21}\\
&+153752055595254480N^{20}+1447319691038873910N^{19}\\
&+11190488414741125800N^{18}
+70039824986203849752N^{17}\\
&+346532977382762175048N^{16}
+1312005508728640499520N^{15}\\
&+3793333103365811265840N^{14}
+11393657028699311769812N^{13}\\
&+69064262866672092293976N^{12}
+522392080045370037989824N^{11}\\
&+3043147936884611344150440N^{10}
+13083027545954381661929712N^9\\
&+42830681734033508053184784N^8
+111317114438048692678778001N^7\\
&+243604947576316366029810024N^6
+480646577255806671959180232N^5\\
&+873519583247772703908748992N^4
+1364219785253014229635837200N^3\\
&+1619826464400833718564278400N^2
+1257746087215564538508000000N\\
&+477839802926302619362560000\bigr).
\end{aligned}
\]
Set
\[
\Delta:=m_{N+17}^{(N)}m_{N+15}^{(N)}-(m_{N+16}^{(N)})^2.
\]
For \(N=6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22\), direct
substitution of the hook-length formulas gives
\[
\begin{gathered}
\Delta=330238622677566348606835529204103156175,\quad
\Delta=872102308788948712677575406149626520015600,\\
\Delta=1120897840405634629076647391458073994424748271,\quad
\Delta=1004492636116933773371918637879014179232292817386,\\
\Delta=774506482407336683257416909538219299797550880712831,\\
\Delta=579001310665093245201283446824017396486304717642532350,\\
\Delta=446378739088979085292195662379817155835016450582843989695,\\
\Delta=364320932296595309334011745791227707612957134448588564742008,\\
\Delta=317319773562757293214856889348445193189399665915240975886796911,\\
\Delta=295171511125069989150120762186237439520916502129632642598907018446,\\
\Delta=292820764604640553121509660575263824764065214300349450305356457570525,\\
\Delta=309220387351521072341310666724453602963008431628290319054633424308542136,\\
\Delta=346945521370138340723626312748425317453821787075501638824669996145750232199,\\
\Delta=412865205604433925549492731250365517300106849555696632856941796196709799943664,\\
\Delta=520210163899359629670611798655554663543427312809006038795749239385237602521118115,\\
\Delta=692919938959339933478628583702810798535977424127318348436216303866489144750508747024,\\
\Delta=974245434260715062491503730932581213331379514237844661888200654048505103562014649829055,
\end{gathered}
\]
respectively.

Now assume \(N\ge23\), and set \(A:=D_N\). Verification. Substituting the
formulas for \(m_{N+15}^{(N)},m_{N+16}^{(N)},m_{N+17}^{(N)}\), and using
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
gives a quadratic in \(A\) in each parity of \(N\). Also \(A=D_N\ge N^{16}\):
the base case is \(D_{23}=9510425471055777937262\ge23^{16}\), and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge (r+1)^{16}
\]
for \(r\ge23\) once \(D_r\ge r^{16}\), because \(r^{17}\ge(r+1)^{16}\)
on this range. The coefficient of \(A^2\) in \(\Delta\) is
\[
(N+16)\bigl((N+1)\cdots(N+15)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

At \(A=N^{16}\), the derivative of the even-parity quadratic is
\[
\frac{(N+1)\cdots(N+16)}{20922789888000}
\bigl(41845579775999N^{31}+5021469573119792N^{30}
+275343914926059520N^{29}+9139074623077131640N^{28}\\
+205027523273189311296N^{27}+3287325141919492676352N^{26}\\
+38836704327890014886652N^{25}
+343452952244815732926360N^{24}\\
+2286071289948387049605120N^{23}
+11415608515470285615167640N^{22}\\
+42250314702481121141534928N^{21}
+113268173776411389699023424N^{20}\\
+211612917055454927695998330N^{19}
+258012212464983227463314040N^{18}\\
+181574790395932207263473208N^{17}
+54720388832911145409135336N^{16}\\
-44745958376129003424N^{15}-561532991274452622240N^{14}\\
-4393674245553634916612N^{13}-23594314157633245077272N^{12}\\
-92082824025797215613392N^{11}-274790253112528303821112N^{10}\\
-712342786518700680553392N^9-2075992024109602417696800N^8\\
-7327945473374217304663233N^7-24632948784634560967442472N^6\\
-66918252700047442682120904N^5-141369341364532130519210688N^4\\
-230222234358381424024455696N^3
-277922726996275632501072000N^2\\
-221301411393866989730630400N-85819436554715217978624000\bigr)>0,
\]
and the derivative of the odd-parity quadratic is
\[
\frac{(N+1)\cdots(N+16)}{20922789888000}
\bigl(41845579775999N^{31}+5021469573119792N^{30}
+275343914926059520N^{29}+9139074623077131640N^{28}\\
+205027523273189311296N^{27}+3287325141919492676352N^{26}\\
+38836704327890014886652N^{25}
+343452952244815732926360N^{24}\\
+2286071289948387049605120N^{23}
+11415608515470285615167640N^{22}\\
+42250314702481121141534928N^{21}
+113268173776411389699023424N^{20}\\
+211612917055454927695998330N^{19}
+258012212464983227463314040N^{18}\\
+181574790395932207263473208N^{17}
+54720388832911145409135336N^{16}\\
-44745958376129003424N^{15}-561532907583293070240N^{14}\\
-4393664369996807780612N^{13}-23593783137225887637272N^{12}\\
-92065598292882425213392N^{11}-274414147217471764173112N^{10}\\
-706504610655637841257392N^9-2009665142579723463712800N^8\\
-6768816825165911128663233N^7-21126815181270548738834472N^6\\
-50689299034426596516104904N^5-86975174249664747067562688N^4\\
-102764329454183694146439696N^3
-81639350879122681442640000N^2\\
-44807171653422348483398400N-16639666903743521654016000\bigr)>0.
\]
Verification. After substituting \(N=s+6\), the two derivative numerators
have positive coefficients; each shifted numerator has degree \(47\), and the
minimum coefficient is \(41845579775999\). Thus the derivative is positive for
\(A\ge N^{16}\). The lower value \(\Delta|_{A=N^{16}}\) has denominator
\(27360196043587190784000000\) in both parities. After substituting \(N=s+6\),
the two shifted numerators have positive coefficients, degree \(63\), and
minimum coefficient \(27360196043585883109632000\). Hence \(\Delta>0\).
\(\square\)

## Proposition 41 (quadratic transition window)

For \(N\ge25\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\qquad
\left(N+16\le k\le \left\lfloor\frac{N^2}{11}\right\rfloor-2\right).
\]

Proof.
Let
\[
K:=\left\lfloor\frac{N^2}{11}\right\rfloor,\qquad
\theta:=\frac{121}{176}.
\]
For \(a\ge N+1\), Proposition 21 and RSK give
\[
a!-E_a^{(N)}
=\#\{\pi\in S_a:\operatorname{lds}(\pi)>N\}.
\]
If \(\operatorname{lds}(\pi)>N\), then \(\pi\) has a decreasing subsequence
of length \(N+1\). For a fixed \((N+1)\)-element set of positions, exactly
one of the \((N+1)!\) relative value orders is decreasing, hence
\[
0\le a!-E_a^{(N)}
\le a!\frac{\binom a{N+1}}{(N+1)!}.
\]
For \(r\le K\), comparing the binomial transform in Proposition 21 with the
derangement transform gives
\[
\begin{aligned}
\left|m_r^{(N)}-D_r\right|
&\le\sum_{a=N+1}^r\binom ra(a!-E_a^{(N)})\\
&\le
r!\frac{\binom r{N+1}}{(N+1)!}\sum_{b=0}^{r-N-1}\frac1{b!}\\
&\le
3r!\frac{\binom K{N+1}}{(N+1)!}.
\end{aligned}
\]
Since \(K\le N^2/11<(N+1)^2/11\), while
\((N+1)!\ge((N+1)/e)^{N+1}\) and \(e<11/4\), this gives
\[
\frac{\binom K{N+1}}{(N+1)!}
\le \frac{K^{N+1}}{((N+1)!)^2}
\le \theta^{N+1}.
\]
Also \(D_r\ge r!/3\) for \(r\ge2\), because
\[
\frac{D_r}{r!}=\sum_{i=0}^r\frac{(-1)^i}{i!}\ge 1-1+\frac12-\frac16=\frac13.
\]
Thus, for \(r\le K\),
\[
\left|m_r^{(N)}-D_r\right|\le \tau D_r,\qquad
\tau:=9\theta^{N+1}.
\]

Now let \(N+16\le k\le K-2\). Put
\[
e_r:=m_r^{(N)}-D_r\qquad(r=k,k+1,k+2).
\]
The derangement-ratio gap from Proposition 23 gives
\[
D_{k+2}D_k-D_{k+1}^2
=D_kD_{k+1}\left(\frac{D_{k+2}}{D_{k+1}}-\frac{D_{k+1}}{D_k}\right)
\ge \frac7{18}D_kD_{k+1}.
\]
Expanding
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2
\]
around \(D_{k+2}D_k-D_{k+1}^2\), and using
\(|e_r|\le\tau D_r\), gives
\[
\begin{aligned}
&m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2\\
&\quad\ge D_{k+2}D_k-D_{k+1}^2
-\tau(2D_{k+2}D_k+2D_{k+1}^2)
-\tau^2(D_{k+2}D_k+D_{k+1}^2).
\end{aligned}
\]
Since \(D_{j+1}=(j+1)D_j+(-1)^{j+1}\),
\[
D_{k+1}\le(k+2)D_k,\qquad D_{k+2}\le(k+3)D_{k+1}.
\]
Therefore
\[
\begin{aligned}
&m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2\\
&\quad\ge
D_kD_{k+1}
\left[
\frac7{18}
-\tau(4k+10)-\tau^2(2k+5)
\right].
\end{aligned}
\]
It remains to check the bracket. Because \(k\le K\le N^2/11\),
\[
\tau(4k+10)+\tau^2(2k+5)
\le
9\theta^{N+1}\left(\frac{4N^2}{11}+10\right)
+81\theta^{2N+2}\left(\frac{2N^2}{11}+5\right).
\]
Verification. Each of the two terms on the right is decreasing for
\(N\ge25\), since
\[
\theta\frac{4(N+1)^2/11+10}{4N^2/11+10}<1,\qquad
\theta^2\frac{2(N+1)^2/11+5}{2N^2/11+5}<1.
\]
At \(N=25\), direct rational arithmetic gives the right side
\[
<\frac{129}{1000}<\frac7{18}.
\]
Hence the bracket is positive, so
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. The positivity paragraph in Proposition 42 gives
\(m_k^{(N)},m_{k+1}^{(N)}>0\), since \(k\ge N+16\ge41\). Therefore this is
equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 42 (even transition ratios)

For \(N\ge4\) and every even \(k\ge2\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}.
\]

Proof.
First \(m_r^{(N)}>0\) for every \(r\ge2\). By character orthogonality,
\[
m_r^{(N)}=\dim\left((\mathfrak{sl}_N)^{\otimes r}\right)^{SU(N)}.
\]
The Killing form gives a nonzero invariant tensor in degree \(2\), and
\((x,y,z)\mapsto\kappa([x,y],z)\) gives a nonzero invariant tensor in
degree \(3\). Since every \(r\ge2\) is \(2a+3b\) for some
\(a,b\ge0\), tensor products of these two invariant tensors give
\(m_r^{(N)}>0\).

Let \(X=\chi_{\rm adj}(U)\), with \(U\) Haar-distributed in \(SU(N)\), and
write \(k=2q\). The functions \(X^q\) and \(X^{q+1}\) are real
square-integrable, so Cauchy-Schwarz gives
\[
\left(m_{k+1}^{(N)}\right)^2
=\left(E[X^qX^{q+1}]\right)^2
\le E[X^{2q}]E[X^{2q+2}]
=m_k^{(N)}m_{k+2}^{(N)}.
\]
Verification. The exponents are \(2q=k\), \(2q+1=k+1\), and
\(2q+2=k+2\), and the first paragraph gives positive denominators.
Because \(m_k^{(N)}\) and \(m_{k+1}^{(N)}\) are positive, this is
equivalent to \(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 43 (endpoint tail for odd transition ratios)

Fix \(N\ge6\), put \(A=N^2-1\), and let \(\mu_N\) be the distribution of
\(X=\chi_{\rm adj}(U)\) under Haar measure on \(SU(N)\). Define
\[
p_0(N):=\mu_N([3A/4,A]),\qquad
p_1(N):=\mu_N([A/2,5A/8]),
\]
and
\[
K_{\rm end}(N):=
\min\left\{K\ge3:
\frac{p_0(N)p_1(N)A^2}{64}\left(\frac{3A}{8}\right)^K>N^4
\right\}.
\]
Then, for every odd \(k\ge K_{\rm end}(N)\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}.
\]

Proof.
The constants are defined by nonempty positive-measure trace bands. The first
band has positive Haar measure because \(X(e)=A\) and \(X\) is continuous.
For the second band, choose
\[
s=\sqrt{1+\frac{9A}{16}}\in(0,N).
\]
If \(N\) is even, the diagonal matrix with \(N/2\) pairs
\((e^{it},e^{-it})\) has determinant \(1\) and trace \(N\cos t\). If
\(N\) is odd, the diagonal matrix with \((N-1)/2\) such pairs and one
additional \(1\) has determinant \(1\) and trace \(1+(N-1)\cos t\). In
both cases one can choose \(t\) so that the trace is \(s\), hence
\[
X=s^2-1=\frac{9A}{16}\in(A/2,5A/8).
\]
Continuity gives a nonempty open subset of the second band, so
\(p_0(N),p_1(N)>0\). Since \(A\ge35\), \(3A/8>1\), and therefore
\(K_{\rm end}(N)\) is finite.

For any \(k\ge0\),
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2
=\frac12\iint x^ky^k(x-y)^2\,d\mu_N(x)d\mu_N(y).
\]
Now assume \(k\) is odd. The contribution from
\([3A/4,A]\times[A/2,5A/8]\) and the symmetric rectangle is at least
\[
p_0(N)p_1(N)\left(\frac{3A}{4}\right)^k
\left(\frac{A}{2}\right)^k\left(\frac{A}{8}\right)^2
=\frac{p_0(N)p_1(N)A^2}{64}\left(\frac{3A^2}{8}\right)^k.
\]
The same-sign regions outside these two rectangles contribute
nonnegatively. The only possible negative contribution has one variable in
\((0,A]\) and the other in \([-1,0)\). Its absolute value is bounded by
\[
\int_{(0,A]\times[-1,0)}
x^k|y|^k(x+|y|)^2\,d\mu_N(x)d\mu_N(y)
\le (A+1)^2A^k=N^4A^k.
\]
Thus
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2
\ge
A^k\left[
\frac{p_0(N)p_1(N)A^2}{64}\left(\frac{3A}{8}\right)^k
-N^4
\right].
\]
For \(k\ge K_{\rm end}(N)\), the bracket is positive. Hence
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. The positivity paragraph in Proposition 42 gives positive
denominators for \(\rho_k^{(N)}\) and \(\rho_{k+1}^{(N)}\) because
\(k\ge K_{\rm end}(N)\ge3\). The displayed determinant is therefore
equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 44 (sharpened quadratic transition window)

For \(N\ge31\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\qquad
\left(N+16\le k\le \left\lfloor\frac{N^2}{10}\right\rfloor-2\right).
\]

Proof.
Put
\[
K:=\left\lfloor\frac{N^2}{10}\right\rfloor,\qquad
\theta:=\frac34.
\]
The RSK omission estimate in Proposition 41 gives, for \(r\le K\),
\[
\left|m_r^{(N)}-D_r\right|
\le
3r!\frac{\binom K{N+1}}{(N+1)!}.
\]
Since \(K\le N^2/10<(N+1)^2/10\), while
\((N+1)!\ge((N+1)/e)^{N+1}\) and \(e^2<15/2\), this gives
\[
\frac{\binom K{N+1}}{(N+1)!}
\le \frac{K^{N+1}}{((N+1)!)^2}
\le \theta^{N+1}.
\]
As in Proposition 41, \(D_r\ge r!/3\) for \(r\ge2\), so
\[
\left|m_r^{(N)}-D_r\right|\le\tau D_r,\qquad
\tau:=9\theta^{N+1}.
\]

Now let \(N+16\le k\le K-2\). Expanding around the derangement determinant
and using the same derangement-ratio gap and recurrence bounds as in
Proposition 41 gives
\[
\begin{aligned}
&m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2\\
&\quad\ge
D_kD_{k+1}
\left[
\frac7{18}
-\tau(4k+10)-\tau^2(2k+5)
\right].
\end{aligned}
\]
Because \(k\le K\le N^2/10\),
\[
\tau(4k+10)+\tau^2(2k+5)
\le
9\theta^{N+1}\left(\frac{2N^2}{5}+10\right)
+81\theta^{2N+2}\left(\frac{N^2}{5}+5\right).
\]
Verification. Each term on the right is decreasing for \(N\ge31\), since
\[
\theta\frac{2(N+1)^2/5+10}{2N^2/5+10}<1,\qquad
\theta^2\frac{(N+1)^2/5+5}{N^2/5+5}<1.
\]
At \(N=31\), direct rational arithmetic gives the right side
\[
<\frac{357}{1000}<\frac7{18}.
\]
Hence the bracket is positive, so
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. The positivity paragraph in Proposition 42 gives
\(m_k^{(N)},m_{k+1}^{(N)}>0\), since \(k\ge N+16\ge47\). Therefore this is
equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 45 (larger quadratic transition windows)

The following transition inequalities hold:
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\]
for every \(k\) in either of the two ranges
\[
N\ge56,\qquad
N+16\le k\le \left\lfloor\frac{N^2}{9}\right\rfloor-2,
\]
and
\[
N\ge202,\qquad
N+16\le k\le \left\lfloor\frac{N^2}{8}\right\rfloor-2.
\]

Proof.
Let \(c\in\{9,8\}\), put \(K=\lfloor N^2/c\rfloor\), and set
\[
\theta=\frac{15}{2c}.
\]
Thus \(\theta=5/6\) for \(c=9\), and \(\theta=15/16\) for \(c=8\).
The proof of Proposition 44 applies with this \(K\) and \(\theta\): since
\(K\le N^2/c<(N+1)^2/c\), \((N+1)!\ge((N+1)/e)^{N+1}\), and
\(e^2<15/2\), one has
\[
\frac{\binom K{N+1}}{(N+1)!}
\le\frac{K^{N+1}}{((N+1)!)^2}
\le\theta^{N+1}.
\]
Consequently, for \(N+16\le k\le K-2\),
\[
\begin{aligned}
&m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2\\
&\quad\ge
D_kD_{k+1}
\left[
\frac7{18}
-\tau(4k+10)-\tau^2(2k+5)
\right],
\qquad \tau:=9\theta^{N+1}.
\end{aligned}
\]
It remains only to check the bracket. Since \(k\le K\le N^2/c\),
\[
\tau(4k+10)+\tau^2(2k+5)
\le
9\theta^{N+1}\left(\frac{4N^2}{c}+10\right)
+81\theta^{2N+2}\left(\frac{2N^2}{c}+5\right).
\]
Verification. In the two rows below, each displayed upper bound is
decreasing for \(N\ge N_0\), because
\[
\theta\frac{4(N+1)^2/c+10}{4N^2/c+10}<1,\qquad
\theta^2\frac{2(N+1)^2/c+5}{2N^2/c+5}<1.
\]
At the lower endpoints, direct rational arithmetic gives
\[
\begin{array}{c|c|c|c}
c&\theta&N_0&
9\theta^{N_0+1}(4N_0^2/c+10)
+81\theta^{2N_0+2}(2N_0^2/c+5)\\ \hline
9&5/6&56&<388/1000<7/18\\
8&15/16&202&<376/1000<7/18 .
\end{array}
\]
Therefore the bracket is positive throughout the two stated ranges. The
positivity paragraph in Proposition 42 supplies the positive denominators,
so the determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 46 (eighteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+16}^{(N)}\le\rho_{N+17}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+18}^{(N)}m_{N+16}^{(N)}-(m_{N+17}^{(N)})^2.
\]
For \(N=6,7,\ldots,24\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&103118037378555716807221023542784948909350\\
7&368368471165820769740907882543444274381465088\\
8&585740106424241759282709474428257053041819353521\\
9&611922039973400113959633219298262606651362397475964\\
10&529354979263291023676548035611156314214340804113180509\\
11&433633449816186777146235126857400343587215337600068879500\\
12&361458697490377187380063902422353341904244753113850978684852\\
13&316742284561904929381189090876366785374628355415846701800709296\\
14&295064758388323286864055837813724966164503303315943149473398657949\\
15&292802409626245048500541418095947280698812948934106262363658568090844\\
16&309217416903020032106237613922674763970958939711675940773785181218572125\\
17&346945064580251921165883600732989930929879958383288852836677601772335849712\\
18&412865138337388747866404769504106566320537016833054903183514706799180472273149\\
19&520210154352596273126938930753959459228533147165102036706577763271667757471193020\\
20&692919937646549027925662341747635667588147721292329499083437820050264145092998817965\\
21&974245434085013769739339390505401838540269933549890570549910217592660397696655571141440\\
22&1443831733562066350382836606222078221399956415398153576900761831124868070245609432990952829\\
23&2252377504391049231253427228512461326270424513432867650379294344386322988832645469876627132316\\
24&3693899107205942385269006581615395843931902227557882167718785580692083416598810579847480569904056.
\end{array}
\]

Now assume \(N\ge25\), and set \(A:=D_N\). Proposition 30 with
\(j=16,17,18\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{17}\): the base case is
\[
D_{25}=5706255282633466762357224>25^{17},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{17}
\]
for \(r\ge25\) once \(D_r\ge r^{17}\), because
\(r^{18}\ge(r+1)^{17}\) on this range. The coefficient of \(A^2\) is
\[
(N+17)\bigl((N+1)\cdots(N+16)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{17}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(50\), denominator \(355687428096000\), and
minimum coefficient \(711374856191999\). Thus both derivatives are positive
for \(A\ge N^{17}\). At \(A=N^{17}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has degree
\(67\), denominator \(7441973323855715893248000000\), and minimum coefficient
\[
7441973323855694970458112000.
\]
Hence \(\Delta>0\). Since \(m_{N+16}^{(N)}\) and \(m_{N+17}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+16}^{(N)}\le\rho_{N+17}^{(N)}\). \(\square\)

## Proposition 47 (nineteenth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+17}^{(N)}\le\rho_{N+18}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+19}^{(N)}m_{N+17}^{(N)}-(m_{N+18}^{(N)})^2.
\]
For \(N=6,7,\ldots,25\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&33828008558668041059612268358422260969920875\\
7&164706021745838449301509467767189467930413279744\\
8&325947355157931913838486151342127109120067630632279\\
9&398515255155384750123775208461587196368643924045026536\\
10&387541151632832178661091433546118320044027772432760244175\\
11&347980197814847090981498747163135758585106025501551189525500\\
12&313324206392538264876890980154613656086502808057654710283547034\\
13&294291236672046408543394447323495717925666150100696328949192238304\\
14&292642930900306677088040188933981945168229507841124356573310744412111\\
15&309186983830463036364439794829926173261067410920170243269425576069893320\\
16&346939621962508477085229100321533951510142047211780700180962774532166921495\\
17&412864216992149694145182832231594563975760217023724496416015579899136518588304\\
18&520210005513351784781317914032861870424994403412668798287081753746517951448705391\\
19&692919914546961513034919415453662079693856863764299421951840728054739087760033983000\\
20&974245430621526351222321865225194392210802206819125302300741407972158742644415883174711\\
21&1443831733057989601827287413579628278507316627587225991305635168199043031743943827947723840\\
22&2252377504319550272105831876411617733113000844041251550142414924992191952529905034595015659415\\
23&3693899107196024441941949195941025640806782765621254506282947397733559743324794580019021261931514\\
24&6360894262607895949897537033656404651603858868179850144561009673200243986974433062070785726281880464\\
25&11487775038272196548546963783664510428797253825790041278905723390556195265979535383950311074259106515760.
\end{array}
\]

Now assume \(N\ge26\), and set \(A:=D_N\). Proposition 30 with
\(j=17,18,19\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{18}\): the base case is
\[
D_{26}=148362637348470135821287825>26^{18},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{18}
\]
for \(r\ge26\) once \(D_r\ge r^{18}\), because
\(r^{19}\ge(r+1)^{18}\) on this range. The coefficient of \(A^2\) is
\[
(N+18)\bigl((N+1)\cdots(N+17)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{18}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(53\), denominator \(6402373705728000\), and
minimum coefficient \(12804747411455999\). Thus both derivatives are
positive for \(A\ge N^{18}\). At \(A=N^{18}\), after substituting \(N=s+6\),
the two parity numerators for \(\Delta\) have positive coefficients; each
has degree \(71\), denominator
\[
2277243837099849063333888000000,
\]
and minimum coefficient
\[
2277243837099848707646459904000.
\]
Hence \(\Delta>0\). Since \(m_{N+17}^{(N)}\) and \(m_{N+18}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+17}^{(N)}\le\rho_{N+18}^{(N)}\). \(\square\)

## Proposition 48 (twentieth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+18}^{(N)}\le\rho_{N+19}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+20}^{(N)}m_{N+18}^{(N)}-(m_{N+19}^{(N)})^2.
\]
For \(N=6,7,\ldots,26\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&11618034532089060770125006895178670560956600625\\
7&77662574636323741403474319901803045565719014987631\\
8&192424410216017853650964205774858917080298059209794821\\
9&276480113537715674025620620283344567143271430662468483120\\
10&302935715649537216553648049861955691554802657401846586099977\\
11&298362915107662765306669478130444224408473920089633965420271900\\
12&290019790713130255340136617908226149937467657335399392342643886216\\
13&291561935912761280139584251910388335967365451256109031358461520391616\\
14&308939121662383441563894626659055666462542804296864465132898209855869105\\
15&346887268832940451131269290573030548051440921810623103392268026150044852592\\
16&412853897026536016185936607656717423579353583688223787322596392953359336772689\\
17&520208087071893631684312254956654407511124638532826618461601383822855931403814896\\
18&692919575358197927679875753361250624768955985052324527216504166690419606565800554169\\
19&974245373183453843176323032050887071080668309438379390290294227686857011782629270666800\\
20&1443831723687442996604114648925456792769682124386117254370494856309651628198406043756647937\\
21&2252377502839483091403010587643487000789592888471121857648508528947699320625387257170057215936\\
22&3693899106968736463849977140654072639792881064632068035017532977228979467047950112078862047416025\\
23&6360894262573837939260072800232650436799571077901180076984041040663134966747323450379839772686159631\\
24&11487775038267201181647001671341578703191253941900138064724493263832004387544013545890082006693455422860\\
25&21734870372410604846102445448364481709426359810257260642781301770735951307903862428204384719337785677140992\\
26&43035043337374364965539875371421197587930300926384277302415566995128260120569125191468478501311580836821085921.
\end{array}
\]

Now assume \(N\ge27\), and set \(A:=D_N\). Proposition 30 with
\(j=18,19,20\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{19}\): the base case is
\[
D_{27}=4005791208403350787852707224>27^{19},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{19}
\]
for \(r\ge27\) once \(D_r\ge r^{19}\), because
\(r^{20}\ge(r+1)^{19}\) on this range. The coefficient of \(A^2\) is
\[
(N+19)\bigl((N+1)\cdots(N+18)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{19}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(56\), denominator \(121645100408832000\),
and minimum coefficient \(243290200817663999\). Thus both derivatives are
positive for \(A\ge N^{19}\). At \(A=N^{19}\), after substituting \(N=s+6\),
the two parity numerators for \(\Delta\) have positive coefficients; each
has degree \(75\), denominator
\[
778817392288148379660189696000000,
\]
and minimum coefficient
\[
778817392288148373257815990272000.
\]
Hence \(\Delta>0\). Since \(m_{N+18}^{(N)}\) and \(m_{N+19}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+18}^{(N)}\le\rho_{N+19}^{(N)}\). \(\square\)

## Proposition 49 (twenty-first transition ratio)

For every \(N\ge6\),
\[
\rho_{N+19}^{(N)}\le\rho_{N+20}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+21}^{(N)}m_{N+19}^{(N)}-(m_{N+20}^{(N)})^2.
\]
For \(N=6,7,\ldots,27\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&4164303248310455205107862410570539013867725599375\\
7&38486531549673568789461790717538653338985250372725869\\
8&120099078562444683401472002185700722117768875230438094535\\
9&203669721131534499631048902022191123556206162139447339003508\\
10&252081682435044317937356005981767608178715109956345152774734499\\
11&272602770358626265305607457292878620703910411446645452093348396900\\
12&285979249404332071336758561097807013216415129398402159966800024544639\\
13&307363503039176680342672088491473746613151691440459734256098230002881600\\
14&346486551373343466291347703544993110788606512533030961139323914810035739343\\
15&412760452002938669291940890158098347046017056165446891921190803759559898759600\\
16&520187832946967994040217333027361363937826862435360313429395192540422560696406007\\
17&692915450226045476720831276956755927948934074294036438257104730649596134123432315064\\
18&974244576703436063672302053825968796290104261592164694453443232522380539618273246144783\\
19&1443831576829770622833489485585161457275013930949602395832081909796556269231179775388707320\\
20&2252377476822793357854949912285825742625652990239089117925828030869365682875529357980302241559\\
21&3693899102517536786519686816818175922997130124333368486946727350565826099374990318803321826763744\\
22&6360894261835122473552450625773014852678085853654722673782078393866614424455892377797166412316960775\\
23&11487775038147832688078583451281477835574258986819829981070019728891619041028271745091500452288987861645\\
24&21734870372391762629763831297578139588241249859258623478474091114093619823355580441561180890823822299922512\\
25&43035043337371451304993171098909330114852532947615817196674107485689970525856166602469331370701522712029464224\\
26&89082539708364695792726680137994092009043283361920103365823392011322504526317922234825237709550722026594513024729\\
27&192596450849485392426112477963452851493282631287330752823763136613574001275965833109891571821520751769837877044907350.
\end{array}
\]

Now assume \(N\ge28\), and set \(A:=D_N\). Proposition 30 with
\(j=19,20,21\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{20}\): the base case is
\[
D_{28}=112162153494581237110696612735>28^{20},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{20}
\]
for \(r\ge28\) once \(D_r\ge r^{20}\), because
\(r^{21}\ge(r+1)^{20}\) on this range. The coefficient of \(A^2\) is
\[
(N+20)\bigl((N+1)\cdots(N+19)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{20}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(59\), denominator \(2432902008176640000\),
and minimum coefficient \(4865804016353279999\). Thus both derivatives are
positive for \(A\ge N^{20}\). At \(A=N^{20}\), after substituting \(N=s+6\),
the two parity numerators for \(\Delta\) have positive coefficients; each
has degree \(79\), denominator
\[
295950609069496384270872084480000000,
\]
and minimum coefficient
\[
295950609069496384149226984071168000.
\]
Hence \(\Delta>0\). Since \(m_{N+19}^{(N)}\) and \(m_{N+20}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+19}^{(N)}\le\rho_{N+20}^{(N)}\). \(\square\)

## Proposition 50 (twenty-second transition ratio)

For every \(N\ge6\),
\[
\rho_{N+20}^{(N)}\le\rho_{N+21}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+22}^{(N)}m_{N+20}^{(N)}-(m_{N+21}^{(N)})^2.
\]
For \(N=6,7,\ldots,29\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&1553396016157092417387265305925683624437091534663125\\
7&19982893204220226915670654537199230731618694388734656375\\
8&78995409019157319000898038661308750379110907977514392656033\\
9&158818411977832660057890642007664366145090086951965668782843471\\
10&222670743425829172010697676492838134543778396743183932495149929806\\
11&264733340531834944434290670994343603457796060031930014353963704615900\\
12&299737116820008275819320209892062990187342982063839836636779559372250505\\
13&344092273630441128804937809531985586562514163708640080051077872602630697408\\
14&412086727253542264887039253018866360673214172345815638227041639160769030822661\\
15&520014791322360462617696476873963039079256747173497009655025986589132809129331440\\
16&692874303018657475996748560921797859139689794423667640112620761050439736611690653501\\
17&974235415032468960209796941344287413329867352565566229475739291512824383147085523680860\\
18&1443829649060530036518719524989804376971940457753976267352663444083121950669192220144467437\\
19&2252377090561893777744132271352324998239758423323575185600354176675439094982005996104964939484\\
20&3693899028351287196099393339177077024570774356964289468751075524284145384231902468359721908209161\\
21&6360894248115035800039114648444804837074745834919645170013292325096098439659313228113636600017985776\\
22&11487775035691348190746103647065291496996729082434141545786132220587468263239949047499953873318794196650\\
23&21734870371964415067612358473311641502942229945915531815417396645899915468360484442393020077961603824807063\\
24&43035043337298968716834077982374557885622584866942426835261714357055294667577498071874335049565444337184331584\\
25&89082539708352674323042752508118788606504509624480697368344559873610436053993997171234782082437266827002550250151\\
26&192596450849483437722917648555634411585968533883951891990789967618423384497266251126795974118852201405561432351402726\\
27&434497593116438876720613626157046235262318684823605128771520462697577344019421046729214176037457118266998325403956625295\\
28&1021938339009864945602955310999847884459601238129921835902317892480657794719033139733449489443318993708336290777683023913569\\
29&2503748930574169232841741933797851293303554218852221277058321415728111893754771331066430862697366497176850496239514297808147775.
\end{array}
\]

Now assume \(N\ge30\), and set \(A:=D_N\). Proposition 30 with
\(j=20,21,22\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{21}\): the base case is
\[
D_{30}=97581073836835777732377428235481
>10460353203000000000000000000000=30^{21},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{21}
\]
for \(r\ge30\) once \(D_r\ge r^{21}\), because
\(r^{22}\ge(r+1)^{21}\) on this range. The coefficient of \(A^2\) is
\[
(N+21)\bigl((N+1)\cdots(N+20)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{21}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(62\), denominator \(51090942171709440000\),
and minimum coefficient \(102181884343418879999\). Thus both derivatives
are positive for \(A\ge N^{21}\). At \(A=N^{21}\), after substituting
\(N=s+6\), the two parity numerators for \(\Delta\) have positive
coefficients; each has degree \(83\), denominator
\[
124299255809188481393766275481600000000,
\]
and minimum coefficient
\[
124299255809188481391333373473423360000.
\]
Hence \(\Delta>0\). Since \(m_{N+20}^{(N)}\) and \(m_{N+21}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+20}^{(N)}\le\rho_{N+21}^{(N)}\). \(\square\)

## Proposition 51 (twenty-third transition ratio)

For every \(N\ge6\),
\[
\rho_{N+21}^{(N)}\le\rho_{N+22}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+23}^{(N)}m_{N+21}^{(N)}-(m_{N+22}^{(N)})^2.
\]
For \(N=6,7,\ldots,30\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&601517467790849111561440559514188269937261482442059375\\
7&10840298082898314998254347408151662075332290181668436528832\\
8&54597392463047919561806638968746474663783420788982593837048391\\
9&130720224659562946970515233542291251039469585455426864108835066944\\
10&208237811559110255774753142812765221503211290576239137582622371055559\\
11&272605442507777507289589988585438165479286788986633645708978358927290500\\
12&333210159556217156252832099386513289582568009143508410481159638706261318087\\
13&408295427725102157893120909642897275832398168594027945628750268155038062530896\\
14&518837166192212217806568778848015086269204158362128857600971427309550429035544151\\
15&692541909293239396664910500850858637825679931200933125189970589697702884420523259224\\
16&974148889688083630398277946842668576818564207429017280950037059170883874868923076574695\\
17&1443808630819795396609486412914481595541245914459017661393495915908992520580136815294073824\\
18&2252372280307116573173149644150792039331523347781549791192573329339701292766111688743642764751\\
19&3693897982949021267536717881793907587197452994480872563583432626530519686196913944753388992189116\\
20&6360894030944535261129649679135252113865611665704843659446029172050832055637328824551453614580278719\\
21&11487774992327018079760613315660098139963079071332808632075770325738097544383291919278974139054971285780\\
22&21734870363601937720035767270165394166625788386579348871902082486107165071040921387271255925547728093547825\\
23&43035043335735189456725124570475553155873912852039843974850252845756769614822834304101787933580028831356188226\\
24&89082539708068106673927381454414032117699805257025498270666569258769645595741167915645740933751232867770010338516\\
25&192596450849432890078881447730213764341151960667039088508647462728727311694607231342561544402916108210685426551686919\\
26&434497593116430088645807506832247774756979500397480695770268908169280539272253455422397607992443154025480408055665569335\\
27&1021938339009863446615087771200348605503154158064011141013796399040116691656778908730518586397848973291245043037703074124369\\
28&2503748930574168981461683440636725696506923492001386724424673944819489251136825553689758243154789022389257937015192997145300079\\
29&6384559772964131521000544926551670927027415844465088347330634866666124007114923664605671259447075804094935798363101361151060341400\\
30&16931852517900876900047408336592063965039060025093776495299190920418293965377679524979181277941352328966500543665508448880057263756949.
\end{array}
\]

Now assume \(N\ge31\), and set \(A:=D_N\). Proposition 30 with
\(j=21,22,23\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{22}\): the base case is
\[
D_{31}=3025013288941909109703700275299910
>645590698195138073036733040138561=31^{22},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{22}
\]
for \(r\ge31\) once \(D_r\ge r^{22}\), because
\(r^{23}\ge(r+1)^{22}\) on this range. The coefficient of \(A^2\) is
\[
(N+22)\bigl((N+1)\cdots(N+21)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{22}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(65\), denominator
\(1124000727777607680000\), and minimum coefficient
\(2248001455555215359999\). Thus both derivatives are positive for
\(A\ge N^{22}\). At \(A=N^{22}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(87\), denominator
\[
57426256183845078403920019272499200000000,
\]
and minimum coefficient
\[
57426256183845078403868928330327490560000.
\]
Hence \(\Delta>0\). Since \(m_{N+21}^{(N)}\) and \(m_{N+22}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+21}^{(N)}\le\rho_{N+22}^{(N)}\). \(\square\)

## Proposition 52 (twenty-fourth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+22}^{(N)}\le\rho_{N+23}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+24}^{(N)}m_{N+22}^{(N)}-(m_{N+23}^{(N)})^2.
\]
For \(N=6,7,\ldots,31\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&241235142874756187787250564866640135491848651641884985625\\
7&6128394822262633184046355117809754112860257751626144684012864\\
8&39543502945423603066280271555632007749852210347034959407382851109\\
9&113263343307808015853453699895982579131451140052026999470336355447536\\
10&205655076223731328351301001923473230582319908786308797194474006510050505\\
11&296971760526924980280362522164784963560933267445709052664457436290717891100\\
12&392087644205509007971308171018096920842638096704083445751130532441358652580909\\
13&512584476313433348485635209205000203096612043013868979826304113479231724698039804\\
14&690402706744189590174338559167053706857472684991910919927535341782649968647659934961\\
15&973486733130064963717122815394865755105671644902923288744807296416314823573118318080956\\
16&1443620320932770094469128960169812465527730786977982175013602699212716970988526451967780769\\
17&2252322471932263048846263902780969720724082702623746692326411883383737891701105002873033615040\\
18&3693885607568822710108933818006069910654124178647755523042351854191467503115642868655945757393249\\
19&6360891118946817231543178165728090326860287788976903220920201346803034936457376719645501089451667804\\
20&11487774338954807720372320355973844448861968710939797941095159379112461406975031797328384022700910946765\\
21&21734870223008334492197915773663938720530831514510192502616003694589304348331278396493855673531175432254583\\
22&43035043306578737226745592294333293203410835313720222628276937945287552927709552977115657609417050695479784025\\
23&89082539702216186208520180657398071408456779082621043492610109878043261653419108557949483605163011614202446690924\\
24&192596450848292002554595197103558484035159984232905476096384614328215758546756634503343447533656675796918683441748604\\
25&434497593116213342300323139637022168642148198915272185981031693392402091438348991503921020016205677691321919200527862015\\
26&1021938339009823208266090002093993567096639080871145869244755614613998484515045800645869918761700998254164246320899538719569\\
27&2503748930574161663664779139310013432864082807050008317207079997396640660543579535448023588252112123234661956975323202050053551\\
28&6384559772964130214457526882859451522919509701741962945420462677179983645832590933533505577448915892757385864637207021318477835621\\
29&16931852517900876670581631786157954264323627991472961749779513286231304480261747749882504061727005466911523166026628371036530923744400\\
30&46664185539334816714854426094630940487075443909226084851689108275592339006751653996945798481642541961656936116212613392702195227690129305\\
31&133552899013576245547806008673896217063561515176338763482495385449776081949501858389161720672132569725006014869877969671395926713442751922623.
\end{array}
\]

Now assume \(N\ge32\), and set \(A:=D_N\). Proposition 30 with
\(j=22,23,24\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{23}\): the base case is
\[
D_{32}=96800425246141091510518408809597121
>41538374868278621028243970633760768=32^{23},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{23}
\]
for \(r\ge32\) once \(D_r\ge r^{23}\), because
\(r^{24}\ge(r+1)^{23}\) on this range. The coefficient of \(A^2\) is
\[
(N+23)\bigl((N+1)\cdots(N+22)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{23}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(68\), denominator
\(25852016738884976640000\), and minimum coefficient
\(51704033477769953279999\). Thus both derivatives are positive for
\(A\ge N^{23}\). At \(A=N^{23}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(91\), denominator
\[
29057685629025609672383529751884595200000000,
\]
and minimum coefficient
\[
29057685629025609672382405751156817592320000.
\]
Hence \(\Delta>0\). Since \(m_{N+22}^{(N)}\) and \(m_{N+23}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+22}^{(N)}\le\rho_{N+23}^{(N)}\). \(\square\)

## Proposition 53 (twenty-fifth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+23}^{(N)}\le\rho_{N+24}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+25}^{(N)}m_{N+23}^{(N)}-(m_{N+24}^{(N)})^2.
\]
For \(N=6,7,\ldots,32\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&99988556423933094809622988485252764317959145223199722286875\\
7&3602140834398025969293129379752595257680290439852360416148519136\\
8&29938382499226016919344031061980798518920420102670606906249669826911\\
9&103051316433965383815518405461915245054567452025979628920698145769118080\\
10&213979561455305408724613018182131249177619114072042900230752047443202580407\\
11&341512703473835796587412884330459786855377144253090839831797882560030601071900\\
12&487405615810799166943363677739437904942692927701875818288294888063853685892187691\\
13&679668451804618415007897251154925616259051285406606986177285452919549238518122542296\\
14&969449920389864192694635604427442559414717485414847346309898469435467687925920845261495\\
15&1442252751958608395259060806641390336394810249756053211714243995115236953152917692710322664\\
16&2251898386238294500323344162732298993652846216712951386722326961376947263069355192611465514311\\
17&3693763691293834814450687193768437071961892037305540922958004116827453012213426979064003460946000\\
18&6360858291496591923777033413191164637301512253631117388560174304277962917890781294274814146448704351\\
19&11487765989599508583395488495415593200121769863191947852014132550422128235744949803380881528472006144840\\
20&21734868202916018780280213825938421261651968998983683218846737229562559676633152149966861764816848616453457\\
21&43035042838870982876557809251658679314895049900217934311135404946956639939094356723911309481439117232374032006\\
22&89082539598065184810163755253469246914409202584735469161552858014468464971831640695900432279013307373709651616975\\
23&192596450825887792660574321391553646883154612075630502349650493954188115192611658681530481461440636514183830799303906\\
24&434497593111540162453516957201681145261097736414698516407098192957627716185157077922539887318520385222635456744164199165\\
25&1021938339008874914880743579047535536077272595870080005788825744181116304902553723241374327017266301602486423376674287720767\\
26&2503748930573973915064226468714403712780977861915053908687794958393834932758130611114188705837769454489373017108148217138862591\\
27&6384559772964093854448154447786704736130473284412348772294457902668504934494648856253714250422243966485983052073994551888200923699\\
28&16931852517900869767081879632109842931250730651518097911410957498973978135872886394361421976826847669006489128188279013600673641151729\\
29&46664185539334815427225662203675140670466215261073922258014842195177389996269178559838722140627111704725097486585047198560168722265566500\\
30&133552899013576245311445974353179083622495379330455698789730604843872537434918132909065128823855722135880974516282636091029552438775404854507\\
31&396652110070321449253715312819127210475238272442856107621604714034106772008513474998389142493351658115325357004725993492851186696144831261823146\\
32&1221688499016590063828948162691707302906328034956958898375660364653813407437692339402727307448997824483390360199377301736367220273202498937345496604.
\end{array}
\]

Now assume \(N\ge33\), and set \(A:=D_N\). Proposition 30 with
\(j=23,24,25\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{24}\): the base case is
\[
D_{33}=3194414033122656019847107490716704992
>2781855434090103443811378243892171521=33^{24},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{24}
\]
for \(r\ge33\) once \(D_r\ge r^{24}\), because
\(r^{25}\ge(r+1)^{24}\) on this range. The coefficient of \(A^2\) is
\[
(N+24)\bigl((N+1)\cdots(N+23)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{24}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(71\), denominator
\(620448401733239439360000\), and minimum coefficient
\(1240896803466478878719999\). Thus both derivatives are positive for
\(A\ge N^{24}\). At \(A=N^{24}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(95\), denominator
\[
16039842467222136539155708423040296550400000000,
\]
and minimum coefficient
\[
16039842467222136539155682571023557665423360000.
\]
Hence \(\Delta>0\). Since \(m_{N+23}^{(N)}\) and \(m_{N+24}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+23}^{(N)}\le\rho_{N+24}^{(N)}\). \(\square\)

## Proposition 54 (twenty-sixth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+24}^{(N)}\le\rho_{N+25}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+26}^{(N)}m_{N+24}^{(N)}-(m_{N+25}^{(N)})^2.
\]
For \(N=6,7,\ldots,34\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&42751231924207826133562424194438598813969941947317665054887250\\
7&2196617514531631701138397372091684568021587828097979545814022831004\\
8&23639137727326453080630820489930876465384200348551277260829827142657809\\
9&98224309015165421189173795994099816315627246213085043118084967134955191296\\
10&234036993736573782069416108978450269330827760386762582396891634790751239840253\\
11&413725070944298563820013280922485203622626847495177847657140978182933832018339100\\
12&638899543179629587143555470269445736606472652396107222477881505600433831658107106349\\
13&950278517785645230108005819655696737010307088694449564542821614456415719711282109055536\\
14&1434342781165073043546079990196175141799907596264953599933966230267104466731085702390877981\\
15&2248971028146687505016222073492431813574379807360248586090668374929096143526628633865047620464\\
16&3692775621969854227017548605729690554005275502327599249278066229178157461245094566319236911430781\\
17&6360550101436856039150444431199660792445671700376192533747817924495804637318372304575777810933693820\\
18&11487676208086272203957389450754429076541992241712997328260822050918852348058990648048953829276923604909\\
19&21734843560221517832696398455327056735464113293634335657452085178436587806886288914586926610202638520829360\\
20&43035036419609984757438941944318587072145084802090624028619012584068078960283113194497896807034453637387159244\\
21&89082538001291162480346344054126479916655148823477722732839612370531552210376191909775435663473717287081550017292\\
22&192596450444616350955579811528211124462322213974053335479385072034820084237787302669681580369169410759652035238258825\\
23&434497593023757357408308744795150346564793906298946817790308669848297780593156395157960629765271874189801535366297344911\\
24&1021938338989310786410657696662579263653824228943345985258958725264176400602130320147730448915021574340374711486115082878541\\
25&2503748930569738813934339475098249645360160701091896048418772290082150172950903556877189110119736960033962321826552246905181055\\
26&6384559772963200721545886998694414775247233942210847616822177963270638629255268126116445534973650711711828537807609715510646485021\\
27&16931852517900685790310998815808842060847059435247621536190801045379038129957342710550215006146971005363233506827245420984682110645351\\
28&46664185539334778323146259014249870205622883366305952573721164045928600292259462476723438928477164984133481985320669638710922538464967964\\
29&133552899013576237969699132412631639590216541167550770306511249080432101358043467128609770966711309125121715527583261445018353549928423885055\\
30&396652110070321447825778741359192004230623449878030155320912323307639622559061919264729957049388384566162662876007206953045150033228477144835913\\
31&1221688499016590063555496297597933271818485051535672185314894381399864140006106788866350363302838203733488083479561315836619825385655919061937256412\\
32&3899629688860955483713993378320472638808520578349960702866270217115752044672025238677742309115473719685653651559007425205059405463812459976984041893136\\
33&12892175751374318829323875973169061195705086418808299772606464943718808025325802189777437765159971504928861218702406212959668379202807715037083361076495356\\
34&44117025421202919033978262393500063273095697095550385231345640126217511996697243227396757068105463697597058167340507511701434089547033185649025297170807844490.
\end{array}
\]

Now assume \(N\ge35\), and set \(A:=D_N\). Proposition 30 with
\(j=24,25,26\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{25}\): the base case is
\[
D_{35}=3801352699415960663618057913952878940514
>399669593472470313551127910614013671875=35^{25},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{25}
\]
for \(r\ge35\) once \(D_r\ge r^{25}\), because
\(r^{26}\ge(r+1)^{25}\) on this range. The coefficient of \(A^2\) is
\[
(N+25)\bigl((N+1)\cdots(N+24)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{25}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(74\), denominator
\(15511210043330985984000000\), and minimum coefficient
\(31022420086661971967999999\). Thus both derivatives are positive for
\(A\ge N^{25}\). At \(A=N^{25}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(99\), denominator
\[
9623905480333281923493425053824177930240000000000,
\]
and minimum coefficient
\[
9623905480333281923493424433375776197000560640000.
\]
Hence \(\Delta>0\). Since \(m_{N+24}^{(N)}\) and \(m_{N+25}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+24}^{(N)}\le\rho_{N+25}^{(N)}\). \(\square\)

## Proposition 55 (twenty-seventh transition ratio)

For every \(N\ge6\),
\[
\rho_{N+25}^{(N)}\le\rho_{N+26}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+27}^{(N)}m_{N+25}^{(N)}-(m_{N+26}^{(N)})^2.
\]
For \(N=6,7,\ldots,35\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&18822512533488054546721505757402039139997960284080293087763859775\\
7&1386996928158242193420090063523408043391859222539331804745328761705196\\
8&19425001300285288381858404035951634736043391511005307311193041377763323263\\
9&97866818523759406923618647152746929014296480736430362170792486505737104203520\\
10&268507009872634595296387281396111917059161825898204355171288168828648493489449535\\
11&526957060190918003098307812315200567988087993713643647052771657601407297480166328400\\
12&881528050717157795224593516129379004991213752331978387484857360114929246588109520495151\\
13&1398742270292083918263267087455264096820161719807841037473374360801530450973971146890598088\\
14&2232884653918051259211374028106560226021875374093922852645618639607908318633248764683455565255\\
15&3686283536193215958621692773412900539949306637969966192875886083102090205282074791477087317713656\\
16&6358169097019563359977261468669203563332920118238232745521881091524000177531530140786364926942871271\\
17&11486871760693181569997556220233669900548059547537238245552722986206914403274229830904412505339650843916\\
18&21734590415838668515543592077217548061304019874202088997379743421077253097750707063476238579934279251763063\\
19&43034961552560587542658587915603301827906524400388358311534100192340168633100740558453454599768033100439084516\\
20&89082517034725958821310215381036375260820044647154975107963638265675888712113333586798835359364163098917693999106\\
21&192596444849274884650067288544421123204638846636768815390789059921305191178182452889939111947488360369081078101407994\\
22&434497591593122359060604171523853056579969428388918610650889619318853310998533874542513398779087815152522311903329462711\\
23&1021938338637227691467801104734881407665299798882768581756172845439992558914694155064484153938973373141370763818593138517979\\
24&2503748930486001660601337440273731250847007215140083927043276079628102615178843690428102357873697021779190960155547331793900199\\
25&6384559772943886869596339038223437623363335645955626894666105445292435916970208704780821780940540242013167957116077414541970010285\\
26&16931852517896352343401548254656140028477281517476615525829602223407182566735622173386870422589397815944202282797803142586562837605751\\
27&46664185539333829920856352856980835675421537980884457607941056210392770827662148408294855823355083549675935653764027545383542097558557700\\
28&133552899013576035016974414920441958669224371460306451971119146657891047346345105551598139801169562113490546923214846158607742369185391832810\\
29&396652110070321405268277824644906537539978756674019372577916623294641720755812939523541444463600863905248923292585997820179670570911453891711664\\
30&1221688499016590054794069326877474407759563530429622720421289268257706916511795328227826652831319714995670161498998955552986446580765136106092095495\\
31&3899629688860955481940026889298607063773663144282850819887319414906698808936829876422380098302238641679250448181550124725675934722887338838768106049014\\
32&12892175751374318828970067057743834572827078728359076775773725499234094601973110947281947017116059406761317093436961804237860165212590069934057090163680788\\
33&44117025421202919033908654000781510774316285286813779133085499938447162361314388939529213893463897101354023007532560460394595693186743374382110554838652306992\\
34&156174269991058333380275623243784568391795317101719292189894707984729702322195560546844492938608076722886837035374230878179774601354970109672256133636805240039009\\
35&571597828167273500171856840708396057867388278603196492485444411823361802289417516574828719012380644842932669019633188715158618722033424173622088189508118984853535700.
\end{array}
\]

Now assume \(N\ge36\), and set \(A:=D_N\). Proposition 30 with
\(j=25,26,27\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{26}\): the base case is
\[
D_{36}=136848697178974583890250084902303641858505
>29098125988731506183153025616435306561536=36^{26},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{26}
\]
for \(r\ge36\) once \(D_r\ge r^{26}\), because
\(r^{27}\ge(r+1)^{26}\) on this range. The coefficient of \(A^2\) is
\[
(N+26)\bigl((N+1)\cdots(N+25)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{26}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(77\), denominator
\(403291461126605635584000000\), and minimum coefficient
\(806582922253211271167999999\). Thus both derivatives are positive for
\(A\ge N^{26}\). At \(A=N^{26}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(103\), denominator
\[
6255538562216633250270726284985715654656000000000000,
\]
and minimum coefficient
\[
6255538562216633250270726269474505611325014016000000.
\]
Hence \(\Delta>0\). Since \(m_{N+25}^{(N)}\) and \(m_{N+26}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+25}^{(N)}\le\rho_{N+26}^{(N)}\). \(\square\)

## Proposition 56 (twenty-eighth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+26}^{(N)}\le\rho_{N+27}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+28}^{(N)}m_{N+26}^{(N)}-(m_{N+27}^{(N)})^2.
\]
For \(N=6,7,\ldots,36\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&8520047457098094903429543443701071845221196482395128501935584380625\\
7&905189013128731228908333624925777973385584981652604194752757673362246396\\
8&16579124543422310709091530163860751737242618882033095083338003136727837713461\\
9&101721896415569276052699976048311145557264422987237677581107903171897244858857920\\
10&322489251994535590282523816910361628558320128570700467727598593532230959912435195945\\
11&704340276631983214573288588346642659771369488264902455284560669878127904273194186888000\\
12&1278087370446724992664584382623496910918955244751049402457989879183006987710656630154858625\\
13&2164190908263157473311299582455535812237309433355114208095409524934908679209001868380975927324\\
14&3652346149621104354818019442141968119919176774322228127413243673378039539942292855601870257040265\\
15&6343258063134231537413194809024369958095947177116250253766554112006322864384470970465020419631607164\\
16&11480939170219977632542052660339623320641998681872290017155607669721851952642802731985276796498063913425\\
17&21732422675805153891223426543192872159161674117960452310695694226246374318376746912945067470511386316861468\\
18&43034225814473565259118337738401267298702297975065806964421800126099296863489370587247053383633232573377053521\\
19&89082282917308548986135650325136287998376233380429590636305581676437550093354998095515498350990446577976054769071\\
20&192596374460991661547223045691442038952290313772426026560087610660707575224384897812555332619526322035273157287740664\\
21&434497571467298242034463265224142350143236671696181641927177765511360757296242487312401204411279595423880786900222399095\\
22&1021938333134190694786766666048980301090414675648852576467847931084008936270561802744390479142731939401499122671953701732585\\
23&2503748929040175502270093076428127911936482863778558147403132032290791623242098618299144968245283438107591097315132301330526775\\
24&6384559772577375533184222351562861226700176020871750647492815289506608733523615812818440926518080019741533520196365459973298559601\\
25&16931852517806385775101393054686938556525272671473774677736964266525524413782599635076132154104476016280137929006313726345903738635639\\
26&46664185539312377713175623416748809197830626473502196868002573126596036622014551607385397371254146025969772075828901006590498103700988545\\
27&133552899013571052201735978775657702116728243835341524056326035879307856030352460322394285568256272927009553757359366493396188202509338992976\\
28&396652110070320275033110899967068161308625730130114482874866256362273332000555424207223343794646381877955850017351036395529443461664619091065696\\
29&1221688499016589803884344193662065963897113505155856576444052984644055959361702624465716451538378807693195750278441118172772924625613936343175320752\\
30&3899629688860955427315903523526577737455852659396554025615041522901890751317740620768933637529507430341701955195555416926572403089334922995057932046905\\
31&12892175751374318817287264577339117847339200514597923445529064132396664940563231430762305550438247437431484065884191368883660441116118816710082439518281015\\
32&44117025421202919031449944918113895820111552486425769770837457431238816784606273221551291322930844412975727122861270713813510972938146137007426669157618237804\\
33&156174269991058333379765704614464320578153317674398027328059347603132324784678669863485890227960821243542005750056514812149342467468401834277612027566729658322944\\
34&571597828167273500171752487386054014327287724421351867408205318762488534492320381643491151979505216428310531006460194055694840238939915660024937501193555381389501297\\
35&2161782986128628377649951024434157365426518001690305306132975334224128470917717711525770944496999694436254218499674957366037472286068065270421194132130526787364474575159\\
36&8443924343818422443100788770874135565850227412002531583348802002451192884495963050231203160906865896541550471122340789939533410458769876009996072411626183069228126697458720.
\end{array}
\]

Now assume \(N\ge37\), and set \(A:=D_N\). Proposition 30 with
\(j=26,27,28\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{27}\): the base case is
\[
D_{37}=5063401795622059603939253141385234748764684
>2195060778453276448585190172728648509121533=37^{27},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{27}
\]
for \(r\ge37\) once \(D_r\ge r^{27}\), because
\(r^{28}\ge(r+1)^{27}\) on this range. The coefficient of \(A^2\) is
\[
(N+27)\bigl((N+1)\cdots(N+26)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{27}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(80\), denominator
\(10888869450418352160768000000\), and minimum coefficient
\(21777738900836704321535999999\). Thus both derivatives are positive for
\(A\ge N^{27}\). At \(A=N^{27}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(107\), denominator
\[
4391388070676076541690049852059972389568512000000000000,
\]
and minimum coefficient
\[
4391388070676076541690049851656680928441906364416000000.
\]
Hence \(\Delta>0\). Since \(m_{N+26}^{(N)}\) and \(m_{N+27}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+26}^{(N)}\le\rho_{N+27}^{(N)}\). \(\square\)

## Proposition 57 (twenty-ninth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+27}^{(N)}\le\rho_{N+28}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+29}^{(N)}m_{N+27}^{(N)}-(m_{N+28}^{(N)})^2.
\]
For \(N=6,7,\ldots,37\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&3959156656862021238454578298420304511758448541973133393464265624651375\\
7&609567669826793475742860615308288987432886838435318685191419281973837331476\\
8&14670362890250756962664722022735578529707394053833955627985598806439365473452839\\
9&110084155898695194800161801205576393197955039721603478987512552237698439272942763936\\
10&404706510537808637651779480956876814855394991113096666645350748077364842061345973422631\\
11&986182361752417049643083100624341733714153584638212267489494024341216016116174424631863520\\
12&1943997797241645895222032648973398715968775384962691490235723297044678122959516832338098705079\\
13&3514697613234349169865546618887381329611623454877878627161590345884161204425036040228223829503756\\
14&6269021104223235721416282720014688400287680734344801434341549185124102083396781967877108096625354199\\
15&11445484627906256613524759907694016900096184090803477552240154358510278947991574585678631825254645142220\\
16&21717143430007749990327550093318867258182246131216694379519348927098860439367326039068347419366249271858343\\
17&43028196847738916945279239782698195482889735303851004768359799089433335100326557603715604033841991936387238264\\
18&89080079084526215424916609145990295936284554211821987283600851069029068241492747719324829360945524503530006419079\\
19&192595620978278257776599334966958779642727933005657107463806114494876648058516941232499834620824272102097045356133841\\
20&434497328589585936761149720311659398201383168324966093248340145565960896624089437591923241595622751001542936867634775395\\
21&1021938258826520492439482968514622667768104128660227386932694281111164215674281407254673402608234400098805216789788383679829\\
22&2503748907338763667513100827230749691653058837195199647292343200992718412438130345306261175628140122749633692108295726634832631\\
23&6384559766497738274955981604379571248155290900804298120777753742780097326792873708675818266964433236183011681322257245957399816025\\
24&16931852516165626446321049192303084902330652722536315313170058949572241584872657294506959345934746418189810534935867268339183096485399\\
25&46664185538884226412931590500884088279768878454851635505124289343441770470779144221118683858185347704795479240269012645256622690120898527\\
26&133552899013462672688054465715579842873432479825686133368522262742260597318153122121136986956725797724627425036951284208307029891596373759623\\
27&396652110070293585378898595081344158495008811787580622967418808589110713511073216848307260352871012567087018656323351025436030355420917048043760\\
28&1221688499016583393358162605675785101731650615315127095988653513492739738959159953918171105325660459119489620295177000652765765273729981636864363929\\
29&3899629688860953922131507081679358713698240265966498528961516882670344499077512739299098267453718942434456799160055628935841379991546324529297219797816\\
30&12892175751374318471093616074430867268913574657053351474814306731879208941517760762790607437113341071523368868981351696838966705278577094540513307841955431\\
31&44117025421202918953308514559661146811973182968877994929724924864895508947166888377002576154112195098522078568136563793662399214852912380548978030863285732007\\
32&156174269991058333362427709440587732260074255934393762173062306724547254629429818208418415991583627463405161173969749185719706461129705354512744264798442213671934\\
33&571597828167273500167965215285874113949379052224929855292480831756309818933956718529995302089901439450329202703635868326595476682236432071718709815100694115635287456\\
34&2161782986128628377649135457394817803891343438759649092976788583709993279263178766948561438310897559058974863802024566416708170330396621000100573946862043823510708333611\\
35&8443924343818422443100615413935271926595006359370405099731780040183011364056628242419604672868564239462440566822000985052795022441534543574817424084188292171361624059499396\\
36&34045902954275879290582360451819177053273633229119552874337976009332085184304735361386655506130308656951793692151066120710372093905423805740847550691020107807952857964809930647\\
37&141630956289787657848822766776531444650886183043310278486580830685559648939184899011889021985704534603037806628031029150212098513903781376717874395799309900519572601476360663971236.
\end{array}
\]

Now assume \(N\ge38\), and set \(A:=D_N\). Proposition 30 with
\(j=27,28,29\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{28}\): the base case is
\[
D_{38}=192409268233638264949691619372638920453057993
>171372331490336173605981741257818705594679296=38^{28},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{28}
\]
for \(r\ge38\) once \(D_r\ge r^{28}\), because
\(r^{29}\ge(r+1)^{28}\) on this range. The coefficient of \(A^2\) is
\[
(N+28)\bigl((N+1)\cdots(N+27)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{28}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(83\), denominator
\(304888344611713860501504000000\), and minimum coefficient
\(609776689223427721003007999999\). Thus both derivatives are positive for
\(A\ge N^{28}\). At \(A=N^{28}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(111\), denominator
\[
3319889381431113865517677688157339126513795072000000000000,
\]
and minimum coefficient
\[
3319889381431113865517677688146450257063376719839232000000.
\]
Hence \(\Delta>0\). Since \(m_{N+27}^{(N)}\) and \(m_{N+28}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+27}^{(N)}\le\rho_{N+28}^{(N)}\). \(\square\)

## Proposition 58 (thirtieth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+28}^{(N)}\le\rho_{N+29}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+30}^{(N)}m_{N+28}^{(N)}-(m_{N+29}^{(N)})^2.
\]
For \(N=6,7,\ldots,39\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&1886124885775560466352433001144220041990743958429260767302556861303410725\\
7&422915988503202581669643114745088055552791756606301592553458239079368937254684\\
8&13435793665997150715721955191697363959204511543279991641195301078998313202720720401\\
9&123819469002281333923065644213842721947391607137258781557423509001937337094598767138800\\
10&529726792907340145395117275058263895326197168677164538854432510747899493774043872779259397\\
11&1443971967295404967202145901681259964713021311145681443644857989981387609404465492024363807984\\
12&3097139212187978136705292625187256774402303169495654065865779475484511550232110324122972664132121\\
13&5982766790030760679165432028406426603840822315062075087839278546706378441353528525345612690012210044\\
14&11277192230029844979242885071695263535188419903617782153895055886739720602376427704046773091858638772325\\
15&21629907297809186768952807209199399960815718724399674888469076026732137003157542531178425190269503331880732\\
16&42987535028683689194605516707029142001918122488568664025897008516236811806359019284105135366596762520243912381\\
17&89062777332263356913404248145500892144757650946153658084332635853459608788268189911875187990937652016813609366576\\
18&192588818782885551923567709397286874923076551069348269398928918293838847490212031835752033747746223482556018793209126\\
19&434494833125360116538462743755419889686318480189981911736674457925594862286228529016718114278530948408383831998279805975\\
20&1021937397522117066783616430157462811836347391179768514833150033894724744828117488505901379729998645161479216167502448937137\\
21&2503748625722319344610089564761558427812454464293146115073159394100020707959226316996329079946091817934134777462712138804537775\\
22&6384559678757349713847493229412133043750717166757076101619871538868847510944246077046124697873516056719434618271979353403828441069\\
23&16931852489985894972468008381320025138222642749568199417467820031825280597014876452450100832958152985964858837555500099756346276094775\\
24&46664185531370665285894219798002349125767125162211607770500115295627377388821319791824920723281201858201814426773797198868066386111681005\\
25&133552899011380621996129787681428456863873243208045197371003810156390684579245445859872339912465205934190412282537638733899069436169222254231\\
26&396652110069734663699202340404222077244819917275564311449508641789399430173739417629994196933885528254463390814534156776034331213700123564081806\\
27&1221688499016437610732364407992174726044889562135840477994648487614841137051481580837370328957856804852828505840756237257647258658875365903623023600\\
28&3899629688860916880160973401679624278999674886980182661785524545948992700872402650415185043867522931176553149845399555117078390229815733048103861363341\\
29&12892175751374309280735121258017936885153675085090956224729597549513035699503918484468025891196915921221957718486761672253389539125372017604821102346070620\\
30&44117025421202916722138333432951981975477788153021678462578449870595708319141730910356177993057160152865318684937664116031404089947076636009596861349457745342\\
31&156174269991058332831401488110355923344035193395915661520227838859489797355864381918810656910127323883185978342103711017431935995352309122904230131032077782099719\\
32&571597828167273500043850161137708302469118318848414199413147360063750443307842689628311015030329378669185367847529767145917942138019675774674794784772217423200250256\\
33&2161782986128628377620603520230508525579687655510408688101655363550262103844443880570787944283332618465242501682138780188638493696862112282106415980303314502115211054279\\
34&8443924343818422443094155174438167141868935963100270261464864811854183193489947673767887942282557017436092281929130529163941628351245289566108415209168246949230296814204745\\
35&34045902954275879290580917887981813276606289975959300523456857288898122175931026812878033738784149342633667331639360757736584621706338975243334031532587768127432487584253849024\\
36&141630956289787657848822448722203862531851422072999628159035955918808591830297032102765823879926861624733018026570437751107832761323104879046755458998658950365265094955879813750649\\
37&607596802483189052171449631746159506286065959649598681006775137865720856282599748643424200114407277936665798509404269643666445802174773810292433871610239302575825348241746445984194044\\
38&2686793060580661988702150569843058070501462623946806896746894352017741323465032371652959431732089176188271626080572780177311906743708706016247044354100801517639760227339419674298839461925\\
39&12241029184005496020526998062514467511721167950650307087743974726891860738840237728774244135844320516446265674152513671854182431266804563634604100102053704705548628389158955341340372036573680.
\end{array}
\]

Now assume \(N\ge40\), and set \(A:=D_N\). Proposition 30 with
\(j=28,29,30\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{29}\): the base case is
\[
D_{40}=300158458444475693321518926221316715906770469041
>28823037615171174400000000000000000000000000000=40^{29},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{29}
\]
for \(r\ge40\) once \(D_r\ge r^{29}\), because
\(r^{30}\ge(r+1)^{29}\) on this range. The coefficient of \(A^2\) is
\[
(N+29)\bigl((N+1)\cdots(N+28)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{29}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(86\), denominator
\(8841761993739701954543616000000\), and minimum coefficient
\(17683523987479403909087231999999\). Thus both derivatives are positive
for \(A\ge N^{29}\). At \(A=N^{29}\), after substituting \(N=s+6\), the two
parity numerators for \(\Delta\) have positive coefficients; each has
degree \(115\), denominator
\[
2695750177722064458800354282783759370729201598464000000000000,
\]
and minimum coefficient
\[
2695750177722064458800354282783454482384589884603498496000000.
\]
Hence \(\Delta>0\). Since \(m_{N+28}^{(N)}\) and \(m_{N+29}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+28}^{(N)}\le\rho_{N+29}^{(N)}\). \(\square\)

## Proposition 59 (thirty-first transition ratio)

For every \(N\ge6\),
\[
\rho_{N+29}^{(N)}\le\rho_{N+30}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+31}^{(N)}m_{N+29}^{(N)}-(m_{N+30}^{(N)})^2.
\]
For \(N=6,7,\ldots,40\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&920027928161318416918832570046508702822537621019331842837606970037330385975\\
7&301866901145569002619982014464386386119645364413621190521710617840612795521085596\\
8&12715842339970004119714465659309541555541854673934330748271603112617227363769784108375\\
9&144503652485988988040867552064280945550682177291093131754116915331923691420748179990311000\\
10&721961770138223500972529542885137292160205012843107557745180065428941325659379035595154106559\\
11&2207397544408925375602215607248575112188672228985679852437790007915704643541592938207992605980992\\
12&5160618445226671833406283864758433545132568333597985940867372003068407846939279546338994274940946559\\
13&10659759169462506503307510750495534917475313320734675548623465891612386370614924643195545921757570965552\\
14&21234731733704322277606282584166708404393872701310585570934144214975890415590204975007459517767705464499919\\
15&42765512512959818697433347211528068609525118578550720024904632454324064841920016722868874987985203816978550464\\
16&88951001765319860692795642678833434174928086408809199422004051556986020371989835199226001143951590373527704413095\\
17&192537600375152728521055773787271093064763278624034237616436937970500569552996838737635055441131931330313054228703484\\
18&434473203984127682963448356637108710287882723703489285684514960340089164298301736345019653366926152737372521808529246431\\
19&1021928894013758766073549399757172814419902587481768744548928866783035518918725248669075347367762568038376139382245682158317\\
20&2503745486936785410894322388835613675230381304766271729712508539600842662777148263449254113641127610316763717230744028552367719\\
21&6384558583279577185132736251140857626014505133010764674063943027156216559929711276343133330385514565143723114501721056767862774375\\
22&16931852126292741149163391312719765930098273978705223212991357271099993842262907077830068021876528182683910270041843071883168504905631\\
23&46664185415919956353568269631973766267829668015866058864492010779144058832331093227178081661292284619949017418811414044756434889348194299\\
24&133552898976181888668371639519702622155444927863375839560840140160502524685465642358439810056267864311579535832282090470224528227830750424207\\
25&396652110059387583164207428243346180436332383315429165350782103446745745280412735683994608421122882368943184450656611702166730350619106718445402\\
26&1221688499013494852018749168967108337467279121252429708477631099949401326156895554349467429251372084981799993794599751991093412801002180464565561419\\
27&3899629688860104703424967109710877434128705955906492613340298801871510606498169979266824433928860932455566567982161569939541157018980956242820338823320\\
28&12892175751374091172332366123815227100212570714283592983353016674908853423612299679234177547137995006755640817654556920990717044123645702090633602567968231\\
29&44117025421202859592379118718041355609601100640643302792386417503181593783086034123245099197396042001764303075276566857964855881632350868219351452371306234326\\
30&156174269991058318204392901069343991314166972196672654012292948716722735989709687756106616533536652096419951872134695449233977224963330418986026235602139204493761\\
31&571597828167273496376148626998398909540467885015896824080530625420176178215774787842005775236906299891784262169039914735455106195255155058592264538933335508890269870\\
32&2161782986128628376718327786788701319379561500843901895692347354684965218711368063026610475823860493559907706623159816629438148028727219080877403818283731272627278997004\\
33&8443924343818422442876041859640909922468561755308715284297251146415825319075379928542723113032409282317057125343683687027610672332702065844766550971068046471123188483970530\\
34&34045902954275879290529031755547613950232992581388991180529862819317949919775545112854006392085654146045773060680254489183529391224278086800163495488463561257472479169840104039\\
35&141630956289787657848810286293770606924888865531830392429339358727215545050645400921123054037381051910228785863266983808535389738028267107726593331637333547048126105239198427361984\\
36&607596802483189052171446819106799782644559004630880365925585905545674877120968789460005562586633196533026825886240380862940932064728443201376766042623712794824340191132923326011304551\\
37&2686793060580661988702149927430577359258746905194469885213798292096794613526216450527662541195756465600292940143064827925165696670910685283453150997496331273191035656815994675241512612276\\
38&12241029184005496020526997917449481981289958746550449770763919499893692008019169306199896431560817969417538902928665203130021775608626298718197853490516063415589465636427026760852611006942375\\
39&57434908931353787328312674891493468951377465886105560137427854029295284614540870435057170640466353817747210111343614880241702616167665224312601628119876023458949653699119930006569419991463914448\\
40&277410610138438792795750219878492256391025105975054184466126972582644018263248867082831987384791459010422009377221421306341224077966867930009765456091112934591285932249127419483731934214974836005479.
\end{array}
\]

Now assume \(N\ge41\), and set \(A:=D_N\). Proposition 30 with
\(j=29,30,31\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{30}\): the base case is
\[
D_{41}=12306496796223503426182275975073985352177589230680
>2418330769289520463963038742566759257517431737201=41^{30},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{30}
\]
for \(r\ge41\) once \(D_r\ge r^{30}\), because
\(r^{31}\ge(r+1)^{30}\) on this range. The coefficient of \(A^2\) is
\[
(N+30)\bigl((N+1)\cdots(N+29)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{30}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(89\), denominator
\[
265252859812191058636308480000000,
\]
and minimum coefficient
\[
530505719624382117272616959999999.
\]
Thus both derivatives are positive for \(A\ge N^{30}\). At \(A=N^{30}\),
after substituting \(N=s+6\), the two parity numerators for \(\Delta\) have
positive coefficients; each has degree \(119\), denominator
\[
2345302654618196079156308226021870652534405390663680000000000000,
\]
and minimum coefficient
\[
2345302654618196079156308226021861810772411650961725456384000000.
\]
Hence \(\Delta>0\). Since \(m_{N+29}^{(N)}\) and \(m_{N+30}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+29}^{(N)}\le\rho_{N+30}^{(N)}\). \(\square\)

## Proposition 60 (thirty-second transition ratio)

For every \(N\ge6\),
\[
\rho_{N+30}^{(N)}\le\rho_{N+31}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+32}^{(N)}m_{N+30}^{(N)}-(m_{N+31}^{(N)})^2.
\]
For \(N=6,7,\ldots,41\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&458976551728687329318428163850912841326539954024699228360769641792341651578225\\
7&221375744920179019039317075368555904801141338655528610927129724799507731540211600732\\
8&12417934672363808030736975500012127530558926193933202252014956011102027807963865216856421\\
9&174707897415190922768284822092534902935623697445125218903186923880994210260085324505733505564\\
10&1022894051575798577464813195600388631392457720677222576583139235879131510109743757702760103993521\\
11&3517615550443289314340907645037205503713905256431369328707637191967153373371506675438440452254548736\\
12&8980307326081280825112620507175529016984490139987088537055203122815424607520776304279258312380283662877\\
13&19854296873941142520336314044721839441227989356335467042609423210555652369998216268682108103960149795953152\\
14&41804842217507519917768885529919948067802778128463001421380541356188352180714373662628240835408500409858664985\\
15&88366767925321263294557263476210713699337148761735483488362656611497525469297939564572741094356618158824501767936\\
16&192220332958274789487898993482122898933007084650472773726035944948889935643920626312382500529355668441749723939706625\\
17&434316845384719509712601447894677403408911353835593480857661960778258202965350174260954471595264491640982881139292516703\\
18&1021858059360692877219686887871400688773697957154634011566371154375254820856123194941910102328042176098863991539976512835769\\
19&2503715678834601913745375984178222128264110042797612087998806423484812909370744870521590705620992171223102007010465392116225279\\
20&6384546830359128661429999924978097138478770842231423091514441857285355113001730221282162219389637336221982042180119977487221115381\\
21&16931847752736269963545192784203704506753386996905389706775223314884349921058153444292212498073346888991380799427666586110316570449375\\
22&46664183870384086328831037603489807886075640441228310723572095437524629362549589688769901037818703158607314128105361257039328785268526249\\
23&133552898454776969044465585521276161536437495650620392745408535101964309606016338015701499332220686407431492162529653084504935663588200494847\\
24&396652109890689753553970181634928545288102656212709245727200474133619760732683205552563909187887904880598935005311131488360454696619315421816449\\
25&1221688498960939779906969230145191872209989516051253936287270487990102396659599022605025659004746758366128406456421630807221563388053679969963711900\\
26&3899629688844284503434363611532785821513818798033217417217275764240408465365192327513981912103271397935374516069389557193160345435285546646561607814401\\
27&12892175751369475382643827758785006756858950207384806895656880176414739245238664702970122720596141828576848111752948186072214468734244409526867655845731484\\
28&44117025421201550676792064167496518697660116656102416868488543274147890452767182211843569434995388885924279104365899534538936179104797917059797331170125061265\\
29&156174269991057956564692338191388471719647762712111721301649656467151324643989972589728068749940982486576406449810125229687457055409998370581675796672627057870095\\
30&571597828167273398810485153319455645819803424930906351474840585997150037281160257113745874349965268431747129064544275504583860838528265742287749122235360340661239145\\
31&2161782986128628350964777181854206639474945899266668862983748665261770618055108364456339421945848722125823591786345872711872510168937866347217030846863200482876959551500\\
32&8443924343818422436212936565131232971302508443096671434497493224329390476553234131621310922895309785521182035059374439236984815352941222789202362726891561586423661541264868\\
33&34045902954275879288836551085106293534888416140390370519465375193325871307387578912487824633468228379740757345780330835172750989482402949970658336169628814392781067213762715964\\
34&141630956289787657848387595544972157988111713288978927942953536291188329596951200476719715008359827666969636459650046704224934255613957712055209531985731412415310874783187479820441\\
35&607596802483189052171342884033035444761232540058884976404641706507637057001569380100806441996017982646857682753441650206758142357365591889700466552589518846124936397823547012651576256\\
36&2686793060580661988702124734398752260112045312336596333748189705777886210941025680956082071055863361939881157696036147020840685069942548841850437095759027079035968943569210539579902323549\\
37&12241029184005496020526991890856865755506618577374913687982062834908229345936744041650154126693172524817174559113460236879136528240367291593555802720044733771869841159869460580139573763586780\\
38&57434908931353787328312673467228466656149361955280037775325980236010078144491124302116467923201001990702067997114554396339257446137475038435059021206939976206164081380356579346862880646061505825\\
39&277410610138438792795750219545636411489289189261289872889124180270825220033871280340582039853325824272153329145461413794516411414551189385914495209726145555022804498975384744315423853156852258303856\\
40&1378730732388040800194878592753884858091953940310562123070688129521650568994125332332728604457436112782745469207250822293962492293953410370409437922600398070775484851566636183620972851904915346059787141\\
41&7048071503967664570596219366541660495311129748959929078166599386966841734163381858936839512115058758427505411624624850815276015271043474589940791515304068088950729749499254073421641226706317689889605642684.
\end{array}
\]

Now assume \(N\ge42\), and set \(A:=D_N\). Proposition 30 with
\(j=30,31,32\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{31}\): the base case is
\[
D_{42}=516872865441387143899655590953107384791458747688561
>209280132851277227501988209234026802122409700753408=42^{31},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{31}
\]
for \(r\ge42\) once \(D_r\ge r^{31}\), because
\(r^{32}\ge(r+1)^{31}\) on this range. The coefficient of \(A^2\) is
\[
(N+31)\bigl((N+1)\cdots(N+30)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{31}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(92\), denominator
\[
8222838654177922817725562880000000,
\]
and minimum coefficient
\[
16445677308355845635451125759999999.
\]
Thus both derivatives are positive for \(A\ge N^{31}\). At \(A=N^{31}\),
after substituting \(N=s+6\), the two parity numerators for \(\Delta\) have
positive coefficients; each has degree \(123\), denominator
\[
2181131468794922353615366650200339706856997013317222400000000000000,
\]
and minimum coefficient
\[
2181131468794922353615366650200339441604137201126163763691520000000.
\]
Hence \(\Delta>0\). Since \(m_{N+30}^{(N)}\) and \(m_{N+31}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+30}^{(N)}\le\rho_{N+31}^{(N)}\). \(\square\)

## Proposition 61 (thirty-third transition ratio)

For every \(N\ge6\),
\[
\rho_{N+31}^{(N)}\le\rho_{N+32}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+33}^{(N)}m_{N+31}^{(N)}-(m_{N+32}^{(N)})^2.
\]
For \(N=6,7,\ldots,42\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&233922223255560929189339214432068057029253279862469750985534430724301160747120875\\
7&166594899844930038646286645793578571299395639108289699560774306477426661311749266408592\\
8&12496311816540681963530557045817434336297182672367075831887708448532648150962285051853851485\\
9&218499663038413555148650065780554609514990867425127576411286635219196719951214782985888782451592\\
10&1504328812591996056476689421419875895867209700895531512367148109602526412064583411759576140730185999\\
11&5834703075852637865448016129009181565583471747219995813962595946062963811393957776737717504454611277120\\
12&16297695402178327839755758667378087187880309430179641375937608566292919336508514536716972662422559895270617\\
13&38607707524356286410604046544989178970384454824734550866427913042666264001418529683191212389724085810072003584\\
14&85950226299097604366811921803194549384796966432247296169274677208381615340659532748001793208240842724451725604495\\
15&190631505369885443737160117458269048917176112592450513342254821697968525734653308924045949661537764951188805249475584\\
16&433387317006534248998575157019956793575991811140269257243698321855433023705842988713868413151745188508031283821165283375\\
17&1021365970834261929435592992089660782463268863843140910110006315551741187143924360070886081895932627384322565037342841176382\\
18&2503476810396661415882342286912647918364251472181447137968848136730505406486917522223753033667780586509591205962773742196616375\\
19&6384439363964186351241577271249730627417610594188177969882159134580132012024438040893648626259983937012574550296417136148837458396\\
20&16931802541837431511146358634799622841390520660741517415146090812151755269590307252112151853178660822147875381285568466647713157067619\\
21&46664165951725071614464175551156324043476999059384867591025215314969065587905721578162954868891355548475321574753102248307423744004596250\\
22&133552891721844480054553572923478130204659248176547165374260517476936238299894151631779485279045565876885535211125735191485070928682536178543\\
23&396652107479167153168618343446230970540067834028386794490041995984709423668579107807825161057059635194420996072720192956389906024725797432587766\\
24&1221688498133758715612535815763794081737733601848507046318524835043949358628138730513720108947244220763250237410173204804227533136029932613528652663\\
25&3899629688571446421412109580605869336174434140802729782033103807934855348883094182138635708140777121630266093344030952849895692936294629351235495718844\\
26&12892175751282627103487535378510847809306865235282520700084916817820416600398921181369364160533050994978069417333855649321100175697267933947646160867510355\\
27&44117025421174787062840344994707251135313506716451975722952927770368357911465765518812283945724781192872957150814688365737869431525416397052332329703795394152\\
28&156174269991049949415159064977550632619507763234249371405639558723862254437039266619661406863507459642213137821040024206588637964583113990836039863160192795531975\\
29&571597828167271067218399440081168442270698976507373070008012133133883453380504689240290656664762511770322706764692995231551093351998545664907795930839984031009342900\\
30&2161782986128627688674487022941484764264085795126614346894176622044781664752185330447863570782443307684690567863633426138786649244110884759590518145325249619492465122105\\
31&8443924343818422252327124346418501812521192060562901976321293041153812808547032176101080612448534038551518600575644441509018546291962326205363457209599533731669742493178300\\
32&34045902954275879238838781952237869914717029009153862598012645143793185682505682890913560628423517074335763040725808762980164777113191144138643425430538373554000767108953305897\\
33&141630956289787657835052982489896541299849963371241659857557366529535271858658792601814063749034247988010290255019769641224805060741489387389084381188038777839795982465934792272440\\
34&607596802483189052167849087786066300568149226281340652553417065409516803903596233671582415152616596908237280175624505961964561658039728720728101492630085291210213168030777969092878439\\
35&2686793060580661988701224190707481849509062656184826859748840495842141843076096622177873518157505941083819536735992495979152482796134866829948335976445595982052923978402797935132824791824\\
36&12241029184005496020526763249958377060918685384695429700602994614508702053248640654947103066060834364796869187307205678148560287657901843237670440849713123496403565513969018435947053834071435\\
37&57434908931353787328312616220387154445347277267653963162664968460863519622139632843493316682859817303912961109726068921703793397336923276837496259490934139548546739905364517343056248484856794608\\
38&277410610138438792795750205395400539554144450533760707286933345074731390379325562262143402089472868385023469208019446075742662911172035815273406948016136386103325984623055238402696544425782893748975\\
39&1378730732388040800194878589297516042821361764900317566981038235215111619992536037510537561980998193006243775777525153097916717386533961582507906747144268789872913814261110759204457104666411848984613304\\
40&7048071503967664570596219365706600444815293082566090993084936332242221474325113860197177375489744475265426942536386682736331456368828849564543673882957945349222382112199310103746877910535001657217053715419\\
41&37044663824854044983053728990433722638318178556588045365208875223269745876486114986967724830450842663301757104291035183337305555744909675245292078267363961058818417854608390963074852900597359305980676553340840\\
42&200115273981861550998456244007375674580302211967558006858041357896286878157400718365640083756027046576217963708225486029846256677361359793772805304824346062965028938905640123148726486904364868016308146243572451143.
\end{array}
\]

Now assume \(N\ge43\), and set \(A:=D_N\). Proposition 30 with
\(j=31,32,33\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{32}\): the base case is
\[
D_{43}=22225533213979647187685190410983617546032726150608122
>18663392025969146670022260524972832947786731075670401=43^{32},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{32}
\]
for \(r\ge43\) once \(D_r\ge r^{32}\), because
\(r^{33}\ge(r+1)^{32}\) on this range. The coefficient of \(A^2\) is
\[
(N+32)\bigl((N+1)\cdots(N+31)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{32}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(95\), denominator
\[
263130836933693530167218012160000000,
\]
and minimum coefficient
\[
526261673867387060334436024319999999.
\]
Thus both derivatives are positive for \(A\ge N^{32}\). At \(A=N^{32}\),
after substituting \(N=s+6\), the two parity numerators for \(\Delta\) have
positive coefficients; each has degree \(127\), denominator
\[
2163682417044562974786443716998736989202141037210684620800000000000000,
\]
and minimum coefficient
\[
2163682417044562974786443716998736980979302383032761803074437120000000.
\]
Hence \(\Delta>0\). Since \(m_{N+31}^{(N)}\) and \(m_{N+32}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+31}^{(N)}\le\rho_{N+32}^{(N)}\). \(\square\)

## Proposition 62 (thirty-fourth transition ratio)

For every \(N\ge6\),
\[
\rho_{N+32}^{(N)}\le\rho_{N+33}^{(N)}.
\]

Proof.
Set
\[
\Delta:=m_{N+34}^{(N)}m_{N+32}^{(N)}-(m_{N+33}^{(N)})^2.
\]
For \(N=6,7,\ldots,44\), direct substitution in the hook-length formula of
Proposition 30 gives the following positive integer values:
\[
\begin{array}{c|r}
N&\Delta\\ \hline
6&121677291308954583482672208535289611205332649072340790713426562741984303799077413750\\
7&128502275504958124230416237432758656257839071237453826207435618364125469177746172445371136\\
8&12941581376776516387873641414642198076560033651610454696575126853796558543214384921286422842000\\
9&282288812418574841224667090051616962063444831365794168836209982346702349215147414939885603613626096\\
10&2293125156996670090668247266731209786267991064052910740243495241732627723410135372541808138603975079950\\
11&10059528646875886357030701458276159735994444843892782421465359864509113194406233221747804828532905250381552\\
12&30805444061035708409539209150116037860299924865900721518959479196123311877619009158422023098684575996195887804\\
13&78284450912312767451200877038068336501102744604423638115547611364464634047935713107038822941270834467276087499776\\
14&184344745696245757943207648768274332936587214457870881621482300524007708472108964265623918623839114205143893937514906\\
15&428923794965198883298132538994065809494553082293728481004971368375182604892053685222193722434879394729864840166012081920\\
16&1018556035644853372833377459909287132788792203876881652362915433239129602501242507594619301441308449569470377537397895632625\\
17&2501880736988850291657526740249115216045027406584381331498322869938025112679161682446570249076342309379956410859857893095658556\\
18&6383610156185760367682029032113421899490492872908026418751667581115481815058688632537853955538155623470396086201350988251711048725\\
19&16931404141610385259453807430960361464396532236661190723031033296790605503409109651276584741128638183389917550339889851315440649170304\\
20&46663987311172433890930415424515796578130791962482266005785934689223730513188820871112881140354413736338086752963587500369783105978078005\\
21&133552816393037332273094993555280807780450302618487435543732302762317344567684527603893277843698777596673619136790920341970489886671637604300\\
22&396652077413306823450089785473541658048345348604990440720691994029581700104050344252813553200025780005164864634686672764597946115083292342080901\\
23&1221688486712209117901234899842277981291857073033886977475540770491807782418681083126129926716395913234743715838286255555297389303954953425288544924\\
24&3899629684421944435923154682054780093502482861489261004620804341093181582581844913261274690188567613097944215616198138476203539217728304993015830816856\\
25&12892175749834857007079130884248725879501439288138146592330335850925448904469986167765531225043589880506548121553021146865100495201027933714166675534080316\\
26&44117025420687900959858871324637878911581050044468904956623855791245655025803656368618565117520210529977915986986503268920587029388217230685196263197404883858\\
27&156174269990891611365021773282349266531648303891344487959264635673573828151078994398057620453934069503733166302804810563385130323950877498685156124165442668592944\\
28&571597828167221130193122924895998579103204025068583448925180943793912437248192718939665218469768283909050676829599662666260407982515992330102742007204925406630333025\\
29&2161782986128612375846437791865017880327501801352018051737290490870662855443887187397031693025986504007471271178126415431120963581373577189051470887066278704722626154000\\
30&8443924343818417676338822490051067439483678365583451112921689155709178681089317163751403577845252025575404994761742105648444929108046315019983779007914209457397529386616390\\
31&34045902954275877903431843179042573847663979070231233655069621339839773135632610672570303213692656618849241773897861196908233645095768852701015123120647179763868040793511701020\\
32&141630956289787657453760708354941893586197150802036027545420746244482261481529997953381426082445695065309081143234429378789495676698191647329000360407548108271445777179446271103261\\
33&607596802483189052061150362064021066941131157981556659663173952118141576118657303482550503508489850994825156027217440956685207017251357449299990796632166334011877463407131823576866800\\
34&2686793060580661988671915773396767412306184204910882791863822048539160701842795184888184998683778116715634312320971354881737814312769396605462793263074817325672340397131684591538090490141\\
35&12241029184005496020518849650959793308452273494744334843405924394590540688107427120617994517517493723354704312037839968235146553167142702613268732315418264745449413581055256970798593061686460\\
36&57434908931353787328310513076202748769696091833546088436036625693786771951335989673532301542119030374174347028851818285380811083290662220395262190004542540806932233122734381690084158820647584477\\
37&277410610138438792795749654592771608882446605463244696968551578568523594800986726017883341941983357593674676472441739707574978337389490192387144881461988307922806322861823514067672574898153644964416\\
38&1378730732388040800194878446988561644284317104769484323320772256273337289779746462863851989354796820187859847352701657643823546902386367912750396721726454037238660690627939907370706438106531090436174525\\
39&7048071503967664570596219329397292254220834898811163957184455817313201822246654094166339001093683534686337529113885634876293371308659645489307134636273124348267108478671984917138596294450973846266277853116\\
40&37044663824854044983053728981276576711959931816084929193889252945045112349171512022255753012817260450102646600080712028347133758086229454872311624355910756795847516030087690307822084554769041471389400782563165\\
41&200115273981861550998456244005090941501143129543818530685184570035636946050067845124913508740667675757002278203772772984713136584352337714639594557329324768820445505394878813758868771618139743108494484646671198192\\
42&1110639770599331608041432154240627011222488315463771056959375695152041471702360003023472364045583150914730892557353212259703435606321029961845174225149247581582239434217422737914672709154000105192249278419576025261541\\
43&6330646692416190165836163279174714842369206422659523673719778615017339322696123630071074980528264092581970694530320735223290382530070176697603692377386991433972950986181439694898242670048473305813316100249977272370591100\\
44&37046944444019544850473227509731221016415441371109851589934071252774534337219106346248304849706056478727176286040954082207891480845080997276917929629863562507980343307373851266704612512305753856370420727187059334171542737821.
\end{array}
\]

Now assume \(N\ge45\), and set \(A:=D_N\). Proposition 30 with
\(j=32,33,34\), together with
\[
D_{r+1}=(r+1)D_r+(-1)^{r+1},
\]
expresses \(\Delta\) as a quadratic in \(A\), separately in the two parities
of \(N\). Also \(A\ge N^{33}\): the base case is
\[
D_{45}=44006555763679701431616677013747562741144797778204081604
>3597600662921626628600135655553196556866168975830078125=45^{33},
\]
and
\[
D_{r+1}=r(D_r+D_{r-1})\ge rD_r\ge(r+1)^{33}
\]
for \(r\ge45\) once \(D_r\ge r^{33}\), because
\(r^{34}\ge(r+1)^{33}\) on this range. The coefficient of \(A^2\) is
\[
(N+33)\bigl((N+1)\cdots(N+32)\bigr)^2,
\]
so the derivative with respect to \(A\) is increasing in \(A\).

Verification. Substituting \(A=N^{33}\) into the two parity derivatives,
then substituting \(N=s+6\), gives numerator polynomials with positive
coefficients; each has degree \(98\), denominator
\[
8683317618811886495518194401280000000,
\]
and minimum coefficient
\[
17366635237623772991036388802559999999.
\]
Thus both derivatives are positive for \(A\ge N^{33}\). At \(A=N^{33}\),
after substituting \(N=s+6\), the two parity numerators for \(\Delta\) have
positive coefficients; each has degree \(131\), denominator
\[
2284848632399058501374484565150666260597460935294482959564800000000000000,
\]
and minimum coefficient
\[
2284848632399058501374484565150666260334330098360789429397581987840000000.
\]
Hence \(\Delta>0\). Since \(m_{N+32}^{(N)}\) and \(m_{N+33}^{(N)}\) are
positive by Proposition 42, this is equivalent to
\(\rho_{N+32}^{(N)}\le\rho_{N+33}^{(N)}\). \(\square\)

## Proposition 63 (fixed-band endpoint domination)

Let \(\mu\) be a probability measure on \([-1,A]\), where
\[
A>c>u>b>1.
\]
For an odd integer \(k\), put
\[
\Delta_k:=
\frac12\iint x^ky^k(x-y)^2\,d\mu(x)d\mu(y).
\]
Let
\[
p_B:=\mu([b,u]),\qquad p_C:=\mu([c,A]).
\]
If
\[
p_Bp_Cb^k(c-u)^2>2(c+1)^2,
\]
then \(\Delta_k>0\).

Proof.
The contribution from \([b,u]\times[c,A]\) and the symmetric rectangle is
\[
P:=\int_{[c,A]}\int_{[b,u]}x^ky^k(x-y)^2\,d\mu(y)d\mu(x),
\]
and therefore
\[
P\ge p_Bp_Cb^kc^k(c-u)^2.
\]
The negative part with the positive variable in \([c,A]\) has absolute value
at most
\[
H:=\int_{[c,A]}x^k(x+1)^2\,d\mu(x).
\]
For \(x\ge c\), the function \((x+1)/(x-u)\) is decreasing, so
\[
\frac{(x-u)^2}{(c-u)^2}\ge\frac{(x+1)^2}{(c+1)^2}.
\]
Thus
\[
p_Bb^k(x-u)^2\ge
p_Bb^k(c-u)^2\frac{(x-u)^2}{(c-u)^2}
>2(x+1)^2.
\]
Here the strict inequality follows from the displayed hypothesis, since
\(p_C\le1\). Hence \(P>2H\), and after paying for this high-positive
negative part, at least \(P/2\) remains.

The remaining negative part has the positive variable in \((0,c)\). Its
absolute value is bounded by
\[
L\le(c+1)^2c^k.
\]
The displayed hypothesis also gives
\[
\frac P2\ge
\frac12p_Bp_Cb^kc^k(c-u)^2
>(c+1)^2c^k\ge L.
\]
The same-sign regions contribute nonnegatively, so \(\Delta_k>0\).
\(\square\)

## Proposition 64 (endpoint overlap for \(N\ge30\))

For \(N\ge30\) and every odd \(k\ge89\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}.
\]

Proof.
Let \(U\) be Haar-distributed in \(SU(N)\),
\[
Y:=|\operatorname{tr}U|^2,\qquad X:=Y-1.
\]
The random variable \(Y\) is invariant under scalar phases, so its
\(SU(N)\) and \(U(N)\) moments agree. By Proposition 21,
\[
E[Y^r]=r!\qquad(0\le r\le30),
\]
because \(N\ge30\).

First define
\[
Q(y):=2418-15425y+15867y^2-3333y^3+246y^4-6y^5
\]
and
\[
S:=2418+15425\cdot5+15867\cdot5^2+3333\cdot5^3
+246\cdot5^4+6\cdot5^5=1065343.
\]
For \(y\in[5/2,5]\),
\[
0\le (y-5/2)(5-y)\le25/16,\qquad |Q(y)|\le S.
\]
Outside \([5/2,5]\), the factor \((y-5/2)(5-y)\) is nonpositive. Therefore
\[
{\bf 1}_{[5/2,5]}(y)\ge
\frac{(y-5/2)(5-y)Q(y)^2}{(25/16)S^2}
\qquad(y\ge0).
\]
Expanding the numerator and using \(E[Y^r]=r!\) for \(r\le12\) gives
\[
E[(Y-5/2)(5-Y)Q(Y)^2]=15280964.
\]
Hence
\[
\mu_N([3/2,4])
\ge
\delta_B:=
\frac{244495424}{25\cdot1065343^2}.
\]

For the tail, apply Paley-Zygmund to \(Z=Y^{15}\). Since
\[
E[Z]=15!,\qquad E[Z^2]=30!,
\]
and \(6^{15}<15!\), one has
\[
\mu_N([5,A])
=P(Y\ge6)
\ge
\delta_C:=
\left(1-\frac{6^{15}}{15!}\right)^2\frac{(15!)^2}{30!}
=\frac{471756427}{178410035553750000}.
\]

Direct rational arithmetic gives
\[
\delta_B\delta_C\left(\frac32\right)^{89}
=
\frac{
582582919223096898907695464995062465288803084365782528009
}{
5439830057284128878259185134333467458181005312000000000
}
>72.
\]
Thus for every \(k\ge89\),
\[
\delta_B\delta_C\left(\frac32\right)^k(5-4)^2>2(5+1)^2.
\]
Proposition 63 applies with
\[
b=3/2,\qquad u=4,\qquad c=5,\qquad A=N^2-1.
\]
Consequently
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. Proposition 42 gives positive denominators for
\(\rho_k^{(N)}\) and \(\rho_{k+1}^{(N)}\), since \(k\ge89\). Therefore the
displayed determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 65 (lower-rank endpoint overlap for \(N\ge24\))

For \(N\ge24\) and every odd \(k\ge71\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}.
\]

Proof.
Let
\[
Y:=|\operatorname{tr}U|^2,\qquad X:=Y-1,
\]
with \(U\) Haar-distributed in \(SU(N)\). As in Proposition 64,
\[
E[Y^r]=r!\qquad(0\le r\le24),
\]
because \(N\ge24\).

Set \(z:=y-7/2\), and define
\[
\begin{aligned}
R(z):={}&925503648+348262482z-142832255z^2-40966535z^3\\
&+8723312z^4+880405z^5-242416z^6+14686z^7-274z^8.
\end{aligned}
\]
Let
\[
T:=925503648+348262482+142832255+40966535+8723312
+880405+242416+14686+274=1467426013.
\]
For \(y\in[5/2,9/2]\), one has \(|z|\le1\),
\[
0\le(y-5/2)(9/2-y)=1-z^2\le1,\qquad |R(z)|\le T.
\]
Outside \([5/2,9/2]\), the factor \((y-5/2)(9/2-y)\) is nonpositive.
Therefore
\[
{\bf 1}_{[5/2,9/2]}(y)
\ge
\frac{(y-5/2)(9/2-y)R(y-7/2)^2}{T^2}
\qquad(y\ge0).
\]
Expanding the numerator and using \(E[Y^r]=r!\) for \(r\le18\) gives
\[
E[(Y-5/2)(9/2-Y)R(Y-7/2)^2]
=\frac{122923229137548665407}{65536}.
\]
Hence
\[
\mu_N([3/2,7/2])
\ge
\eta_B:=
\frac{122923229137548665407}{65536\cdot1467426013^2}.
\]

For the tail, apply Paley-Zygmund to \(Z=Y^{12}\). Since
\[
E[Z]=12!,\qquad E[Z^2]=24!,
\]
and \(5^{12}<12!\), one has
\[
\mu_N([4,A])
=P(Y\ge5)
\ge
\eta_C:=
\left(1-\frac{5^{12}}{12!}\right)^2\frac{(12!)^2}{24!}
=\frac{88255484124721}{992717442773183102976}.
\]
Direct rational arithmetic gives
\[
\eta_B\eta_C\left(\frac32\right)^{71}
=
\frac{
1379660402919710663977809665618779978693971095956987607403200541
}{
5601897465565742964999518305129780652282929769510225699667968
}
>200.
\]
Thus for every \(k\ge71\),
\[
\eta_B\eta_C\left(\frac32\right)^k(4-7/2)^2>2(4+1)^2.
\]
Proposition 63 applies with
\[
b=3/2,\qquad u=7/2,\qquad c=4,\qquad A=N^2-1.
\]
Consequently
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. Proposition 42 gives positive denominators for
\(\rho_k^{(N)}\) and \(\rho_{k+1}^{(N)}\), since \(k\ge71\). Therefore the
displayed determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 66 (finite overlap into the \(N\ge24\) endpoint range)

The transition inequality
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\]
holds for the following finite list of odd transition ranges:
\[
\begin{array}{c|c}
N&k\\ \hline
24&57,59,61,63,65,67,69\\
25&59,61,63,65,67,69\\
26&61,63,65,67,69\\
27&65,67,69.
\end{array}
\]

Proof.
Fix one of the listed pairs \((N,k)\), and put \(K:=k+2\). By
Proposition 21 and RSK, for \(a\ge N+1\),
\[
a!-E_a^{(N)}
=\#\{\pi\in S_a:\operatorname{lds}(\pi)>N\}.
\]
If \(\operatorname{lds}(\pi)>N\), then \(\pi\) has a decreasing subsequence
of length \(N+1\), so the union bound over positions gives
\[
0\le a!-E_a^{(N)}\le a!\frac{\binom a{N+1}}{(N+1)!}.
\]
Therefore, for \(r\le K\),
\[
\begin{aligned}
\left|m_r^{(N)}-D_r\right|
&\le\sum_{a=N+1}^r\binom ra(a!-E_a^{(N)})\\
&\le r!\frac{\binom K{N+1}}{(N+1)!}
\sum_{b=0}^{r-N-1}\frac1{b!}\\
&\le
3r!\frac{\binom K{N+1}}{(N+1)!}.
\end{aligned}
\]
As in Proposition 41, \(D_r\ge r!/3\), hence
\[
\left|m_r^{(N)}-D_r\right|\le\tau D_r,\qquad
\tau:=9\frac{\binom K{N+1}}{(N+1)!}.
\]
The same determinant expansion used in Proposition 41 gives
\[
\begin{aligned}
&m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2\\
&\quad\ge
D_kD_{k+1}
\left[
\frac7{18}
-\tau(4k+10)-\tau^2(2k+5)
\right].
\end{aligned}
\]
The bracket is a finite rational check. Over the listed odd \(k\)'s, its
minimum in each row is:
\[
\begin{array}{c|c|c}
N&\text{minimum attained at }k&\text{minimum bracket}\\ \hline
24&69&>3/10\\
25&69&>3/10\\
26&69&>3/10\\
27&69&>3/10.
\end{array}
\]
Hence the determinant is positive in every listed case. Verification.
Proposition 42 gives positive denominators for \(\rho_k^{(N)}\) and
\(\rho_{k+1}^{(N)}\), since all listed \(k\)'s are at least \(57\). Therefore
the determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 67 (endpoint overlap for \(6\le N\le23\))

For \(6\le N\le23\) and every odd \(k\ge53\),
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\]

Proof.
Let
\[
Y:=|\operatorname{tr}U|^2,\qquad X:=Y-1,
\]
with \(U\) Haar-distributed in \(SU(N)\). Put \(z:=y-43/10\), and define
\[
Q(z):=
959324176457+276017837985z-58083597219z^2-11664458272z^3
+1230335979z^4.
\]
On \(y\in[29/10,57/10]\), one has \(|z|\le7/5\) and
\[
0\le(y-29/10)(57/10-y)\le49/25.
\]
Also
\[
|Q(z)|\le
S:=959324176457+276017837985\frac75
+58083597219\left(\frac75\right)^2
+11664458272\left(\frac75\right)^3
+1230335979\left(\frac75\right)^4
=\frac{935204207737834}{625}.
\]
Outside \([29/10,57/10]\), the factor
\((y-29/10)(57/10-y)\) is nonpositive. Therefore
\[
{\bf 1}_{[29/10,57/10]}(y)\ge
\frac{(y-29/10)(57/10-y)Q(y-43/10)^2}{(49/25)S^2}
\qquad(y\ge0).
\]
Using Proposition 21 and the hook formula
\[
E[Y^r]=E_r^{(N)}
=\sum_{\lambda\vdash r,\ \ell(\lambda)\le N}(f^\lambda)^2,\qquad
f^\lambda
=r!\frac{\prod_{i<j}(\lambda_i-\lambda_j+j-i)}
{\prod_i(\lambda_i+\ell(\lambda)-i)!},
\]
the integrand has degree \(10\). For \(N\ge10\), Proposition 21 gives the
stable values \(E[Y^r]=r!\) for every \(r\le10\), and expansion gives
\[
E[(Y-29/10)(57/10-Y)Q(Y-43/10)^2]
=
\frac{44525786465344542780431067563007}{10000000000}.
\]
For \(6\le N\le9\), the same hook-formula computation gives a value at
least this large; the minimum excess occurs at \(N=9\) and equals
\[
1513726621221888441>0.
\]
Hence
\[
\mu_N([19/10,47/10])
\ge
\delta_B:=
\frac{44525786465344542780431067563007}
{43884276324717505323728973379833856}.
\]

Set
\[
\delta_C:=\frac{471756427}{178410035553750000}.
\]
For the tail, apply Paley-Zygmund to \(Z=Y^{15}\). For \(N\ge15\),
\(E[Y^{15}]=15!\), while \(E[Y^{30}]\le30!\), so
\[
P(Y\ge6)\ge
\left(1-\frac{6^{15}}{15!}\right)^2
\frac{(15!)^2}{30!}
=\delta_C.
\]
For \(6\le N\le14\), the same finite hook computation gives
\[
P(Y\ge6)\ge
\left(1-\frac{6^{15}}{E[Y^{15}]}\right)^2
\frac{E[Y^{15}]^2}{E[Y^{30}]}
\]
and the minimum excess over \(\delta_C\) occurs at \(N=14\), where it is
\[
\frac{
2602599137945210797524667371658649
}{
56336839730247629366109305276233265935886250000
}
>0.
\]
Direct rational arithmetic gives
\[
\delta_B\delta_C\left(\frac{19}{10}\right)^{52}
\left(\frac3{10}\right)^2
=
\frac{
3457507041948414895333087349928700632838006234393012741595904511533668451042683761282837131757375075222791
}{
45785937423061049735271877451609038317122560000000000000000000000000000000000000000000000000000000000
}
>72.
\]
Thus for every \(k\ge53\),
\[
\delta_B\delta_C\left(\frac{19}{10}\right)^k
\left(5-\frac{47}{10}\right)^2>2(5+1)^2.
\]
Proposition 63 applies with
\[
b=19/10,\qquad u=47/10,\qquad c=5,\qquad A=N^2-1.
\]
Consequently
\[
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2>0.
\]
Verification. Proposition 42 gives positive denominators for
\(\rho_k^{(N)}\) and \(\rho_{k+1}^{(N)}\), since \(k\ge53\). Therefore the
displayed determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

## Proposition 68 (final finite overlap for \(6\le N\le18\))

The transition inequality
\[
\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\]
holds for the following finite list of odd transition ranges:
\[
\begin{array}{c|c}
N&k\\ \hline
6&39,41,43,45,47,49,51\\
7&41,43,45,47,49,51\\
8&41,43,45,47,49,51\\
9&43,45,47,49,51\\
10&43,45,47,49,51\\
11&45,47,49,51\\
12&45,47,49,51\\
13&47,49,51\\
14&47,49,51\\
15&49,51\\
16&49,51\\
17&51\\
18&51.
\end{array}
\]

Proof.
For each listed pair \((N,k)\), compute
\[
\Delta_{N,k}:=
m_{k+2}^{(N)}m_k^{(N)}-(m_{k+1}^{(N)})^2
\]
from Proposition 21 using the same hook-formula evaluation of
\(E_r^{(N)}\) as in Proposition 67. The exact finite check gives the
following row minima:
\[
\begin{array}{c|c|r}
N&\text{minimum attained at }k&\min_k\Delta_{N,k}\\ \hline
6&39&64535719854392029544773326762011774756925778563684823690382903928732901693552607059275\\
7&41&81983250005769227285498956748387927182179586638327455931621250522805530412550262461119689354176\\
8&41&13776730806844857067459723475572806115606561762537749689970379383403774177954665948696550730822960\\
9&43&516741713808588211599059495979604717390970752707423594032209037987274754148154730954160916056953963450368\\
10&43&3618232844122202584832213951704951536993240226634987294481850293593859749135402230512808610277609410270817\\
11&45&33399774449175848055011934469656070323910964218709116241112806497510841469515062121126196910235665905718776943228\\
12&45&60567314448754365018106205779950325082006402099956453971960484870790769677108940276606446856248401017048182106706\\
13&47&363239305623168787675421494307053710143238482046635190860778910231693491256558565289022469174629997636148032811272314736\\
14&47&412017453922855017028903546300156146196642966415729406173097639746151298064337858733695094575907116122702045750587973615\\
15&49&2454353316530929119018984347889491575535943771589993590419065745590973632964833686716148957357041276645873773505319638231954239\\
16&49&2493119549736314583261850567832642447348332533596249440198341404887202602996116002199771825733723134557950156991608933086839975\\
17&51&16910085949149500902700908171041557824229854633053201817574276758325371717389977268687410623596430698509784622587525716043509684517756\\
18&51&16928441762086713834411725710317433241255361224654439216646167863265666371228166615581072054953650607845481890940623126734747743109775.
\end{array}
\]
Every listed minimum is positive, so \(\Delta_{N,k}>0\) for every listed
pair. Verification. Proposition 42 gives positive denominators for
\(\rho_k^{(N)}\) and \(\rho_{k+1}^{(N)}\), since every listed \(k\) is at
least \(39\). Therefore the determinant inequality is equivalent to
\(\rho_k^{(N)}\le\rho_{k+1}^{(N)}\). \(\square\)

The companion script `verify_sun_repair_certificates.py` recomputes the
endpoint constants and the row minima in Propositions 64, 65, 66, 67, and
68 from the hook formula of Proposition 21, using exact rational arithmetic.

## Corollary 69 (SUN-RatioMono discharge for \(N\ge6\))

For every \(N\ge6\), Assumption SUN-RatioMono holds:
\[
\rho_3^{(N)}=\frac92,\qquad
\frac92\le\rho_k^{(N)}\le\rho_{k+1}^{(N)}
\quad(k\ge3).
\]

Proof.
Proposition 23 gives the displayed initial equality and the inequalities
for \(3\le k\le N-2\). Proposition 42 gives the inequality for every even
\(k\ge2\). It remains to cover odd \(k\ge N-1\).

For \(N-1\le k\le N+33\), the transition inequalities follow from the
uniform transition-strip lemma in `paper.tex`. Hence assume
\(k\ge N+34\). If \(6\le N\le18\), then Proposition 68 covers the listed
odd \(k\)'s up to \(51\), and Proposition 67 covers every odd \(k\ge53\).
If \(19\le N\le23\), then \(k\ge N+34\ge53\), so Proposition 67 applies.
If \(24\le N\le27\), Proposition 66 covers the remaining odd \(k\)'s below
\(71\), and Proposition 65 covers every odd \(k\ge71\). If \(N\ge28\) and
\(k\le69\), then \(k\ge N+34>N+16\), and Proposition 41 applies; the upper
bound in Proposition 41 is at least \(69\) at \(N=28\) and increases for
\(N\ge28\). Finally, if \(N\ge24\) and \(k\ge71\), Proposition 65 applies.
These cases cover all odd \(k\ge N-1\). \(\square\)

## Corollary 70 (finite-rank SU(N) tail for \(N\ge6\))

For \(N\ge6\) and \(N<n+2\),
\[
Q_3^{SU(N),{\rm adj}}(n)\ge0.
\]

Proof.
If \(n\) is even, the integrand is pointwise nonnegative. If \(n\) is odd,
then \(N<n+2\) and \(N\ge6\) imply \(n\ge5\). Corollary 69 discharges
Assumption SUN-RatioMono, so Proposition 22 gives
\[
Q_3^{SU(N),{\rm adj}}(n)>0.
\]
\(\square\)

## Theorem 71 (SU(N) adjoint Q_3 closure)

For every \(N\ge2\) and every \(n\ge0\),
\[
Q_3^{SU(N),{\rm adj}}(n)\ge0.
\]

Proof.
Theorem 6 proves the claim for \(N=2\). Corollary 12 proves the claim for
\(N=3\), Corollary 16 for \(N=4\), and Corollary 20 for \(N=5\). Now let
\(N\ge6\). If \(N\ge n+2\), Theorem 5 applies. If \(N<n+2\), Corollary 70
applies. \(\square\)
