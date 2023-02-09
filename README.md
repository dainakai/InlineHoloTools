# InlineHoloTools
このパッケージはインラインホログラフィ(Gabor, 1948)に基づくホログラムの記録・再生のための関数群を提供します。

This package provides a toolset for light propagation calculations in inline holography (Gabor, 1948). 

## ホログラムの記録 Recording holograms
$(x,y,z)$ の右手系デカルト座標で、 $z$ 軸正方向にコヒーレントな平行波が伝搬しているとき、この平行波の振幅を $B$、波長を $\lambda$ で表現する。 $xy$ 平面 $z=0$ で物体分布を表現する透過関数を定義し、これを $A_{0}(x,y)$ とおく。すなわち、物体が存在している点では $A$ は値として $0$ を取り、それ以外の点では $1$ を取る。透過関数を定義した $z=0$
面でのコヒーレントな平行波の位相を $0$ とすると、この面での平行波 $\psi_{z=0-\mathrm{dz}}(x,y)$　は $\psi_{z=0-\mathrm{dz}}(x,y) = B$ と表現できるから、透過関数で表現された物体を透過した瞬間の波の複素振幅 $\psi_{0}(x,y)$ は 以下で表される。
$$\psi_{0}(x,y)=BA_{0}(x,y)$$
この光波面がさらに $z_0$ 伝搬して $z=z_0$ における複素振幅は、Rayleigh-Sommerfeld回折積分に基づく角スペクトル法によって、以下で表される。
$$\psi_{z_0}(x,y)=\mathcal{F}^{-1}\lbrace \mathcal{F}\{\psi_0(x,y)\} \cdot G_{z_0} \rbrace $$
ここで、 $\mathcal{F}$ は $x,y$ に関する二次元フーリエ変換を表し、 $G_{z_0}$ は以下に定義する光伝搬を表現する伝達関数である。
$$G_{z_0}(\alpha,\beta)=\exp \left[\frac{2\pi i}{\lambda} \sqrt{1-\alpha^2 - \beta^2}\right]$$
$\alpha, \beta$ は $x,y$ それぞれに対応する周波数空間の変数である。このようにして透過関数の任意の $z$ 軸方向の光伝搬を表現できる。記録されるホログラムは、光伝搬した複素振幅の強度によって取得できる。
$$I(x,y) = ||\psi_{z_0}(x,y)||$$