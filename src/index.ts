"use strict";
export const QRScanner = ()=> {

    let MODULE = {
        initiate: undefined
    };

    function QRCodeDecoder(scope) {
        let _ay: any;
        if (!scope.qrcode) {
            let _aa = {
                _ab: undefined,
                _af: undefined,
                _ah: undefined
            };
            _aa._ab = function (f, e) {
                let d = qrcode.width;
                let b = qrcode.height;
                let c = true;
                for (let g = 0; g < e.length && c; g += 2) {
                    let a = Math.floor(e[g]);
                    let h = Math.floor(e[g + 1]);
                    if (a < -1 || a > d || h < -1 || h > b) {
                        throw"Error._ab "
                    }
                    c = false;
                    if (a == -1) {
                        e[g] = 0;
                        c = true
                    } else {
                        if (a == d) {
                            e[g] = d - 1;
                            c = true
                        }
                    }
                    if (h == -1) {
                        e[g + 1] = 0;
                        c = true
                    } else {
                        if (h == b) {
                            e[g + 1] = b - 1;
                            c = true
                        }
                    }
                }
                c = true;
                for (let g = e.length - 2; g >= 0 && c; g -= 2) {
                    let a = Math.floor(e[g]);
                    let h = Math.floor(e[g + 1]);
                    if (a < -1 || a > d || h < -1 || h > b) {
                        throw"Error._ab "
                    }
                    c = false;
                    if (a == -1) {
                        e[g] = 0;
                        c = true
                    } else {
                        if (a == d) {
                            e[g] = d - 1;
                            c = true
                        }
                    }
                    if (h == -1) {
                        e[g + 1] = 0;
                        c = true
                    } else {
                        if (h == b) {
                            e[g + 1] = b - 1;
                            c = true
                        }
                    }
                }
            };
            _aa._af = function (b, d, a) {
                // @ts-ignore
                let k = new _ac(d);
                let j = new Array(d << 1);
                for (let f = 0; f < d; f++) {
                    let g = j.length;
                    let i = f + 0.5;
                    for (let h = 0; h < g; h += 2) {
                        j[h] = (h >> 1) + 0.5;
                        j[h + 1] = i
                    }
                    a._ad(j);
                    _aa._ab(b, j);
                    try {
                        for (let h = 0; h < g; h += 2) {
                            let e = b[Math.floor(j[h]) + qrcode.width * Math.floor(j[h + 1])];
                            if (e) {
                                k._dq(h >> 1, f)
                            }
                        }
                    } catch (c) {
                        throw"Error._ab"
                    }
                }
                return k
            };
            _aa._ah = function (h, o, l, k, q, p, b, a, f, e, n, m, s, r, d, c, j, i) {
                let g = _ae._ag(l, k, q, p, b, a, f, e, n, m, s, r, d, c, j, i);
                return _aa._af(h, o, g)
            };

            // @ts-ignore
            const _a1:any=(b, a)=>{
                this.count = b;
                this._fc = a;
                this.__defineGetter__("Count", function () {
                    return this.count
                });
                this.__defineGetter__("_dm", function () {
                    return this._fc
                })
            };

            // @ts-ignore
            const _a2: any = (a, c, b) => {
                this._bm = a;
                if (b) {
                    this._do = [c, b]
                } else {
                    this._do = new Array(c)
                }
                this.__defineGetter__("_bo", function () {
                    return this._bm
                });
                this.__defineGetter__("_dn", function () {
                    return this._bm * this._fo
                });
                this.__defineGetter__("_fo", function () {
                    let e = 0;
                    for (let d = 0; d < this._do.length; d++) {
                        e += this._do[d].length
                    }
                    return e
                });
                this._fb = function () {
                    return this._do
                }
            };


            const _a3: any = (k, l, h, g, f, e) => {
                this._bs = k;
                this._ar = l;
                this._do = [h, g, f, e];
                let j = 0;
                let b = h._bo;
                let a = h._fb();
                for (let d = 0; d < a.length; d++) {
                    let c = a[d];
                    j += c.Count * (c._dm + b)
                }
                this._br = j;
                this.__defineGetter__("_fd", function () {
                    return this._bs
                });
                this.__defineGetter__("_as", function () {
                    return this._ar
                });
                this.__defineGetter__("_dp", function () {
                    return this._br
                });
                this.__defineGetter__("_cr", function () {
                    return 17 + 4 * this._bs
                });
                this._aq = function () {
                    let q = this._cr;
                    // @ts-ignore
                    let o = new _ac(q);
                    o._bq(0, 0, 9, 9);
                    o._bq(q - 8, 0, 8, 9);
                    o._bq(0, q - 8, 9, 8);
                    let n = this._ar.length;
                    for (let m = 0; m < n; m++) {
                        let p = this._ar[m] - 2;
                        for (let r = 0; r < n; r++) {
                            if ((m == 0 && (r == 0 || r == n - 1)) || (m == n - 1 && r == 0)) {
                                continue
                            }
                            o._bq(this._ar[r] - 2, p, 5, 5)
                        }
                    }
                    o._bq(6, 9, 1, q - 17);
                    o._bq(9, 6, q - 17, 1);
                    if (this._bs > 6) {
                        o._bq(q - 11, 0, 3, 6);
                        o._bq(0, q - 11, 6, 3)
                    }
                    return o
                };
                this._bu = function (i) {
                    return this._do[i.ordinal()]
                }
            };
            _a3._bv = [31892, 34236, 39577, 42195, 48118, 51042, 55367, 58893, 63784, 68472, 70749, 76311, 79154, 84390, 87683, 92361, 96236, 102084, 102881, 110507, 110734, 117786, 119615, 126325, 127568, 133589, 136944, 141498, 145311, 150283, 152622, 158308, 161089, 167017];
            _a3.VERSIONS = _ay();
            _a3._av = function (a) {
                if (a < 1 || a > 40) {
                    throw"bad arguments"
                }
                return _a3.VERSIONS[a - 1]
            };
            _a3._at = function (b) {
                if (b % 4 != 1) {
                    throw"Error _at"
                }
                try {
                    return _a3._av((b - 17) >> 2)
                } catch (a) {
                    throw"Error _av"
                }
            };
            _a3._aw = function (d) {
                let b = 4294967295;
                let f = 0;
                for (let c = 0; c < _a3._bv.length; c++) {
                    let a = _a3._bv[c];
                    if (a == d) {
                        return this._av(c + 7)
                    }
                    let e = _ax._gj(d, a);
                    if (e < b) {
                        f = c + 7;
                        b = e
                    }
                }
                if (b <= 3) {
                    return this._av(f)
                }
                return null
            };


            const _ae: any = (i, f, c, h, e, b, g, d, a) => {
                this.a11 = i;
                this.a12 = h;
                this.a13 = g;
                this.a21 = f;
                this.a22 = e;
                this.a23 = d;
                this.a31 = c;
                this.a32 = b;
                this.a33 = a;
                this._ad = function (v) {
                    let s = v.length;
                    let z = this.a11;
                    let w = this.a12;
                    let u = this.a13;
                    let q = this.a21;
                    let p = this.a22;
                    let o = this.a23;
                    let m = this.a31;
                    let k = this.a32;
                    let j = this.a33;
                    for (let n = 0; n < s; n += 2) {
                        let t = v[n];
                        let r = v[n + 1];
                        let l = u * t + o * r + j;
                        v[n] = (z * t + q * r + m) / l;
                        v[n + 1] = (w * t + p * r + k) / l
                    }
                };
                this._fp = function (m, k) {
                    let q = m.length;
                    for (let l = 0; l < q; l++) {
                        let j = m[l];
                        let p = k[l];
                        let o = this.a13 * j + this.a23 * p + this.a33;
                        m[l] = (this.a11 * j + this.a21 * p + this.a31) / o;
                        k[l] = (this.a12 * j + this.a22 * p + this.a32) / o
                    }
                };
                this._fr = function () {
                    return new _ae(this.a22 * this.a33 - this.a23 * this.a32, this.a23 * this.a31 - this.a21 * this.a33, this.a21 * this.a32 - this.a22 * this.a31, this.a13 * this.a32 - this.a12 * this.a33, this.a11 * this.a33 - this.a13 * this.a31, this.a12 * this.a31 - this.a11 * this.a32, this.a12 * this.a23 - this.a13 * this.a22, this.a13 * this.a21 - this.a11 * this.a23, this.a11 * this.a22 - this.a12 * this.a21)
                };
                this.times = function (j) {
                    return new _ae(this.a11 * j.a11 + this.a21 * j.a12 + this.a31 * j.a13, this.a11 * j.a21 + this.a21 * j.a22 + this.a31 * j.a23, this.a11 * j.a31 + this.a21 * j.a32 + this.a31 * j.a33, this.a12 * j.a11 + this.a22 * j.a12 + this.a32 * j.a13, this.a12 * j.a21 + this.a22 * j.a22 + this.a32 * j.a23, this.a12 * j.a31 + this.a22 * j.a32 + this.a32 * j.a33, this.a13 * j.a11 + this.a23 * j.a12 + this.a33 * j.a13, this.a13 * j.a21 + this.a23 * j.a22 + this.a33 * j.a23, this.a13 * j.a31 + this.a23 * j.a32 + this.a33 * j.a33)
                }
            };

            _ae._ag = function (p, e, o, d, n, c, m, b, h, q, l, f, a, j, i, r) {
                let g = this._be(p, e, o, d, n, c, m, b);
                let k = this._bf(h, q, l, f, a, j, i, r);
                return k.times(g)
            };
            _ae._bf = function (d, p, c, m, b, k, a, j) {
                let h = j - k;
                let f = p - m + k - j;
                if (h == 0 && f == 0) {
                    return new _ae(c - d, b - c, d, m - p, k - m, p, 0, 0, 1)
                } else {
                    let q = c - b;
                    let o = a - b;
                    let l = d - c + b - a;
                    let i = m - k;
                    let e = q * h - o * i;
                    let n = (l * h - o * f) / e;
                    let g = (q * f - l * i) / e;
                    return new _ae(c - d + n * c, a - d + g * a, d, m - p + n * m, j - p + g * j, p, n, g, 1)
                }
            };
            _ae._be = function (f, h, d, g, b, e, a, c) {
                return this._bf(f, h, d, g, b, e, a, c)._fr()
            };

            const _bg: any = (b, a) => {
                this.bits = b;
                this.points = a
            };

            const Detector: any = (a) => {
                this.image = a;
                this._am = null;
                this._bi = function (m, l, c, b) {
                    let d = Math.abs(b - l) > Math.abs(c - m);
                    if (d) {
                        let r = m;
                        m = l;
                        l = r;
                        r = c;
                        c = b;
                        b = r
                    }
                    let j = Math.abs(c - m);
                    let i = Math.abs(b - l);
                    let p = -j >> 1;
                    let u = l < b ? 1 : -1;
                    let f = m < c ? 1 : -1;
                    let e = 0;
                    for (let h = m, g = l; h != c; h += f) {
                        let t = d ? g : h;
                        let s = d ? h : g;
                        if (e == 1) {
                            if (this.image[t + s * qrcode.width]) {
                                e++
                            }
                        } else {
                            if (!this.image[t + s * qrcode.width]) {
                                e++
                            }
                        }
                        if (e == 3) {
                            let o = h - m;
                            let n = g - l;
                            return Math.sqrt((o * o + n * n))
                        }
                        p += i;
                        if (p > 0) {
                            if (g == b) {
                                break
                            }
                            g += u;
                            p -= j
                        }
                    }
                    let k = c - m;
                    let q = b - l;
                    return Math.sqrt((k * k + q * q))
                };
                this._bh = function (i, g, h, f) {
                    let b = this._bi(i, g, h, f);
                    let e = 1;
                    let d = i - (h - i);
                    if (d < 0) {
                        e = i / (i - d);
                        d = 0
                    } else {
                        if (d >= qrcode.width) {
                            e = (qrcode.width - 1 - i) / (d - i);
                            d = qrcode.width - 1
                        }
                    }
                    let c = Math.floor(g - (f - g) * e);
                    e = 1;
                    if (c < 0) {
                        e = g / (g - c);
                        c = 0
                    } else {
                        if (c >= qrcode.height) {
                            e = (qrcode.height - 1 - g) / (c - g);
                            c = qrcode.height - 1
                        }
                    }
                    d = Math.floor(i + (d - i) * e);
                    b += this._bi(i, g, d, c);
                    return b - 1
                };
                this._bj = function (c, d) {
                    let b = this._bh(Math.floor(c.X), Math.floor(c.Y), Math.floor(d.X), Math.floor(d.Y));
                    let e = this._bh(Math.floor(d.X), Math.floor(d.Y), Math.floor(c.X), Math.floor(c.Y));
                    if (isNaN(b)) {
                        return e / 7
                    }
                    if (isNaN(e)) {
                        return b / 7
                    }
                    return (b + e) / 14
                };
                this._bk = function (d, c, b) {
                    return (this._bj(d, c) + this._bj(d, b)) / 2
                };
                this.distance = function (d, b) {
                    let e = d.X - b.X;
                    let c = d.Y - b.Y;
                    return Math.sqrt((e * e + c * c))
                };
                this._bx = function (g, f, d, e) {
                    let b = Math.round(this.distance(g, f) / e);
                    let c = Math.round(this.distance(g, d) / e);
                    let h = ((b + c) >> 1) + 7;
                    switch (h & 3) {
                        case 0:
                            h++;
                            break;
                        case 2:
                            h--;
                            break;
                        case 3:
                            throw"Error"
                    }
                    return h
                };
                this._bl = function (g, f, d, j) {
                    let k = Math.floor(j * g);
                    let h = Math.max(0, f - k);
                    let i = Math.min(qrcode.width - 1, f + k);
                    if (i - h < g * 3) {
                        throw"Error"
                    }
                    let b = Math.max(0, d - k);
                    let c = Math.min(qrcode.height - 1, d + k);
                    let e = new _ak(this.image, h, b, i - h, c - b, g, this._am);
                    return e.find()
                };
                this.createTransform = function (l, h, k, b, g) {
                    let j = g - 3.5;
                    let i;
                    let f;
                    let e;
                    let c;
                    if (b != null) {
                        i = b.X;
                        f = b.Y;
                        e = c = j - 3
                    } else {
                        i = (h.X - l.X) + k.X;
                        f = (h.Y - l.Y) + k.Y;
                        e = c = j
                    }
                    let d = _ae._ag(3.5, 3.5, j, 3.5, e, c, 3.5, j, l.X, l.Y, h.X, h.Y, i, f, k.X, k.Y);
                    return d
                };
                this._bz = function (e, b, d) {
                    let c = _aa;
                    return c._af(e, d, b)
                };
                this._cd = function (q) {
                    let j = q._gq;
                    let h = q._gs;
                    let n = q._gp;
                    let d = this._bk(j, h, n);
                    if (d < 1) {
                        throw"Error"
                    }
                    let r = this._bx(j, h, n, d);
                    let b = _a3._at(r);
                    let k = b._cr - 7;
                    let l = null;
                    if (b._as.length > 0) {
                        let f = h.X - j.X + n.X;
                        let e = h.Y - j.Y + n.Y;
                        let c = 1 - 3 / k;
                        let t = Math.floor(j.X + c * (f - j.X));
                        let s = Math.floor(j.Y + c * (e - j.Y));
                        for (let p = 4; p <= 16; p <<= 1) {
                            l = this._bl(d, t, s, p);
                            break
                        }
                    }
                    let g = this.createTransform(j, h, n, l, r);
                    let m = this._bz(this.image, g, r);
                    let o;
                    if (l == null) {
                        o = [n, j, h]
                    } else {
                        o = [n, j, h, l]
                    }
                    return new _bg(m, o)
                };
                this.detect = function () {
                    let b = new _cc()._ce(this.image);
                    return this._cd(b)
                }
            };

            let _ca = 21522;
            let _cb = [[21522, 0], [20773, 1], [24188, 2], [23371, 3], [17913, 4], [16590, 5], [20375, 6], [19104, 7], [30660, 8], [29427, 9], [32170, 10], [30877, 11], [26159, 12], [25368, 13], [27713, 14], [26998, 15], [5769, 16], [5054, 17], [7399, 18], [6608, 19], [1890, 20], [597, 21], [3340, 22], [2107, 23], [13663, 24], [12392, 25], [16177, 26], [14854, 27], [9396, 28], [8579, 29], [11994, 30], [11245, 31]];
            let _ch = [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4];

            const _ax: any = (a) => {
                this._cf = _cg.forBits((a >> 3) & 3);
                this._fe = (a & 7);
                this.__defineGetter__("_cg", function () {
                    return this._cf
                });
                this.__defineGetter__("_dx", function () {
                    return this._fe
                });
                this.GetHashCode = function () {
                    return (this._cf.ordinal() << 3) | this._fe
                };
                this.Equals = function (c) {
                    let b = c;
                    return this._cf == b._cf && this._fe == b._fe
                }
            };

            _ax._gj = function (d, c) {
                d ^= c;
                return _ch[d & 15] + _ch[(_ew(d, 4) & 15)] + _ch[(_ew(d, 8) & 15)] + _ch[(_ew(d, 12) & 15)] + _ch[(_ew(d, 16) & 15)] + _ch[(_ew(d, 20) & 15)] + _ch[(_ew(d, 24) & 15)] + _ch[(_ew(d, 28) & 15)]
            };
            _ax._ci = function (a) {
                let b = _ax._cj(a);
                if (b != null) {
                    return b
                }
                return _ax._cj(a ^ _ca)
            };
            _ax._cj = function (d) {
                let b = 4294967295;
                let a = 0;
                for (let c = 0; c < _cb.length; c++) {
                    let g = _cb[c];
                    let f = g[0];
                    if (f == d) {
                        return new _ax(g[1])
                    }
                    let e = this._gj(d, f);
                    if (e < b) {
                        a = g[1];
                        b = e
                    }
                }
                if (b <= 3) {
                    return new _ax(a)
                }
                return null
            };

            const _cg: any = (a, c, b) => {
                this._ff = a;
                this.bits = c;
                this.name = b;
                this.__defineGetter__("Bits", function () {
                    return this.bits
                });
                this.__defineGetter__("Name", function () {
                    return this.name
                });
                this.ordinal = function () {
                    return this._ff
                }
            };

            _cg.forBits = function (a) {
                if (a < 0 || a >= FOR_BITS.length) {
                    throw"bad arguments"
                }
                return FOR_BITS[a]
            };
            let L = new _cg(0, 1, "L");
            let M = new _cg(1, 0, "M");
            let Q = new _cg(2, 3, "Q");
            let H = new _cg(3, 2, "H");
            let FOR_BITS = [M, L, H, Q];

            const _ac: any = (d, a) => {
                if (!a) {
                    a = d
                }
                if (d < 1 || a < 1) {
                    throw"Both dimensions must be greater than 0"
                }
                this.width = d;
                this.height = a;
                let c = d >> 5;
                if ((d & 31) != 0) {
                    c++
                }
                this.rowSize = c;
                this.bits = new Array(c * a);
                for (let b = 0; b < this.bits.length; b++) {
                    this.bits[b] = 0
                }
                this.__defineGetter__("Width", function () {
                    return this.width
                });
                this.__defineGetter__("Height", function () {
                    return this.height
                });
                this.__defineGetter__("Dimension", function () {
                    if (this.width != this.height) {
                        throw"Can't call getDimension() on a non-square matrix"
                    }
                    return this.width
                });
                this._ds = function (e, g) {
                    let f = g * this.rowSize + (e >> 5);
                    return ((_ew(this.bits[f], (e & 31))) & 1) != 0
                };
                this._dq = function (e, g) {
                    let f = g * this.rowSize + (e >> 5);
                    this.bits[f] |= 1 << (e & 31)
                };
                this.flip = function (e, g) {
                    let f = g * this.rowSize + (e >> 5);
                    this.bits[f] ^= 1 << (e & 31)
                };
                this.clear = function () {
                    let e = this.bits.length;
                    for (let f = 0; f < e; f++) {
                        this.bits[f] = 0
                    }
                };
                this._bq = function (g, j, f, m) {
                    if (j < 0 || g < 0) {
                        throw"Left and top must be nonnegative"
                    }
                    if (m < 1 || f < 1) {
                        throw"Height and width must be at least 1"
                    }
                    let l = g + f;
                    let e = j + m;
                    if (e > this.height || l > this.width) {
                        throw"The region must fit inside the matrix"
                    }
                    for (let i = j; i < e; i++) {
                        let h = i * this.rowSize;
                        for (let k = g; k < l; k++) {
                            this.bits[h + (k >> 5)] |= 1 << (k & 31)
                        }
                    }
                }
            };

            const _dl: any = (a, b) => {
                this._dv = a;
                this._dw = b;
                this.__defineGetter__("_du", function () {
                    return this._dv
                });
                this.__defineGetter__("Codewords", function () {
                    return this._dw
                })
            };

            _dl._gn = function (c, h, r) {
                if (c.length != h._dp) {
                    throw"bad arguments"
                }
                let k = h._bu(r);
                let e = 0;
                let d = k._fb();
                for (let q = 0; q < d.length; q++) {
                    e += d[q].Count
                }
                let l = new Array(e);
                let n = 0;
                for (let o = 0; o < d.length; o++) {
                    let f = d[o];
                    for (let q = 0; q < f.Count; q++) {
                        let m = f._dm;
                        let s = k._bo + m;
                        l[n++] = new _dl(m, new Array(s))
                    }
                }
                let t = l[0]._dw.length;
                let b = l.length - 1;
                while (b >= 0) {
                    let v = l[b]._dw.length;
                    if (v == t) {
                        break
                    }
                    b--
                }
                b++;
                let g = t - k._bo;
                let a = 0;
                for (let q = 0; q < g; q++) {
                    for (let o = 0; o < n; o++) {
                        l[o]._dw[q] = c[a++]
                    }
                }
                for (let o = b; o < n; o++) {
                    l[o]._dw[g] = c[a++]
                }
                let p = l[0]._dw.length;
                for (let q = g; q < p; q++) {
                    for (let o = 0; o < n; o++) {
                        let u = o < b ? q : q + 1;
                        l[o]._dw[u] = c[a++]
                    }
                }
                return l
            };

            const _cl: any = (a) => {
                let b = a.Dimension;
                if (b < 21 || (b & 3) != 1) {
                    throw"Error _cl"
                }
                this._au = a;
                this._cp = null;
                this._co = null;
                this._dk = function (d, c, e) {
                    return this._au._ds(d, c) ? (e << 1) | 1 : e << 1
                };
                this._cm = function () {
                    if (this._co != null) {
                        return this._co
                    }
                    let g = 0;
                    for (let e = 0; e < 6; e++) {
                        g = this._dk(e, 8, g)
                    }
                    g = this._dk(7, 8, g);
                    g = this._dk(8, 8, g);
                    g = this._dk(8, 7, g);
                    for (let c = 5; c >= 0; c--) {
                        g = this._dk(8, c, g)
                    }
                    this._co = _ax._ci(g);
                    if (this._co != null) {
                        return this._co
                    }
                    let f = this._au.Dimension;
                    g = 0;
                    let d = f - 8;
                    for (let e = f - 1; e >= d; e--) {
                        g = this._dk(e, 8, g)
                    }
                    for (let c = f - 7; c < f; c++) {
                        g = this._dk(8, c, g)
                    }
                    this._co = _ax._ci(g);
                    if (this._co != null) {
                        return this._co
                    }
                    throw"Error _cm"
                };
                this._cq = function () {
                    if (this._cp != null) {
                        return this._cp
                    }
                    let h = this._au.Dimension;
                    let f = (h - 17) >> 2;
                    if (f <= 6) {
                        return _a3._av(f)
                    }
                    let g = 0;
                    let e = h - 11;
                    for (let c = 5; c >= 0; c--) {
                        for (let d = h - 9; d >= e; d--) {
                            g = this._dk(d, c, g)
                        }
                    }
                    this._cp = _a3._aw(g);
                    if (this._cp != null && this._cp._cr == h) {
                        return this._cp
                    }
                    g = 0;
                    for (let d = 5; d >= 0; d--) {
                        for (let c = h - 9; c >= e; c--) {
                            g = this._dk(d, c, g)
                        }
                    }
                    this._cp = _a3._aw(g);
                    if (this._cp != null && this._cp._cr == h) {
                        return this._cp
                    }
                    throw"Error _cq"
                };
                this._gk = function () {
                    let q = this._cm();
                    let o = this._cq();
                    let c = _dx._gl(q._dx);
                    let f = this._au.Dimension;
                    c._dj(this._au, f);
                    let k = o._aq();
                    let n:any = true;
                    let r = new Array(o._dp);
                    let m = 0;
                    let p = 0;
                    let h = 0;
                    for (let e = f - 1; e > 0; e -= 2) {
                        if (e == 6) {
                            e--
                        }
                        for (let l = 0; l < f; l++) {
                            let g = n ? f - 1 - l : l;
                            for (let d = 0; d < 2; d++) {
                                if (!k._ds(e - d, g)) {
                                    h++;
                                    p <<= 1;
                                    if (this._au._ds(e - d, g)) {
                                        p |= 1
                                    }
                                    if (h == 8) {
                                        r[m++] = p;
                                        h = 0;
                                        p = 0
                                    }
                                }
                            }
                        }
                        // @ts-ignore
                        n ^= true
                    }
                    if (m != o._dp) {
                        throw"Error _gk"
                    }
                    return r
                }
            };

            let _dx: any = {};
            _dx._gl = function (a) {
                if (a < 0 || a > 7) {
                    throw"bad arguments"
                }
                return _dx._dy[a]
            };

            const _fg: any = () => {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return ((b + a) & 1) == 0
                }
            };

            const _fh:any=() =>{
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return (b & 1) == 0
                }
            };

            const _fi:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return a % 3 == 0
                }
            };

            const _fj:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return (b + a) % 3 == 0
                }
            };

            const _fk:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return (((_ew(b, 1)) + (a / 3)) & 1) == 0
                }
            };

            const _fl:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (c, b) {
                    let a = c * b;
                    return (a & 1) + (a % 3) == 0
                }
            };

            const _fm:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (c, b) {
                    let a = c * b;
                    return (((a & 1) + (a % 3)) & 1) == 0
                }
            };

            const _fn:any=()=> {
                this._dj = function (c, d) {
                    for (let b = 0; b < d; b++) {
                        for (let a = 0; a < d; a++) {
                            if (this._fw(b, a)) {
                                c.flip(a, b)
                            }
                        }
                    }
                };
                this._fw = function (b, a) {
                    return ((((b + a) & 1) + ((b * a) % 3)) & 1) == 0
                }
            };

            _dx._dy = [new _fg(), new _fh(), new _fi(), new _fj(), new _fk(), new _fl(), new _fm(), new _fn()];

            const _db:any=(a)=>{
                this._fa = a;
                this.decode = function (j, f) {
                    let c = new _bp(this._fa, j);
                    let p = new Array(f);
                    for (let g = 0; g < p.length; g++) {
                        p[g] = 0
                    }
                    let m = false;
                    let d = true;
                    for (let g = 0; g < f; g++) {
                        let q = c.evaluateAt(this._fa.exp(m ? g + 1 : g));
                        p[p.length - 1 - g] = q;
                        if (q != 0) {
                            d = false
                        }
                    }
                    if (d) {
                        return
                    }
                    let b = new _bp(this._fa, p);
                    let l = this._eb(this._fa._ba(f, 1), b, f);
                    let o = l[0];
                    let n = l[1];
                    let k = this._ey(o);
                    let e = this._di(n, k, m);
                    for (let g = 0; g < k.length; g++) {
                        let h = j.length - 1 - this._fa.log(k[g]);
                        if (h < 0) {
                            throw"ReedSolomonException Bad error location"
                        }
                        j[h] = _az._bd(j[h], e[g])
                    }
                };
                this._eb = function (z, y, f) {
                    if (z._ec < y._ec) {
                        let w = z;
                        z = y;
                        y = w
                    }
                    let B = z;
                    let k = y;
                    let o = this._fa.One;
                    let j = this._fa.Zero;
                    let e = this._fa.Zero;
                    let i = this._fa.One;
                    while (k._ec >= Math.floor(f / 2)) {
                        let x = B;
                        let g = o;
                        let v = e;
                        B = k;
                        o = j;
                        e = i;
                        if (B.Zero) {
                            throw"r_{i-1} was zero"
                        }
                        k = x;
                        let m = this._fa.Zero;
                        let p = B._ex(B._ec);
                        let h = this._fa.inverse(p);
                        while (k._ec >= B._ec && !k.Zero) {
                            let c = k._ec - B._ec;
                            let A = this._fa.multiply(k._ex(k._ec), h);
                            m = m._bd(this._fa._ba(c, A));
                            k = k._bd(B._dc(c, A))
                        }
                        j = m.multiply1(o)._bd(g);
                        i = m.multiply1(e)._bd(v)
                    }
                    let u = i._ex(0);
                    if (u == 0) {
                        throw"ReedSolomonException sigmaTilde(0) was zero"
                    }
                    let d = this._fa.inverse(u);
                    let n = i.multiply2(d);
                    let l = k.multiply2(d);
                    return [n, l]
                };
                this._ey = function (f) {
                    let g = f._ec;
                    if (g == 1) {
                        return new Array(f._ex(1))
                    }
                    let b = new Array(g);
                    let d = 0;
                    for (let c = 1; c < 256 && d < g; c++) {
                        if (f.evaluateAt(c) == 0) {
                            b[d] = this._fa.inverse(c);
                            d++
                        }
                    }
                    if (d != g) {
                        throw"Error locator degree does not match number of roots"
                    }
                    return b
                };
                this._di = function (f, h, g) {
                    let k = h.length;
                    let l = new Array(k);
                    for (let e = 0; e < k; e++) {
                        let b = this._fa.inverse(h[e]);
                        let c = 1;
                        for (let d = 0; d < k; d++) {
                            if (e != d) {
                                c = this._fa.multiply(c, _az._bd(1, this._fa.multiply(h[d], b)))
                            }
                        }
                        l[e] = this._fa.multiply(f.evaluateAt(b), this._fa.inverse(c));
                        if (g) {
                            l[e] = this._fa.multiply(l[e], b)
                        }
                    }
                    return l
                }
            };

            const _bp:any=(f, e)=>{
                if (e == null || e.length == 0) {
                    throw"bad arguments"
                }
                this._fa = f;
                let c = e.length;
                if (c > 1 && e[0] == 0) {
                    let d = 1;
                    while (d < c && e[d] == 0) {
                        d++
                    }
                    if (d == c) {
                        this._dd = f.Zero._dd
                    } else {
                        this._dd = new Array(c - d);
                        for (let b = 0; b < this._dd.length; b++) {
                            this._dd[b] = 0
                        }
                        for (let a = 0; a < this._dd.length; a++) {
                            this._dd[a] = e[d + a]
                        }
                    }
                } else {
                    this._dd = e
                }
                this.__defineGetter__("Zero", function () {
                    return this._dd[0] == 0
                });
                this.__defineGetter__("_ec", function () {
                    return this._dd.length - 1
                });
                this.__defineGetter__("Coefficients", function () {
                    return this._dd
                });
                this._ex = function (g) {
                    return this._dd[this._dd.length - 1 - g]
                };
                this.evaluateAt = function (h) {
                    if (h == 0) {
                        return this._ex(0)
                    }
                    let l = this._dd.length;
                    if (h == 1) {
                        let g = 0;
                        for (let k = 0; k < l; k++) {
                            g = _az._bd(g, this._dd[k])
                        }
                        return g
                    }
                    let j = this._dd[0];
                    for (let k = 1; k < l; k++) {
                        j = _az._bd(this._fa.multiply(h, j), this._dd[k])
                    }
                    return j
                };
                this._bd = function (g) {
                    if (this._fa != g._fa) {
                        throw"GF256Polys do not have same _az _fa"
                    }
                    if (this.Zero) {
                        return g
                    }
                    if (g.Zero) {
                        return this
                    }
                    let o = this._dd;
                    let n = g._dd;
                    if (o.length > n.length) {
                        let j = o;
                        o = n;
                        n = j
                    }
                    let h = new Array(n.length);
                    let k = n.length - o.length;
                    for (let m = 0; m < k; m++) {
                        h[m] = n[m]
                    }
                    for (let l = k; l < n.length; l++) {
                        h[l] = _az._bd(o[l - k], n[l])
                    }
                    return new _bp(f, h)
                };
                this.multiply1 = function (o) {
                    if (this._fa != o._fa) {
                        throw"GF256Polys do not have same _az _fa"
                    }
                    if (this.Zero || o.Zero) {
                        return this._fa.Zero
                    }
                    let q = this._dd;
                    let g = q.length;
                    let l = o._dd;
                    let n = l.length;
                    let p = new Array(g + n - 1);
                    for (let m = 0; m < g; m++) {
                        let h = q[m];
                        for (let k = 0; k < n; k++) {
                            p[m + k] = _az._bd(p[m + k], this._fa.multiply(h, l[k]))
                        }
                    }
                    return new _bp(this._fa, p)
                };
                this.multiply2 = function (g) {
                    if (g == 0) {
                        return this._fa.Zero
                    }
                    if (g == 1) {
                        return this
                    }
                    let j = this._dd.length;
                    let k = new Array(j);
                    for (let h = 0; h < j; h++) {
                        k[h] = this._fa.multiply(this._dd[h], g)
                    }
                    return new _bp(this._fa, k)
                };
                this._dc = function (l, g) {
                    if (l < 0) {
                        throw"bad arguments"
                    }
                    if (g == 0) {
                        return this._fa.Zero
                    }
                    let j = this._dd.length;
                    let k = new Array(j + l);
                    for (let h = 0; h < k.length; h++) {
                        k[h] = 0
                    }
                    for (let h = 0; h < j; h++) {
                        k[h] = this._fa.multiply(this._dd[h], g)
                    }
                    return new _bp(this._fa, k)
                };
                this.divide = function (l) {
                    if (this._fa != l._fa) {
                        throw"GF256Polys do not have same _az _fa"
                    }
                    if (l.Zero) {
                        throw"Divide by 0"
                    }
                    let j = this._fa.Zero;
                    let o = this;
                    let g = l._ex(l._ec);
                    let n = this._fa.inverse(g);
                    while (o._ec >= l._ec && !o.Zero) {
                        let m = o._ec - l._ec;
                        let h = this._fa.multiply(o._ex(o._ec), n);
                        let i = l._dc(m, h);
                        let k = this._fa._ba(m, h);
                        j = j._bd(k);
                        o = o._bd(i)
                    }
                    return [j, o]
                }
            };

            const _az:any=(b)=> {
                this._gh = new Array(256);
                this._gi = new Array(256);
                let a = 1;
                for (let e = 0; e < 256; e++) {
                    this._gh[e] = a;
                    a <<= 1;
                    if (a >= 256) {
                        a ^= b
                    }
                }
                for (let e = 0; e < 255; e++) {
                    this._gi[this._gh[e]] = e
                }
                let d = new Array(1);
                d[0] = 0;
                this.zero = new _bp(this, new Array(d));
                let c = new Array(1);
                c[0] = 1;
                this.one = new _bp(this, new Array(c));
                this.__defineGetter__("Zero", function () {
                    return this.zero
                });
                this.__defineGetter__("One", function () {
                    return this.one
                });
                this._ba = function (j, f) {
                    if (j < 0) {
                        throw"bad arguments"
                    }
                    if (f == 0) {
                        let zero;
                        return zero
                    }
                    let h = new Array(j + 1);
                    for (let g = 0; g < h.length; g++) {
                        h[g] = 0
                    }
                    h[0] = f;
                    return new _bp(this, h)
                };
                this.exp = function (f) {
                    return this._gh[f]
                };
                this.log = function (f) {
                    if (f == 0) {
                        throw"bad arguments"
                    }
                    return this._gi[f]
                };
                this.inverse = function (f) {
                    if (f == 0) {
                        throw"System.ArithmeticException"
                    }
                    return this._gh[255 - this._gi[f]]
                };
                this.multiply = function (g, f) {
                    if (g == 0 || f == 0) {
                        return 0
                    }
                    if (g == 1) {
                        return f
                    }
                    if (f == 1) {
                        return g
                    }
                    return this._gh[(this._gi[g] + this._gi[f]) % 255]
                }
            };

            _az._bb = new _az(285);
            _az._bc = new _az(301);
            _az._bd = function (d, c) {
                return d ^ c
            };
            let Decoder = {
                rsDecoder: undefined, decode: undefined,

                correctErrors(d: any, g: any) {
                    
                }
            };
            Decoder.rsDecoder = new _db(_az._bb);
            Decoder.correctErrors = function (g, b) {
                let d = g.length;
                let f = new Array(d);
                for (let e = 0; e < d; e++) {
                    f[e] = g[e] & 255
                }
                let a = g.length - b;
                try {
                    Decoder.rsDecoder.decode(f, a)
                } catch (c) {
                    throw c
                }
                for (let e = 0; e < b; e++) {
                    g[e] = f[e]
                }
            };
            Decoder.decode = function (q) {
                let b = new _cl(q);
                let o = b._cq();
                let c = b._cm()._cg;
                let p = b._gk();
                let a = _dl._gn(p, o, c);
                let f = 0;
                for (let k = 0; k < a.length; k++) {
                    f += a[k]._du
                }
                let e = new Array(f);
                let n = 0;
                for (let h = 0; h < a.length; h++) {
                    let m = a[h];
                    let d = m.Codewords;
                    let g = m._du;
                    Decoder.correctErrors(d, g);
                    for (let k = 0; k < g; k++) {
                        e[n++] = d[k]
                    }
                }
                let l = new QRCodeDataBlockReader(e, o._fd, c.Bits);
                return l
            };
            let qrcode = {
                binarize: undefined,
                width: undefined,
                imagedata: undefined,
                height: undefined,
                maxImgSize: undefined,
                debug: undefined,
                qrCodeSymbol: undefined,
                callback: undefined,
                _eo: undefined,
                vidSuccess: undefined,
                localstream: undefined,
                webkit: false,
                video: undefined,
                moz: false,
                canvas_qr2: undefined,
                gUM: undefined,
                qrcontext2: undefined,
                vidError: undefined,
                setWebcam: undefined, result: undefined, isUrl: undefined, getPixel: undefined, decode_url: undefined,
                decode_utf8: undefined, grayscale: undefined, _er: undefined,


                captureToCanvas: () => {
                },
                decode: undefined,
                process(a: any) {
                    return undefined;
                },
                grayScaleToBitmap(grayscale: any) {

                },
                _em(f: any) {

                }
            };
            qrcode.imagedata = null;
            qrcode.width = 0;
            qrcode.height = 0;
            qrcode.qrCodeSymbol = null;
            qrcode.debug = false;
            qrcode.maxImgSize = 1024 * 1024;
            qrcode._eo = [[10, 9, 8, 8], [12, 11, 16, 10], [14, 13, 16, 12]];
            qrcode.callback = null;
            qrcode.vidSuccess = function (a) {
                qrcode.localstream = a;
                if (qrcode.webkit) {
                    qrcode.video.src = (window as any).webkitURL.createObjectURL(a)
                } else {
                    if (qrcode.moz) {
                        qrcode.video.mozSrcObject = a;
                        qrcode.video.play()
                    } else {
                        qrcode.video.src = a
                    }
                }
                qrcode.gUM = true;
                qrcode.canvas_qr2 = document.createElement("canvas");
                qrcode.canvas_qr2.id = "QRScanner-canvasEl";
                qrcode.qrcontext2 = qrcode.canvas_qr2.getContext("2d");
                qrcode.canvas_qr2.width = qrcode.video.videoWidth;
                qrcode.canvas_qr2.height = qrcode.video.videoHeight;
                setTimeout(qrcode.captureToCanvas, 500)
            };
            qrcode.vidError = function (a) {
                qrcode.gUM = false;
            };
            qrcode.captureToCanvas = function () {
                if (qrcode.gUM) {
                    try {
                        if (qrcode.video.videoWidth == 0) {
                            setTimeout(qrcode.captureToCanvas, 500);
                            return
                        } else {
                            qrcode.canvas_qr2.width = qrcode.video.videoWidth;
                            qrcode.canvas_qr2.height = qrcode.video.videoHeight
                        }
                        qrcode.qrcontext2.drawImage(qrcode.video, 0, 0);
                        try {
                            qrcode.decode()
                        } catch (a) {
                            setTimeout(qrcode.captureToCanvas, 500)
                        }
                    } catch (a) {
                        setTimeout(qrcode.captureToCanvas, 500)
                    }
                }
            };
            qrcode.setWebcam = function (c) {
                let d: any = navigator;
                qrcode.video = document.getElementById(c);
                let a: any = true;
                if (navigator.mediaDevices && navigator.mediaDevices.enumerateDevices) {
                    try {
                        navigator.mediaDevices.enumerateDevices().then(function (e) {
                            e.forEach(function (f) {
                                if (f.kind === "videoinput") {
                                    if (f.label.toLowerCase().search("back") > -1) {
                                        a = [{sourceId: f.deviceId}]
                                    }
                                }
                            })
                        })
                    } catch (b) {
                    }
                } else {
                    console.log("no navigator.mediaDevices.enumerateDevices")
                }
                if (d.getUserMedia) {
                    d.getUserMedia({video: a, audio: false}, qrcode.vidSuccess, qrcode.vidError)
                } else {
                    if (d.webkitGetUserMedia) {
                        qrcode.webkit = true;
                        d.webkitGetUserMedia({video: a, audio: false}, qrcode.vidSuccess, qrcode.vidError)
                    } else {
                        if (d.mozGetUserMedia) {
                            qrcode.moz = true;
                            d.mozGetUserMedia({video: a, audio: false}, qrcode.vidSuccess, qrcode.vidError)
                        }
                    }
                }
            };
            qrcode.decode = function (d) {
                if (arguments.length == 0) {
                    let a: any;
                    let b: any;
                    if (qrcode.canvas_qr2) {
                        b = qrcode.canvas_qr2;
                        a = qrcode.qrcontext2
                    } else {
                        b = document.getElementById("QRScanner-canvasEl") || document.getElementById("qr-canvas");
                        a = b.getContext("2d")
                    }
                    qrcode.width = b.width;
                    qrcode.height = b.height;
                    qrcode.imagedata = a.getImageData(0, 0, qrcode.width, qrcode.height);
                    qrcode.result = qrcode.process(a);
                    if (qrcode.callback != null) {
                        qrcode.callback(qrcode.result)
                    }
                    return qrcode.result
                } else {
                    let c = new Image();
                    c.crossOrigin = "Anonymous";
                    c.onload = function () {
                        let g: any = document.getElementById("out-canvas");
                        if (g != null) {
                            let j = g.getContext("2d");
                            j.clearRect(0, 0, 320, 240);
                            j.drawImage(c, 0, 0, 320, 240)
                        }
                        let i = document.createElement("canvas");
                        let h = i.getContext("2d");
                        let f = c.height;
                        let l = c.width;
                        if (c.width * c.height > qrcode.maxImgSize) {
                            let k = c.width / c.height;
                            f = Math.sqrt(qrcode.maxImgSize / k);
                            l = k * f
                        }
                        i.width = l;
                        i.height = f;
                        h.drawImage(c, 0, 0, i.width, i.height);
                        qrcode.width = i.width;
                        qrcode.height = i.height;
                        try {
                            qrcode.imagedata = h.getImageData(0, 0, i.width, i.height)
                        } catch (m) {
                            qrcode.result = "Cross domain image reading not supported in your browser! Save it to your computer then drag and drop the file!";
                            if (qrcode.callback != null) {
                                qrcode.callback(qrcode.result)
                            }
                            return
                        }
                        try {
                            qrcode.result = qrcode.process(h)
                        } catch (m) {
                            qrcode.result = "error decoding QR Code"
                        }
                        if (qrcode.callback != null) {
                            qrcode.callback(qrcode.result)
                        }
                    };
                    c.onerror = function () {
                        if (qrcode.callback != null) {
                            qrcode.callback("Failed to load the image")
                        }
                    };
                    c.src = d
                }
            };
            qrcode.isUrl = function (a) {
                let b = /(ftp|http|https):\/\/(\w+:{0,1}\w*@)?(\S+)(:[0-9]+)?(\/|\/([\w#!:.?+=&%@!\-\/]))?/;
                return b.test(a)
            };
            qrcode.decode_url = function (b) {
                let d = "";
                try {
                    d = escape(b)
                } catch (c) {
                    d = b
                }
                let a = "";
                try {
                    a = decodeURIComponent(d)
                } catch (c) {
                    a = d
                }
                return a
            };
            qrcode.decode_utf8 = function (a) {
                if (qrcode.isUrl(a)) {
                    return qrcode.decode_url(a)
                } else {
                    return a
                }
            };
            qrcode.process = function (q) {
                let a = new Date().getTime();
                let c = qrcode.grayScaleToBitmap(qrcode.grayscale());
                if (qrcode.debug) {
                    for (let m = 0; m < qrcode.height; m++) {
                        for (let n = 0; n < qrcode.width; n++) {
                            let o = (n * 4) + (m * qrcode.width * 4);
                            qrcode.imagedata.data[o] = c[n + m * qrcode.width] ? 0 : 0;
                            qrcode.imagedata.data[o + 1] = c[n + m * qrcode.width] ? 0 : 0;
                            qrcode.imagedata.data[o + 2] = c[n + m * qrcode.width] ? 255 : 0
                        }
                    }
                    q.putImageData(qrcode.imagedata, 0, 0)
                }
                let h = new Detector(c);
                let p = h.detect();
                if (qrcode.debug) {
                    for (let m = 0; m < p.bits.Height; m++) {
                        for (let n = 0; n < p.bits.Width; n++) {
                            let o = (n * 4 * 2) + (m * 2 * qrcode.width * 4);
                            qrcode.imagedata.data[o] = p.bits._ds(n, m) ? 0 : 0;
                            qrcode.imagedata.data[o + 1] = p.bits._ds(n, m) ? 0 : 0;
                            qrcode.imagedata.data[o + 2] = p.bits._ds(n, m) ? 255 : 0
                        }
                    }
                    q.putImageData(qrcode.imagedata, 0, 0)
                }
                let k = Decoder.decode(p.bits);
                let g = k.DataByte;
                let l = "";
                for (let f = 0; f < g.length; f++) {
                    for (let e = 0; e < g[f].length; e++) {
                        l += String.fromCharCode(g[f][e])
                    }
                }
                let d = new Date().getTime();
                let b = d - a;
                return qrcode.decode_utf8(l)
            };
            qrcode.getPixel = function (b, d) {
                if (qrcode.width < b) {
                    throw"point error"
                }
                if (qrcode.height < d) {
                    throw"point error"
                }
                let a = (b * 4) + (d * qrcode.width * 4);
                let c = (qrcode.imagedata.data[a] * 33 + qrcode.imagedata.data[a + 1] * 34 + qrcode.imagedata.data[a + 2] * 33) / 100;
                return c
            };
            qrcode.binarize = function (d) {
                let c = new Array(qrcode.width * qrcode.height);
                for (let e = 0; e < qrcode.height; e++) {
                    for (let b = 0; b < qrcode.width; b++) {
                        let a = qrcode.getPixel(b, e);
                        c[b + e * qrcode.width] = a <= d ? true : false
                    }
                }
                return c
            };
            qrcode._em = function (d) {
                let c = 4;
                let k = Math.floor(qrcode.width / c);
                let j = Math.floor(qrcode.height / c);
                let f = new Array(c);
                for (let g = 0; g < c; g++) {
                    f[g] = new Array(c);
                    for (let e = 0; e < c; e++) {
                        f[g][e] = [0, 0]
                    }
                }
                for (let o = 0; o < c; o++) {
                    for (let a = 0; a < c; a++) {
                        f[a][o][0] = 255;
                        for (let l = 0; l < j; l++) {
                            for (let n = 0; n < k; n++) {
                                let h = d[k * a + n + (j * o + l) * qrcode.width];
                                if (h < f[a][o][0]) {
                                    f[a][o][0] = h
                                }
                                if (h > f[a][o][1]) {
                                    f[a][o][1] = h
                                }
                            }
                        }
                    }
                }
                let m = new Array(c);
                for (let b = 0; b < c; b++) {
                    m[b] = new Array(c)
                }
                for (let o = 0; o < c; o++) {
                    for (let a = 0; a < c; a++) {
                        m[a][o] = Math.floor((f[a][o][0] + f[a][o][1]) / 2)
                    }
                }
                return m
            };
            qrcode.grayScaleToBitmap = function (f) {
                let k: any = qrcode._em(f);
                let b = k.length;
                let e = Math.floor(qrcode.width / b);
                let d = Math.floor(qrcode.height / b);
                let h = new ArrayBuffer(qrcode.width * qrcode.height);
                let c:any = new Uint8Array(h);
                for (let j = 0; j < b; j++) {
                    for (let a = 0; a < b; a++) {
                        for (let g = 0; g < d; g++) {
                            for (let i = 0; i < e; i++) {
                                c[e * a + i + (d * j + g) * qrcode.width] = (f[e * a + i + (d * j + g) * qrcode.width] < k[a][j]) ? true : false
                            }
                        }
                    }
                }
                return c
            };
            qrcode.grayscale = function () {
                let e = new ArrayBuffer(qrcode.width * qrcode.height);
                let c = new Uint8Array(e);
                for (let d = 0; d < qrcode.height; d++) {
                    for (let b = 0; b < qrcode.width; b++) {
                        let a = qrcode.getPixel(b, d);
                        c[b + d * qrcode.width] = a
                    }
                }
                return c
            };

            const _ew:any=(a, b)=>{
                if (a >= 0) {
                    return a >> b
                } else {
                    return (a >> b) + (2 << ~b)
                }
            };

            let _gf = 3;
            let _eh = 57;
            let _el = 8;
            let _eg = 2;
            qrcode._er = function (c) {
                function b(m, k) {
                    let n = m.X - k.X;
                    let l = m.Y - k.Y;
                    return Math.sqrt((n * n + l * l))
                }

                function d(k, o, n) {
                    let m = o.x;
                    let l = o.y;
                    return ((n.x - m) * (k.y - l)) - ((n.y - l) * (k.x - m))
                }

                let i = b(c[0], c[1]);
                let f = b(c[1], c[2]);
                let e = b(c[0], c[2]);
                let a, j, h;
                if (f >= i && f >= e) {
                    j = c[0];
                    a = c[1];
                    h = c[2]
                } else {
                    if (e >= f && e >= i) {
                        j = c[1];
                        a = c[0];
                        h = c[2]
                    } else {
                        j = c[2];
                        a = c[0];
                        h = c[1]
                    }
                }
                if (d(a, j, h) < 0) {
                    let g = a;
                    a = h;
                    h = g
                }
                c[0] = a;
                c[1] = j;
                c[2] = h
            };

            const _cz:any=(c, a, b)=> {
                this.x = c;
                this.y = a;
                this.count = 1;
                this._aj = b;
                this.__defineGetter__("_ei", function () {
                    return this._aj
                });
                this.__defineGetter__("Count", function () {
                    return this.count
                });
                this.__defineGetter__("X", function () {
                    return this.x
                });
                this.__defineGetter__("Y", function () {
                    return this.y
                });
                this._ek = function () {
                    this.count++
                };
                this._ev = function (f, e, d) {
                    if (Math.abs(e - this.y) <= f && Math.abs(d - this.x) <= f) {
                        let g = Math.abs(f - this._aj);
                        return g <= 1 || g / this._aj <= 1
                    }
                    return false
                }
            };

            const _es:any=(a)=> {
                this._go = a[0];
                this._gu = a[1];
                this._gr = a[2];
                this.__defineGetter__("_gp", function () {
                    return this._go
                });
                this.__defineGetter__("_gq", function () {
                    return this._gu
                });
                this.__defineGetter__("_gs", function () {
                    return this._gr
                })
            };

            const _cc:any=()=> {
                this.image = null;
                this._cv = [];
                this._ge = false;
                this._al = [0, 0, 0, 0, 0];
                this._am = null;
                this.__defineGetter__("_da", function () {
                    this._al[0] = 0;
                    this._al[1] = 0;
                    this._al[2] = 0;
                    this._al[3] = 0;
                    this._al[4] = 0;
                    return this._al
                });
                this._ao = function (f) {
                    let b = 0;
                    for (let d = 0; d < 5; d++) {
                        let e = f[d];
                        if (e == 0) {
                            return false
                        }
                        b += e
                    }
                    if (b < 7) {
                        return false
                    }
                    let c = Math.floor((b << _el) / 7);
                    let a = Math.floor(c / 2);
                    return Math.abs(c - (f[0] << _el)) < a && Math.abs(c - (f[1] << _el)) < a && Math.abs(3 * c - (f[2] << _el)) < 3 * a && Math.abs(c - (f[3] << _el)) < a && Math.abs(c - (f[4] << _el)) < a
                };
                this._an = function (b, a) {
                    return (a - b[4] - b[3]) - b[2] / 2
                };
                this._ap = function (a, j, d, g) {
                    let c = this.image;
                    let h = qrcode.height;
                    let b = this._da;
                    let f = a;
                    while (f >= 0 && c[j + f * qrcode.width]) {
                        b[2]++;
                        f--
                    }
                    if (f < 0) {
                        return NaN
                    }
                    while (f >= 0 && !c[j + f * qrcode.width] && b[1] <= d) {
                        b[1]++;
                        f--
                    }
                    if (f < 0 || b[1] > d) {
                        return NaN
                    }
                    while (f >= 0 && c[j + f * qrcode.width] && b[0] <= d) {
                        b[0]++;
                        f--
                    }
                    if (b[0] > d) {
                        return NaN
                    }
                    f = a + 1;
                    while (f < h && c[j + f * qrcode.width]) {
                        b[2]++;
                        f++
                    }
                    if (f == h) {
                        return NaN
                    }
                    while (f < h && !c[j + f * qrcode.width] && b[3] < d) {
                        b[3]++;
                        f++
                    }
                    if (f == h || b[3] >= d) {
                        return NaN
                    }
                    while (f < h && c[j + f * qrcode.width] && b[4] < d) {
                        b[4]++;
                        f++
                    }
                    if (b[4] >= d) {
                        return NaN
                    }
                    let e = b[0] + b[1] + b[2] + b[3] + b[4];
                    if (5 * Math.abs(e - g) >= 2 * g) {
                        return NaN
                    }
                    return this._ao(b) ? this._an(b, f) : NaN
                };
                this._ej = function (b, a, e, h) {
                    let d = this.image;
                    let i = qrcode.width;
                    let c = this._da;
                    let g = b;
                    while (g >= 0 && d[g + a * qrcode.width]) {
                        c[2]++;
                        g--
                    }
                    if (g < 0) {
                        return NaN
                    }
                    while (g >= 0 && !d[g + a * qrcode.width] && c[1] <= e) {
                        c[1]++;
                        g--
                    }
                    if (g < 0 || c[1] > e) {
                        return NaN
                    }
                    while (g >= 0 && d[g + a * qrcode.width] && c[0] <= e) {
                        c[0]++;
                        g--
                    }
                    if (c[0] > e) {
                        return NaN
                    }
                    g = b + 1;
                    while (g < i && d[g + a * qrcode.width]) {
                        c[2]++;
                        g++
                    }
                    if (g == i) {
                        return NaN
                    }
                    while (g < i && !d[g + a * qrcode.width] && c[3] < e) {
                        c[3]++;
                        g++
                    }
                    if (g == i || c[3] >= e) {
                        return NaN
                    }
                    while (g < i && d[g + a * qrcode.width] && c[4] < e) {
                        c[4]++;
                        g++
                    }
                    if (c[4] >= e) {
                        return NaN
                    }
                    let f = c[0] + c[1] + c[2] + c[3] + c[4];
                    if (5 * Math.abs(f - h) >= h) {
                        return NaN
                    }
                    return this._ao(c) ? this._an(c, g) : NaN
                };
                this._cu = function (c, f, e) {
                    let d = c[0] + c[1] + c[2] + c[3] + c[4];
                    let n = this._an(c, e);
                    let b = this._ap(f, Math.floor(n), c[2], d);
                    if (!isNaN(b)) {
                        n = this._ej(Math.floor(n), Math.floor(b), c[2], d);
                        if (!isNaN(n)) {
                            let l = d / 7;
                            let m = false;
                            let h = this._cv.length;
                            for (let g = 0; g < h; g++) {
                                let a = this._cv[g];
                                if (a._ev(l, b, n)) {
                                    a._ek();
                                    m = true;
                                    break
                                }
                            }
                            if (!m) {
                                let k = new _cz(n, b, l);
                                this._cv.push(k);
                                if (this._am != null) {
                                    this._am._ep(k)
                                }
                            }
                            return true
                        }
                    }
                    return false
                };
                this._ee = function () {
                    let h = this._cv.length;
                    if (h < 3) {
                        throw"Couldn't find enough finder patterns (found " + h + ")"
                    }
                    if (h > 3) {
                        let b = 0;
                        let j = 0;
                        for (let d = 0; d < h; d++) {
                            let g = this._cv[d]._ei;
                            b += g;
                            j += (g * g)
                        }
                        let a = b / h;
                        this._cv.sort(function (m, l) {
                            let k = Math.abs(l._ei - a);
                            let i = Math.abs(m._ei - a);
                            if (k < i) {
                                return (-1)
                            } else {
                                if (k == i) {
                                    return 0
                                } else {
                                    return 1
                                }
                            }
                        });
                        let e = Math.sqrt(j / h - a * a);
                        let c = Math.max(0.2 * a, e);
                        for (let d = this._cv.length - 1; d >= 0; d--) {
                            let f = this._cv[d];
                            if (Math.abs(f._ei - a) > c) {
                                this._cv.splice(d, 1)
                            }
                        }
                    }
                    if (this._cv.length > 3) {
                        this._cv.sort(function (k, i) {
                            if (k.count > i.count) {
                                return -1
                            }
                            if (k.count < i.count) {
                                return 1
                            }
                            return 0
                        })
                    }
                    return [this._cv[0], this._cv[1], this._cv[2]]
                };
                this._eq = function () {
                    let b = this._cv.length;
                    if (b <= 1) {
                        return 0
                    }
                    let c = null;
                    for (let d = 0; d < b; d++) {
                        let a = this._cv[d];
                        if (a.Count >= _eg) {
                            if (c == null) {
                                c = a
                            } else {
                                this._ge = true;
                                return Math.floor((Math.abs(c.X - a.X) - Math.abs(c.Y - a.Y)) / 2)
                            }
                        }
                    }
                    return 0
                };
                this._cx = function () {
                    let g = 0;
                    let c = 0;
                    let a = this._cv.length;
                    for (let d = 0; d < a; d++) {
                        let f = this._cv[d];
                        if (f.Count >= _eg) {
                            g++;
                            c += f._ei
                        }
                    }
                    if (g < 3) {
                        return false
                    }
                    let e = c / a;
                    let b = 0;
                    for (let d = 0; d < a; d++) {
                        let f = this._cv[d];
                        b += Math.abs(f._ei - e)
                    }
                    return b <= 0.05 * c
                };
                this._ce = function (e) {
                    let o = false;
                    this.image = e;
                    let n = qrcode.height;
                    let k = qrcode.width;
                    let a = Math.floor((3 * n) / (4 * _eh));
                    if (a < _gf || o) {
                        a = _gf
                    }
                    let g = false;
                    let d = new Array(5);
                    for (let h = a - 1; h < n && !g; h += a) {
                        d[0] = 0;
                        d[1] = 0;
                        d[2] = 0;
                        d[3] = 0;
                        d[4] = 0;
                        let b = 0;
                        for (let f = 0; f < k; f++) {
                            if (e[f + h * qrcode.width]) {
                                if ((b & 1) == 1) {
                                    b++
                                }
                                d[b]++
                            } else {
                                if ((b & 1) == 0) {
                                    if (b == 4) {
                                        if (this._ao(d)) {
                                            let c = this._cu(d, h, f);
                                            if (c) {
                                                a = 2;
                                                if (this._ge) {
                                                    g = this._cx()
                                                } else {
                                                    let m = this._eq();
                                                    if (m > d[2]) {
                                                        h += m - d[2] - a;
                                                        f = k - 1
                                                    }
                                                }
                                            } else {
                                                do {
                                                    f++
                                                } while (f < k && !e[f + h * qrcode.width]);
                                                f--
                                            }
                                            b = 0;
                                            d[0] = 0;
                                            d[1] = 0;
                                            d[2] = 0;
                                            d[3] = 0;
                                            d[4] = 0
                                        } else {
                                            d[0] = d[2];
                                            d[1] = d[3];
                                            d[2] = d[4];
                                            d[3] = 1;
                                            d[4] = 0;
                                            b = 3
                                        }
                                    } else {
                                        d[++b]++
                                    }
                                } else {
                                    d[b]++
                                }
                            }
                        }
                        if (this._ao(d)) {
                            let c = this._cu(d, h, k);
                            if (c) {
                                a = d[0];
                                if (this._ge) {
                                    g = this._cx()
                                }
                            }
                        }
                    }
                    let l = this._ee();
                    qrcode._er(l);
                    return new _es(l)
                }
            };

            const _ai:any=(c, a, b)=>{
                this.x = c;
                this.y = a;
                this.count = 1;
                this._aj = b;
                this.__defineGetter__("_ei", function () {
                    return this._aj
                });
                this.__defineGetter__("Count", function () {
                    return this.count
                });
                this.__defineGetter__("X", function () {
                    return Math.floor(this.x)
                });
                this.__defineGetter__("Y", function () {
                    return Math.floor(this.y)
                });
                this._ek = function () {
                    this.count++
                };
                this._ev = function (f, e, d) {
                    if (Math.abs(e - this.y) <= f && Math.abs(d - this.x) <= f) {
                        let g = Math.abs(f - this._aj);
                        return g <= 1 || g / this._aj <= 1
                    }
                    return false
                }
            };

            const _ak:any=(g, c, b, f, a, e, d)=>{
                this.image = g;
                this._cv = [];
                this.startX = c;
                this.startY = b;
                this.width = f;
                this.height = a;
                this._ef = e;
                this._al = [0, 0, 0];
                this._am = d;
                this._an = function (i, h) {
                    return (h - i[2]) - i[1] / 2
                };
                this._ao = function (l) {
                    let k = this._ef;
                    let h = k / 2;
                    for (let j = 0; j < 3; j++) {
                        if (Math.abs(k - l[j]) >= h) {
                            return false
                        }
                    }
                    return true
                };
                this._ap = function (h, q, l, o) {
                    let k = this.image;
                    let p = qrcode.height;
                    let j = this._al;
                    j[0] = 0;
                    j[1] = 0;
                    j[2] = 0;
                    let n = h;
                    while (n >= 0 && k[q + n * qrcode.width] && j[1] <= l) {
                        j[1]++;
                        n--
                    }
                    if (n < 0 || j[1] > l) {
                        return NaN
                    }
                    while (n >= 0 && !k[q + n * qrcode.width] && j[0] <= l) {
                        j[0]++;
                        n--
                    }
                    if (j[0] > l) {
                        return NaN
                    }
                    n = h + 1;
                    while (n < p && k[q + n * qrcode.width] && j[1] <= l) {
                        j[1]++;
                        n++
                    }
                    if (n == p || j[1] > l) {
                        return NaN
                    }
                    while (n < p && !k[q + n * qrcode.width] && j[2] <= l) {
                        j[2]++;
                        n++
                    }
                    if (j[2] > l) {
                        return NaN
                    }
                    let m = j[0] + j[1] + j[2];
                    if (5 * Math.abs(m - o) >= 2 * o) {
                        return NaN
                    }
                    return this._ao(j) ? this._an(j, n) : NaN
                };
                this._cu = function (l, o, n) {
                    let m = l[0] + l[1] + l[2];
                    let t = this._an(l, n);
                    let k = this._ap(o, Math.floor(t), 2 * l[1], m);
                    if (!isNaN(k)) {
                        let s = (l[0] + l[1] + l[2]) / 3;
                        let q = this._cv.length;
                        for (let p = 0; p < q; p++) {
                            let h = this._cv[p];
                            if (h._ev(s, k, t)) {
                                return new _ai(t, k, s)
                            }
                        }
                        let r = new _ai(t, k, s);
                        this._cv.push(r);
                        if (this._am != null) {
                            this._am._ep(r)
                        }
                    }
                    return null
                };
                this.find = function () {
                    let p = this.startX;
                    let s = this.height;
                    let q = p + f;
                    let r = b + (s >> 1);
                    let m = [0, 0, 0];
                    for (let k = 0; k < s; k++) {
                        let o = r + ((k & 1) == 0 ? ((k + 1) >> 1) : -((k + 1) >> 1));
                        m[0] = 0;
                        m[1] = 0;
                        m[2] = 0;
                        let n = p;
                        while (n < q && !g[n + qrcode.width * o]) {
                            n++
                        }
                        let h = 0;
                        while (n < q) {
                            if (g[n + o * qrcode.width]) {
                                if (h == 1) {
                                    m[h]++
                                } else {
                                    if (h == 2) {
                                        if (this._ao(m)) {
                                            let l = this._cu(m, o, n);
                                            if (l != null) {
                                                return l
                                            }
                                        }
                                        m[0] = m[2];
                                        m[1] = 1;
                                        m[2] = 0;
                                        h = 1
                                    } else {
                                        m[++h]++
                                    }
                                }
                            } else {
                                if (h == 1) {
                                    h++
                                }
                                m[h]++
                            }
                            n++
                        }
                        if (this._ao(m)) {
                            let l = this._cu(m, o, q);
                            if (l != null) {
                                return l
                            }
                        }
                    }
                    if (!(this._cv.length == 0)) {
                        return this._cv[0]
                    }
                    throw"Couldn't find enough alignment patterns"
                }
            };

            const QRCodeDataBlockReader:any=(c, a, b) =>{
                this._ed = 0;
                this._cw = 7;
                this.dataLength = 0;
                this.blocks = c;
                this._en = b;
                if (a <= 9) {
                    this.dataLengthMode = 0
                } else {
                    if (a >= 10 && a <= 26) {
                        this.dataLengthMode = 1
                    } else {
                        if (a >= 27 && a <= 40) {
                            this.dataLengthMode = 2
                        }
                    }
                }
                this._gd = function (f) {
                    let k = 0;
                    if (f < this._cw + 1) {
                        let m = 0;
                        for (let e = 0; e < f; e++) {
                            m += (1 << e)
                        }
                        m <<= (this._cw - f + 1);
                        k = (this.blocks[this._ed] & m) >> (this._cw - f + 1);
                        this._cw -= f;
                        return k
                    } else {
                        if (f < this._cw + 1 + 8) {
                            let j = 0;
                            for (let e = 0; e < this._cw + 1; e++) {
                                j += (1 << e)
                            }
                            k = (this.blocks[this._ed] & j) << (f - (this._cw + 1));
                            this._ed++;
                            k += ((this.blocks[this._ed]) >> (8 - (f - (this._cw + 1))));
                            this._cw = this._cw - f % 8;
                            if (this._cw < 0) {
                                this._cw = 8 + this._cw
                            }
                            return k
                        } else {
                            if (f < this._cw + 1 + 16) {
                                let j = 0;
                                let h = 0;
                                for (let e = 0; e < this._cw + 1; e++) {
                                    j += (1 << e)
                                }
                                let g = (this.blocks[this._ed] & j) << (f - (this._cw + 1));
                                this._ed++;
                                let d = this.blocks[this._ed] << (f - (this._cw + 1 + 8));
                                this._ed++;
                                for (let e = 0; e < f - (this._cw + 1 + 8); e++) {
                                    h += (1 << e)
                                }
                                h <<= 8 - (f - (this._cw + 1 + 8));
                                let l = (this.blocks[this._ed] & h) >> (8 - (f - (this._cw + 1 + 8)));
                                k = g + d + l;
                                this._cw = this._cw - (f - 8) % 8;
                                if (this._cw < 0) {
                                    this._cw = 8 + this._cw
                                }
                                return k
                            } else {
                                return 0
                            }
                        }
                    }
                };
                this.NextMode = function () {
                    if ((this._ed > this.blocks.length - this._en - 2)) {
                        return 0
                    } else {
                        return this._gd(4)
                    }
                };
                this.getDataLength = function (d) {
                    let e = 0;
                    while (true) {
                        if ((d >> e) == 1) {
                            break
                        }
                        e++
                    }
                    return this._gd(qrcode._eo[this.dataLengthMode][e])
                };
                this.getRomanAndFigureString = function (h) {
                    let f = h;
                    let g = 0;
                    let j = "";
                    let d = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", " ", "$", "%", "*", "+", "-", ".", "/", ":"];
                    do {
                        if (f > 1) {
                            g = this._gd(11);
                            let i = Math.floor(g / 45);
                            let e = g % 45;
                            j += d[i];
                            j += d[e];
                            f -= 2
                        } else {
                            if (f == 1) {
                                g = this._gd(6);
                                j += d[g];
                                f -= 1
                            }
                        }
                    } while (f > 0);
                    return j
                };
                this.getFigureString = function (f) {
                    let d = f;
                    let e = 0;
                    let g = "";
                    do {
                        if (d >= 3) {
                            e = this._gd(10);
                            if (e < 100) {
                                g += "0"
                            }
                            if (e < 10) {
                                g += "0"
                            }
                            d -= 3
                        } else {
                            if (d == 2) {
                                e = this._gd(7);
                                if (e < 10) {
                                    g += "0"
                                }
                                d -= 2
                            } else {
                                if (d == 1) {
                                    e = this._gd(4);
                                    d -= 1
                                }
                            }
                        }
                        g += e
                    } while (d > 0);
                    return g
                };
                this.get8bitByteArray = function (g) {
                    let e = g;
                    let f = 0;
                    let d = [];
                    do {
                        f = this._gd(8);
                        d.push(f);
                        e--
                    } while (e > 0);
                    return d
                };
                this.getKanjiString = function (j) {
                    let g = j;
                    let i = 0;
                    let h = "";
                    do {
                        i = this._gd(13);
                        let e = i % 192;
                        let f = i / 192;
                        let k = (f << 8) + e;
                        let d = 0;
                        if (k + 33088 <= 40956) {
                            d = k + 33088
                        } else {
                            d = k + 49472
                        }
                        h += String.fromCharCode(d);
                        g--
                    } while (g > 0);
                    return h
                };
                this.__defineGetter__("DataByte", function () {
                    let h = [];
                    let e = 1;
                    let f = 2;
                    let d = 4;
                    let o = 8;
                    do {
                        let l = this.NextMode();
                        if (l == 0) {
                            if (h.length > 0) {
                                break
                            } else {
                                throw"Empty data block"
                            }
                        }
                        if (l != e && l != f && l != d && l != o) {
                            throw"Invalid mode: " + l + " in (block:" + this._ed + " bit:" + this._cw + ")"
                        }
                        let g = this.getDataLength(l);
                        if (g < 1) {
                            throw"Invalid data length: " + g
                        }
                        switch (l) {
                            case e:
                                let m = this.getFigureString(g);
                                let k = new Array(m.length);
                                for (let i = 0; i < m.length; i++) {
                                    k[i] = m.charCodeAt(i)
                                }
                                h.push(k);
                                break;
                            case f:
                                let m2 = this.getRomanAndFigureString(g);
                                let k1 = new Array(m2.length);
                                for (let i = 0; i < m2.length; i++) {
                                    k1[i] = m2.charCodeAt(i)
                                }
                                h.push(k1);
                                break;
                            case d:
                                let n = this.get8bitByteArray(g);
                                h.push(n);
                                break;
                            case o:
                                let m1 = this.getKanjiString(g);
                                h.push(m1);
                                break
                        }
                    } while (true);
                    return h
                })
            }

            scope.qrcode = qrcode
        }
    }

    let initiated = false;
    let WebQR = {
        video: undefined,
        lockLayer: undefined,
        container: undefined,
        timeOut: undefined,
        stream: undefined,
        qrcode: undefined,
        timer: undefined,
        canvas: undefined, ctx: undefined,

        draw(video: undefined) {

        }
    };

    function addCSSRule(selector, rules, index, sheet) {
        index = index || 0;
        sheet = sheet || styleSheet;
        if ("insertRule" in sheet) {
            sheet.insertRule(selector + "{" + rules + "}", index);
        } else if ("addRule" in sheet) {
            sheet.addRule(selector, rules, index);
        }
    }

    let styleSheet = (function () {
        // Create the <style> tag
        let style = document.createElement("style");

        // Add a media (and/or media query) here if you'd like!
        // style.setAttribute("media", "screen")
        // style.setAttribute("media", "only screen and (max-width : 1024px)")

        // WebKit hack :(
        style.appendChild(document.createTextNode(""));

        // Add the <style> element to the page
        document.head.appendChild(style);

        return style.sheet;
    })();

    function createGlobalComponent(options) {

        options = options || {};
        let container = document.createElement('div');
        let lockLayer = document.createElement('div');
        lockLayer.className = 'QRScanner-lock-layer ' + options.lockLayerClassName;
        container.className = 'QRScanner-container ' + options.className;
        let innerHTML = '<canvas id="QRScanner-canvasEl" width="240" height="200"></canvas>';
        innerHTML += '<video id="QRScanner-videoEl" style="display: none;" width="200" height="200"></video>';
        container.innerHTML = innerHTML;
        lockLayer.innerHTML = '.......';

        addCSSRule('.QRScanner-container',
            `position: fixed;
            left: 0;
            top: 0;
            right: 0;
            bottom: 0;
            box-sizing: border-box;
            width: 260px;
            height: 220px;
            z-index: 999999999;
            margin: auto;
            padding: 10px;
            display: block;
            background-color: ${options.bgColor || '#f0f0f0'};
            box-shadow: ${options.shadow || '0px 0px 10px #000'}
        `,undefined,undefined);
        addCSSRule('.QRScanner-lock-layer',
            `position: fixed;
            left: 0;
            top: 0;
            right: 0;
            bottom: 0;
            box-sizing: border-box;
            width: 100%;
            height: 100%;
            z-index: 999999995;
            margin: auto;
            display: block;
            background-color: #000;
            opacity: .3;
        `,undefined,undefined);

        (options.parent || document.body).appendChild(container);
        (options.lockLayerParent || document.body).appendChild(lockLayer);
        WebQR.lockLayer = lockLayer;
        WebQR.container = container
    }

    MODULE.initiate = function (opts: {
        timeout?: number;
        match?: any;
        onResult: any; onError?: any; onTimeout?: any; }) {
        let options = opts ;
        options.onResult = opts.onResult || function (result) {
            console.info('RESULT ::: ', result);
        },
            options.onError = opts.onError || function (err) {
                console.error('ERR :::: ', err);
            },
            options.onTimeout = opts.onTimeout || function () {
                console.warn('TIMEDOUT');
            };

        if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
            options.onError('Media Devices is not supported!');
            return
        }

        function close(result) {
            if (options.match) {
                if (result && !result.match(options.match)) {
                    return
                }
            }
            window.clearInterval(WebQR.timer);
            window.clearTimeout(WebQR.timeOut);
            WebQR.stream.getTracks()[0].stop();
            WebQR.container.style.display = 'none';
            WebQR.lockLayer.style.display = 'none';
            WebQR.ctx.clearRect(0, 0, 240, 200);

            if (result === false) {
                return
            }
            if (result) {
                options.onResult(result)
            }
        }

        window.clearInterval(WebQR.timer);

        if (!initiated) {
            QRCodeDecoder(WebQR);
            createGlobalComponent(options);
            WebQR.video = document.getElementById('QRScanner-videoEl');
            WebQR.canvas = document.getElementById('QRScanner-canvasEl');
            WebQR.ctx = WebQR.canvas.getContext('2d');

            WebQR.draw = function (image) {
                WebQR.ctx.drawImage(image, 0, 0, 240, 200);
            };

            WebQR.qrcode.callback = function (result) {
                close(result)
            };

            initiated = true;
            document.body.addEventListener('keyup', event => {
                if (event.keyCode === 27) {
                    close(undefined)
                }
            })
        } else {
            WebQR.container.style.display = 'block';
            WebQR.lockLayer.style.display = 'block'
        }

        WebQR.lockLayer.onclick = event => {
            close(undefined)
        };

        if (navigator.mediaDevices && navigator.mediaDevices.getUserMedia) {
            navigator.mediaDevices.getUserMedia({video: true}).then(function (stream) {
                WebQR.video.srcObject = stream; // window.URL.createObjectURL(stream);
                WebQR.video.play();
                WebQR.stream = stream;
                WebQR.timer = window.setInterval(function () {
                    WebQR.draw(WebQR.video);
                    try {
                        WebQR.qrcode.decode()
                    } catch (e) {
                        // it throws when cannot find a qrcode
                        // so, we simply ignore it
                    }
                }, 200);
                WebQR.timeOut = window.setTimeout(function () {
                    options.onTimeout();
                    close(false)
                }, options.timeout || 20000)
            });
        }
    };

    return MODULE;
};
