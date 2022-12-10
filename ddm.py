from math import pi, sqrt, log, ceil, floor, exp, sin
from numpy import log


def ftt_01w(tt, w, err=1e-4):
    """Compute f(t|0,1,w) according to Navarro and Fuss (2009)."""

    # calculate number of terms needed for large t
    if pi * tt * err < 1:  # if error threshold is set low enough
        kl = sqrt(-2 * log(pi * tt * err) / (pi ** 2 * tt))  # bound
        kl = max(kl, 1.0 / (pi * sqrt(tt)))  # ensure boundary conditions met
    else:  # if error threshold set too high
        kl = 1.0 / (pi * sqrt(tt))  # set to boundary condition

    # calculate number of terms needed for small t
    if 2 * sqrt(2 * pi * tt) * err < 1:  # if error threshold is set low enough
        ks = 2 + sqrt(-2 * tt * log(2 * sqrt(2 * pi * tt) * err))  # bound
        ks = max(ks, sqrt(tt) + 1)  # ensure boundary conditions are met
    else:  # if error threshold was set too high
        ks = 2  # minimal kappa for that case

    # compute f(tt|0,1,w)
    p = 0  # initialize density
    if ks < kl:  # if small t is better (i.e., lambda<0) ...
        K = ceil(ks)  # round to smallest integer meeting error
        lower = -floor((K - 1) / 2.0)
        upper = ceil((K - 1) / 2.0)
        for k in range(lower, upper + 1):  # loop over k
            p += (w + 2 * k) * exp(-(pow((w + 2 * k), 2)) / 2 / tt)  # increment sum
        p /= sqrt(2 * pi * pow(tt, 3))  # add constant term

    else:  # if large t is better ...
        K = ceil(kl)  # round to smallest integer meeting error
        for k in range(1, K + 1):
            p += (
                    k * exp(-(pow(k, 2)) * (pi ** 2) * tt / 2) * sin(k * pi * w)
            )  # increment sum
        p *= pi  # add constant term

    return p


def pdf(x, k, B, a, err=1e-4):
    """Compute f(t|v,a,z) according to Navarro and Fuss (2009)."""
    # time must be positive
    # x drift time
    # v drift coefficient
    # B decision boundary
    # w start point, the ration of decision boundary

    if x <= 0:
        return 0

    tt = x / B ** 2  # use normalized time

    p = ftt_01w(tt, w=a, err=err)  # get f(t|0,1,w)

    # convert to f(t|v,a,w)
    xx = exp(-k * B * a - (pow(k, 2)) * x / 2.0) / (pow(B, 2))
    r = p * xx
    # print(x, v, a, w,p,r,xx,tt)
    return r

def ddmpdf(k, a, B, ndt, coh, correct, rt):
    """
    该函数接受七个参数, 其中四个是ddm的参数
    k: drift coefficient
    a: initial bias, 这里表示为和B的比例, 范围是(0,1)
    B: decision boundary
    ndt: non-decision time

    另外三个是一个trial的数据
    coh: coherence
    correct: correct 1 or not 0
    rt: reaction time in secs
    """
    if correct==0: # 错误的trial
        return pdf(x=rt-ndt, k=k*coh, B=B, a=a)
    elif correct==1: # 正确的trial
        return pdf(x=rt-ndt, k=-k*coh, B=B, a=1-a) # 我们反转了k的符号，并且改变了initial bias








