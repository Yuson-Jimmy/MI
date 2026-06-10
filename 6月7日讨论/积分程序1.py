import math
import numpy as np
from scipy import integrate
from scipy.special import j1

def compute_M_inf(Rd, Rp, h):
    mu0 = 4 * math.pi * 1e-7 # 真空磁导率（H/m）
    lambda_val = 2e-7 # 假设薄膜磁场穿透深度为200nm
    sigma1 = 0 # 5e4 S/m,  AI说的BSCCO
    omega = 2 * math.pi * 50000 # 文献中驱动频率为50khz
    # h为driven和pickup单匝线圈到薄膜的间距之和，Rd、Rp为线圈的半径，d为薄膜厚度
    d = 6e-9  # 假设薄膜厚度为6nm

    M_inf = 0.0 + 0.0j

    # 预计算公式中不随 q 变化的部分
    lambda_inv_sq = 1.0 / (lambda_val ** 2)
    imag_const = -1j * mu0 * omega * sigma1  # 对应 -i*mu0*omega*sigma1

    def integrand(q):
        # q=0 处的极限处理:
        # J1(x)~x/2, 故分子~q^2; 分母中系数~O(1/q); 整体被积函数~q^3 -> 0
        '''if q < 1e-14:
            return 0.0 + 0.0j'''
                
        # Q^2 = q^2 + lambda^{-2} - i*mu0*omega*sigma1
        Q_sq = q**2 + lambda_inv_sq + imag_const
        Q = np.sqrt(Q_sq)
                
        # 分子: e^{-q*h_{i,j}} * J1(q*R_{d,i}) * J1(q*R_{p,j})
        numerator = np.exp(-q * h) * j1(q * Rd) * j1(q * Rp) * (q ** 2 - Q ** 2) * np.sinh(Q * d) / 2 / Q / q
                
        # 分母: cosh(Q*d) + ((Q^2 + q^2)/(2*q*Q)) * sinh(Q*d)
        coeff = (Q_sq + q**2) / (2.0 * q * Q)
        denominator = np.cosh(Q * d) + coeff * np.sinh(Q * d)
                
        return numerator / denominator
            
    # 分别对实部和虚部进行无穷积分 [0, inf)
    real_integral, _ = integrate.quad(
        lambda q: np.real(integrand(q)), 
        0, np.inf, 
        limit=20000, epsabs=1e-13, epsrel=1e-13
    )
    '''imag_integral, _ = integrate.quad(
        lambda q: np.imag(integrand(q)), 
        0, np.inf, 
        limit=200, epsabs=1e-8, epsrel=1e-8
    )'''
            
    integral_ij = real_integral + 0j
    M_inf += Rd * Rp * integral_ij

    M_inf *= np.pi * mu0
    return M_inf


    
V_driven = 5 # 锁相放大器的震荡输出的有效值为0到5V
R = 1e5 # 用于恒流的电阻为99.4kΩ
I_driven = V_driven / R
M_inf = 0 + 0j # 两线圈之间的互感系数，可为复数
omega = 2 * math.pi * 50000 # 文献中驱动频率为50kHz

#假设是在一个圆柱体上缠绕线圈
L_all = [1e-3, 2e-3, 3e-3, 4e-3] # 假设线圈高度为1、2、3mm
D0_all = [1e-3] # 圆柱体的线圈缠绕部分的直径为0.5到1mm
D = 0.000025 # 铜导线的直径为0.025mm
N1 = 1  # driven线圈缠绕的总圈数
h0_all = [1e-4] #离薄膜最近的线圈距离0.1mm至9mm
N2 = 1

# 假设先缠绕driven线圈
#Rds = [(D0 + 2 * (x - 1) * D) / 2 for x in range(1, N1+1)]
#Rps = [(D0 + (2 * N1 + 1 + 2 * x) * D) / 2 for x in range(N2)]

for L in L_all:
    print("============================================")
    for D0 in D0_all:
        print("******************************************")
        for h0 in h0_all:
            N = int(L / D) #每一层的总匝数
            M_inf = 0
            Rds = [(D0 + 2 * (x - 1) * D) / 2 for x in range(1, N1+1)]
            Rps = [(D0 + (2 * N1 + 1 + 2 * x) * D) / 2 for x in range(N2)]
            for Rd in Rds:
                for Rp in Rps:
                    for n_p in range(N):
                        for n_d in range(N):
                            if n_p < N/2:
                                M_inf += compute_M_inf(Rd, Rp, 2*h0 + D + D*(n_p + n_d))
                            else:
                                M_inf -= compute_M_inf(Rd, Rp, 2*h0 + D + D*(n_p + n_d))
            V_pickup1 = M_inf * omega * I_driven
            print(f'当线圈高度为{1000*L}mm，D0为{1000*D0}mm，h0为{1000*h0}mm 时，\n理想中无限大薄膜时的接收信号：{V_pickup1}')

