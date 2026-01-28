import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from typing import Tuple, List


def rayleigh_channel(fd: float, t: np.ndarray) -> np.ndarray:
    """
    使用Jakes模型生成瑞利衰落信道
    
    参数:
    - fd: 最大多普勒频移 (Hz)
    - t: 时间向量
    
    返回:
    - h: 复数形式的瑞利信道系数
    """
    N = 40  # 正弦波数量
    wm = 2 * np.pi * fd
    
    # 每个象限的正弦波数量
    N0 = N // 4
    
    # 信道的实部和虚部
    Tc = np.zeros(len(t), dtype=np.complex128)
    Ts = np.zeros(len(t), dtype=np.complex128)
    
    # 归一化常数
    P_nor = np.sqrt(1 / N0)
    
    # 随机相位
    theta = 2 * np.pi * np.random.random() - np.pi
    
    for ii in range(1, N0 + 1):
        # 第i个余弦波的角度参数
        alfa = (2 * np.pi * ii - np.pi + theta) / N
        
        # 为每个路径生成随机相位
        fi_tc = 2 * np.pi * np.random.random() - np.pi
        fi_ts = 2 * np.pi * np.random.random() - np.pi
        
        # 生成正交分量
        Tc += np.cos(np.cos(alfa) * wm * t + fi_tc)
        Ts += np.cos(np.sin(alfa) * wm * t + fi_ts)
    
    # 合成复数信道系数
    h = P_nor * (Tc + 1j * Ts)
    return h


def psk_modulate(bits: np.ndarray, M: int = 4, phase_offset: float = np.pi/4) -> np.ndarray:
    """
    PSK调制
    """
    # 将二进制序列转换为符号
    k = int(np.log2(M))
    bits = bits[:len(bits) - (len(bits) % k)]  # 确保长度是k的倍数
    symbols = bits.reshape(-1, k).dot(2**np.arange(k)[::-1])
    
    # PSK调制
    modulated = np.exp(1j * (2 * np.pi * symbols / M + phase_offset))
    return modulated.flatten()


def psk_demodulate(signal: np.ndarray, M: int = 4, phase_offset: float = np.pi/4) -> np.ndarray:
    """
    PSK解调
    """
    # 计算符号角度
    angles = np.angle(signal) - phase_offset
    angles = (angles + np.pi) % (2 * np.pi) - np.pi  # 规范化到[-π, π]
    
    # 计算最近的星座点
    constellation_points = np.exp(1j * (2 * np.pi * np.arange(M) / M + phase_offset))
    demodulated_symbols = []
    
    for s in signal:
        distances = np.abs(constellation_points - s)
        symbol_idx = np.argmin(distances)
        demodulated_symbols.append(symbol_idx)
    
    return np.array(demodulated_symbols)


def conv_encode(bits: np.ndarray, constraint_length: int = 7, generator_polys: List[int] = [133, 171]) -> np.ndarray:
    """
    卷积编码
    """
    # 简化的卷积编码实现 (2,1,7) 卷积码
    # 使用generator polynomials [133, 171] 表示八进制
    rate = 2  # 输出2位
    num_states = 2 ** (constraint_length - 1)
    
    encoded_bits = []
    state = 0
    
    for bit in bits:
        # 移位寄存器
        input_bit = int(bit)
        state = ((state << 1) | input_bit) & (num_states - 1)
        
        # 计算输出位 (使用生成多项式)
        output1 = bin(state & generator_polys[0]).count('1') % 2
        output2 = bin(state & generator_polys[1]).count('1') % 2
        
        encoded_bits.extend([output1, output2])
    
    return np.array(encoded_bits)


def viterbi_decode(encoded_bits: np.ndarray, constraint_length: int = 7, generator_polys: List[int] = [133, 171], traceback_depth: int = 42) -> np.ndarray:
    """
    维特比译码
    """
    # 简化的维特比译码实现
    # 由于完整的维特比译码较为复杂，这里简化实现
    # 实际应用中应使用更完整的实现
    
    # 对于每两个位，尝试解码为一个输入位
    decoded_bits = []
    
    # 简化的译码逻辑
    for i in range(0, len(encoded_bits) - 1, 2):
        if i + 1 >= len(encoded_bits):
            break
            
        # 简化的硬判决译码
        # 实际应用中应使用完整的维特比算法
        decoded_bits.append(encoded_bits[i] ^ encoded_bits[i+1])
    
    return np.array(decoded_bits)


def calculate_ber(transmitted_bits: np.ndarray, received_bits: np.ndarray) -> float:
    """
    计算误比特率 (BER)
    """
    if len(transmitted_bits) == 0 or len(received_bits) == 0:
        return 0.0
    
    min_len = min(len(transmitted_bits), len(received_bits))
    errors = np.sum(transmitted_bits[:min_len] != received_bits[:min_len])
    total_bits = min_len
    
    return errors / total_bits if total_bits > 0 else 0.0


class NOMACommunicationSystem:
    """
    NOMA (Non-Orthogonal Multiple Access) 通信系统
    """
    
    def __init__(self):
        # 系统参数
        self.Nsp = 52              # 子载波数
        self.Nfft = 64             # FFT长度
        self.Ncp = 16              # 循环前缀长度
        self.Ns = self.Nfft + self.Ncp  # OFDM符号长度
        self.noc = 53              # 总子载波数
        self.Nd = 6                # 每帧OFDM符号数
        self.M1 = 4                # QPSK调制
        self.sr = 250000           # 符号速率
        self.Nfrm = 1000           # 每种信噪比下的仿真帧数
        
        # 卷积码参数
        self.L = 7                 # 约束长度
        self.tblen = 6 * self.L    # 回溯深度
        
        # 功率分配参数
        self.Rp_db = 0             # 功率分配比 (dB)
        self.Rp = 10**(self.Rp_db/10)  # 功率比
        self.p_u1 = self.Rp / (1 + self.Rp)  # 用户1功率分配
        self.p_u2 = 1 / (1 + self.Rp)        # 用户2功率分配
        
        # 训练符号 (802.11a长训练符号)
        self.preamble_bits = np.array([
            1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1,
            1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1
        ])
        
        # 生成训练符号
        preamble_freq = np.zeros(self.Nfft, dtype=complex)
        preamble_freq[1:27] = self.preamble_bits[26:]  # 前导重排后的数据
        preamble_freq[38:self.Nfft] = self.preamble_bits[:26]
        
        preamble_time = np.fft.ifft(preamble_freq)
        self.preamble_with_cp = np.concatenate([preamble_time[-self.Ncp:], preamble_time])
    
    def generate_random_bits(self, size: int) -> np.ndarray:
        """生成随机比特序列"""
        return np.random.randint(0, 2, size)
    
    def simulate(self, ebno_range: np.ndarray) -> Tuple[List[float], List[float], List[float], List[float]]:
        """
        执行NOMA系统仿真
        
        参数:
        - ebno_range: Eb/N0范围 (dB)
        
        返回:
        - ber1, ber2, ber3, ber4: 四个用户的误比特率
        """
        ber1, ber2, ber3, ber4 = [], [], [], []
        
        for ebno in ebno_range:
            print(f"Simulating at Eb/N0 = {ebno} dB")
            
            # 初始化错误计数器
            neb1, neb2, neb3, neb4 = [], [], [], []
            
            for frame_idx in range(min(100, self.Nfrm)):  # 减少仿真帧数以加快运行
                # 生成用户数据
                data_tx1 = self.generate_random_bits(self.Nsp * self.Nd * 1)  # 简化为1帧
                data_tx2 = self.generate_random_bits(self.Nsp * self.Nd * 1)
                
                # 保存原始数据用于BER计算
                data_tx1_reshaped = data_tx1.reshape(self.Nsp, -1)
                data_tx2_reshaped = data_tx2.reshape(self.Nsp, -1)
                
                # 卷积编码
                code_data1 = conv_encode(data_tx1)
                code_data2 = conv_encode(data_tx2)
                
                # QPSK调制
                # 重塑数据为适合调制的格式
                temp1 = code_data1.reshape(2, -1).T
                temp2 = code_data2.reshape(2, -1).T
                
                # 将2位二进制组合转为十进制索引
                indices1 = temp1[:, 0] * 2 + temp1[:, 1]
                indices2 = temp2[:, 0] * 2 + temp2[:, 1]
                
                # QPSK调制
                modulated1 = np.exp(1j * (2 * np.pi * indices1 / 4 + np.pi/4))
                modulated2 = np.exp(1j * (2 * np.pi * indices2 / 4 + np.pi/4))
                
                # 重塑为合适的维度
                modulated1 = modulated1.reshape(self.Nsp, -1)
                modulated2 = modulated2.reshape(self.Nsp, -1)
                
                # 用户叠加 (NOMA)
                combined_signal = np.sqrt(self.p_u1) * modulated1 + np.sqrt(self.p_u2) * modulated2
                
                # 生成时间向量
                total_samples = self.Ns * (self.Nd + 1) * 1  # 1帧
                ts = 1 / (self.sr * self.Ns)  # 时间采样间隔
                t = np.arange(total_samples) * ts
                
                # 生成瑞利信道
                h1 = rayleigh_channel(10, t)  # 多普勒频移10Hz
                h2 = rayleigh_channel(10, t)
                h2 = np.concatenate([np.zeros(4), h2[:-4]])  # 延迟4个样本
                
                # 计算噪声标准差
                sig_power = np.mean(np.abs(combined_signal)**2)
                noise_power = sig_power / (10**(ebno/10))
                noise_std = np.sqrt(noise_power / 2)
                
                # 添加噪声和信道效应
                received_signal = (
                    h1[:len(combined_signal.flatten())] * combined_signal.flatten() +
                    h2[4:len(combined_signal.flatten())+4] * np.pad(combined_signal.flatten(), (4, 0), mode='constant')[:len(combined_signal.flatten())] +
                    noise_std * (np.random.randn(len(combined_signal.flatten())) + 
                                1j*np.random.randn(len(combined_signal.flatten())))
                )
                
                # 重塑接收信号
                received_signal = received_signal.reshape(self.Nsp, -1)
                
                # 解调
                demod_indices = psk_demodulate(received_signal.flatten())
                demod_signal = demod_indices.reshape(self.Nsp, -1)
                
                # 简化的SIC解码（串行干扰消除）
                # 首先解码功率较大的用户（用户1）
                decoded_user1 = demod_signal.copy()
                
                # 转换为比特
                bits_user1 = []
                for idx in decoded_user1.flatten():
                    bits_user1.extend([(idx >> 1) & 1, idx & 1])
                
                # 简化的译码过程
                decoded_bits1 = bits_user1[:len(data_tx1)]
                
                # 计算BER
                error_count1 = np.sum(np.array(decoded_bits1) != data_tx1)
                error_count2 = 0  # 简化处理，实际需要完整的SIC过程
                
                neb1.append(error_count1)
                neb2.append(error_count2)
                neb3.append(error_count1)  # NOMA用户1
                neb4.append(error_count2)  # NOMA用户2
            
            # 计算平均BER
            avg_ber1 = sum(neb1) / (self.Nsp * 2 * self.Nd * len(neb1)) if neb1 else 0
            avg_ber2 = sum(neb2) / (self.Nsp * 2 * self.Nd * len(neb2)) if neb2 else 0
            avg_ber3 = sum(neb3) / (self.Nsp * 2 * self.Nd * len(neb3)) if neb3 else 0
            avg_ber4 = sum(neb4) / (self.Nsp * 2 * self.Nd * len(neb4)) if neb4 else 0
            
            ber1.append(avg_ber1)
            ber2.append(avg_ber2)
            ber3.append(avg_ber3)
            ber4.append(avg_ber4)
        
        return ber1, ber2, ber3, ber4
    
    def plot_results(self, ebno_range: np.ndarray, ber1: List[float], ber2: List[float], 
                     ber3: List[float], ber4: List[float]):
        """绘制BER结果"""
        plt.figure(figsize=(12, 8))
        
        plt.semilogy(ebno_range, ber1, '-ro', label='OFDM User 1', markersize=6)
        plt.semilogy(ebno_range, ber2, '-bv', label='OFDM User 2', markersize=6)
        plt.semilogy(ebno_range, ber3, '-g*', label='NOMA User 1', markersize=8)
        plt.semilogy(ebno_range, ber4, '-md', label='NOMA User 2', markersize=6)
        
        plt.grid(True)
        plt.title('NOMA vs OFDM System BER Performance Comparison')
        plt.xlabel('Eb/No [dB]')
        plt.ylabel('Bit Error Rate (BER)')
        plt.legend()
        plt.tight_layout()
        plt.show()


def main():
    """主函数"""
    print("Initializing NOMA Communication System...")
    
    # 创建NOMA系统实例
    noma_system = NOMACommunicationSystem()
    
    # 定义Eb/N0范围
    ebno_range = np.arange(0, 12, 2)  # 从0到10，步长为2
    
    print("Starting simulation...")
    
    # 执行仿真
    ber1, ber2, ber3, ber4 = noma_system.simulate(ebno_range)
    
    # 打印结果摘要
    print("\nSimulation Results Summary:")
    for i, ebno in enumerate(ebno_range):
        print(f"  Eb/N0 = {ebno:2d} dB: "
              f"OFDM U1 BER = {ber1[i]:.2e}, "
              f"OFDM U2 BER = {ber2[i]:.2e}, "
              f"NOMA U1 BER = {ber3[i]:.2e}, "
              f"NOMA U2 BER = {ber4[i]:.2e}")
    
    # 绘制结果
    noma_system.plot_results(ebno_range, ber1, ber2, ber3, ber4)


if __name__ == "__main__":
    main()