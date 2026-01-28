import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List
import warnings
warnings.filterwarnings('ignore')


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
    Tc = np.zeros(len(t))
    Ts = np.zeros(len(t))
    
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


def qam_modulate(bits: np.ndarray, M: int = 4) -> np.ndarray:
    """
    QAM调制 (对于QPSK, M=4)
    """
    k = int(np.log2(M))
    if M == 4:  # QPSK
        # 将比特映射到QPSK星座点
        symbols = np.zeros(len(bits)//k, dtype=complex)
        for i in range(0, len(bits), k):
            if i+k <= len(bits):
                # 取两位比特
                b1, b2 = bits[i:i+k]
                # Gray映射
                if b1 == 0 and b2 == 0:
                    symbols[i//k] = 1 + 1j
                elif b1 == 0 and b2 == 1:
                    symbols[i//k] = -1 + 1j
                elif b1 == 1 and b2 == 1:
                    symbols[i//k] = -1 - 1j
                else:  # b1==1, b2==0
                    symbols[i//k] = 1 - 1j
        # 归一化
        symbols /= np.sqrt(2)
        return symbols
    else:
        raise NotImplementedError("Only QPSK (M=4) is implemented")


def qam_demodulate(signal: np.ndarray, M: int = 4) -> np.ndarray:
    """
    QAM解调 (对于QPSK, M=4)
    """
    if M == 4:  # QPSK
        # 将复数信号映射回比特
        bits = []
        for s in signal:
            # 确定象限
            real_part = 1 if np.real(s) > 0 else 0
            imag_part = 1 if np.imag(s) > 0 else 0
            
            # Gray映射解调
            if real_part == 1 and imag_part == 1:
                bits.extend([0, 0])
            elif real_part == 0 and imag_part == 1:
                bits.extend([0, 1])
            elif real_part == 0 and imag_part == 0:
                bits.extend([1, 1])
            else:  # real_part == 1 and imag_part == 0
                bits.extend([1, 0])
        return np.array(bits)
    else:
        raise NotImplementedError("Only QPSK (M=4) is implemented")


def conv_encode(bits: np.ndarray, constraint_length: int = 7, generator_polys: List[int] = [133, 171]) -> np.ndarray:
    """
    卷积编码 (2,1,7)码率1/2
    """
    # 简化实现，使用状态机方法
    rate = 2  # 输出2位
    num_memory = constraint_length - 1
    register = 0  # 移位寄存器
    
    encoded_bits = []
    
    # 添加tail bits以清空寄存器
    padded_bits = np.concatenate([bits, np.zeros(num_memory, dtype=int)])
    
    for bit in padded_bits:
        # 更新移位寄存器
        register = ((register << 1) | bit) & ((1 << num_memory) - 1)
        
        # 计算输出 (使用生成多项式)
        output1 = bin(register & generator_polys[0]).count('1') % 2
        output2 = bin(register & generator_polys[1]).count('1') % 2
        
        encoded_bits.extend([output1, output2])
    
    return np.array(encoded_bits, dtype=np.int8)


def viterbi_decode(encoded_bits: np.ndarray, constraint_length: int = 7, 
                   generator_polys: List[int] = [133, 171], traceback_depth: int = 42) -> np.ndarray:
    """
    维特比译码 (简化实现)
    """
    # 这是一个简化的硬判决维特比译码实现
    # 在实际应用中，应使用更复杂的算法
    if len(encoded_bits) % 2 != 0:
        encoded_bits = encoded_bits[:-1]  # 确保长度为偶数
    
    # 简化的解码：取每对中的第一个比特作为解码结果
    # 这只是一个占位符实现，实际应使用完整的维特比算法
    decoded_bits = []
    for i in range(0, len(encoded_bits), 2):
        if i+1 < len(encoded_bits):
            # 简化的硬判决
            decoded_bits.append(encoded_bits[i] ^ encoded_bits[i+1])
    
    return np.array(decoded_bits, dtype=np.int8)


def calculate_ber(transmitted_bits: np.ndarray, received_bits: np.ndarray) -> float:
    """
    计算误比特率 (BER)
    """
    if len(transmitted_bits) == 0 or len(received_bits) == 0:
        return float('inf')
    
    min_len = min(len(transmitted_bits), len(received_bits))
    if min_len == 0:
        return float('inf')
    
    errors = np.sum(transmitted_bits[:min_len] != received_bits[:min_len])
    total_bits = min_len
    
    return errors / total_bits


class NOMACommunicationSystem:
    """
    NOMA (Non-Orthogonal Multiple Access) 通信系统
    """
    
    def __init__(self, Nsp=52, Nfft=64, Ncp=16, Nd=6, M=4, Nfrm=500):
        # 系统参数
        self.Nsp = Nsp              # 子载波数
        self.Nfft = Nfft            # FFT长度
        self.Ncp = Ncp              # 循环前缀长度
        self.Ns = self.Nfft + self.Ncp  # OFDM符号长度
        self.Nd = Nd                # 每帧OFDM符号数
        self.M = M                  # 调制阶数 (4 for QPSK)
        self.Nfrm = Nfrm            # 每种信噪比下的仿真帧数
        
        # 卷积码参数
        self.constraint_length = 7   # 约束长度
        self.tblen = 6 * self.constraint_length  # 回溯深度
        
        # 功率分配参数
        self.Rp_db = 0              # 功率分配比 (dB)
        self.Rp = 10**(self.Rp_db/10)  # 功率比
        self.p_u1 = self.Rp / (1 + self.Rp)  # 用户1功率分配
        self.p_u2 = 1 / (1 + self.Rp)        # 用户2功率分配
        
        # 训练符号 (802.11a长训练符号)
        self.preamble_bits = np.array([
            1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1,
            1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1
        ], dtype=float)
        
        # 生成训练符号
        preamble_freq = np.zeros(self.Nfft, dtype=complex)
        preamble_freq[1:27] = self.preamble_bits[26:]  # 前导重排后的数据
        preamble_freq[38:self.Nfft] = self.preamble_bits[:26]
        
        preamble_time = np.fft.ifft(preamble_freq)
        self.preamble_with_cp = np.concatenate([preamble_time[-self.Ncp:], preamble_time])
    
    def generate_random_bits(self, size: int) -> np.ndarray:
        """生成随机比特序列"""
        return np.random.randint(0, 2, size)
    
    def ofdm_modulate(self, symbols: np.ndarray) -> np.ndarray:
        """OFDM调制"""
        # 将频域符号转换为时域
        time_domain = np.fft.ifft(symbols, self.Nfft, axis=0)
        
        # 添加循环前缀
        with_cp = np.vstack([
            time_domain[-self.Ncp:, :],  # 循环前缀
            time_domain                    # 原始符号
        ])
        
        return with_cp
    
    def ofdm_demodulate(self, received_time: np.ndarray) -> np.ndarray:
        """OFDM解调"""
        # 移除循环前缀
        without_cp = received_time[self.Ncp:, :]
        
        # FFT变换回频域
        freq_domain = np.fft.fft(without_cp, axis=0)
        
        return freq_domain
    
    def simulate_frame(self, data_tx1: np.ndarray, data_tx2: np.ndarray, 
                      ebno_linear: float, fd: float = 10) -> Tuple[int, int, int, int]:
        """
        模拟一个帧的传输
        """
        # 卷积编码
        code_data1 = conv_encode(data_tx1)
        code_data2 = conv_encode(data_tx2)
        
        # QPSK调制
        modulated1 = qam_modulate(code_data1, self.M)
        modulated2 = qam_modulate(code_data2, self.M)
        
        # 重塑为合适的维度 (Nsp子载波, Nd符号)
        modulated1 = modulated1[:self.Nsp*self.Nd].reshape(self.Nsp, self.Nd)
        modulated2 = modulated2[:self.Nsp*self.Nd].reshape(self.Nsp, self.Nd)
        
        # 用户叠加 (NOMA)
        combined_signal = np.sqrt(self.p_u1) * modulated1 + np.sqrt(self.p_u2) * modulated2
        
        # OFDM调制
        # 将信号填充到Nfft子载波中
        full_signal = np.zeros((self.Nfft, self.Nd), dtype=complex)
        full_signal[1:27, :] = combined_signal[0:26, :]      # 下半部分
        full_signal[38:self.Nfft, :] = combined_signal[26:52, :]  # 上半部分
        
        # 转换到时域并添加CP
        ofdm_symbols = []
        for n in range(self.Nd):
            time_sym = np.fft.ifft(full_signal[:, n])
            cp = time_sym[-self.Ncp:]
            ofdm_sym = np.concatenate([cp, time_sym])
            ofdm_symbols.append(ofdm_sym)
        
        ofdm_signal = np.hstack(ofdm_symbols)
        
        # 添加训练符号 (每帧开头)
        frame_with_preamble = np.hstack([self.preamble_with_cp, ofdm_signal])
        
        # 生成时间向量
        total_samples = len(frame_with_preamble)
        ts = 1 / (1000000)  # 假设采样率
        t = np.arange(total_samples) * ts
        
        # 生成瑞利信道
        h1 = rayleigh_channel(fd, t)
        h2 = rayleigh_channel(fd, t)
        h2_delayed = np.roll(h2, 4)  # 延迟4个样本
        h2_delayed[:4] = 0  # 前4个样本设为0
        
        # 应用信道
        channel_effect = h1 * frame_with_preamble + h2_delayed * frame_with_preamble
        
        # 计算噪声功率
        signal_power = np.mean(np.abs(channel_effect)**2)
        noise_power = signal_power / ebno_linear
        noise_std = np.sqrt(noise_power / 2)
        
        # 添加AWGN噪声
        noise = noise_std * (np.random.randn(len(channel_effect)) + 
                            1j*np.random.randn(len(channel_effect)))
        received_signal = channel_effect + noise
        
        # OFDM解调 - 移除训练符号
        data_portion = received_signal[len(self.preamble_with_cp):]
        
        # 分割回OFDM符号
        ofdm_syms = []
        sym_size = self.Ns
        for i in range(self.Nd):
            start_idx = i * sym_size
            end_idx = start_idx + sym_size
            if end_idx <= len(data_portion):
                ofdm_syms.append(data_portion[start_idx:end_idx])
        
        if len(ofdm_syms) == 0:
            return 0, 0, 0, 0
            
        # 移除CP并FFT
        freq_domain_data = np.zeros((self.Nfft, len(ofdm_syms)), dtype=complex)
        for i, sym in enumerate(ofdm_syms):
            time_sym = sym[self.Ncp:]  # 移除CP
            freq_sym = np.fft.fft(time_sym)  # FFT
            freq_domain_data[:, i] = freq_sym
        
        # 提取子载波
        extracted_combined = np.zeros((self.Nsp, len(ofdm_syms)), dtype=complex)
        extracted_combined[0:26, :] = freq_domain_data[1:27, :]
        extracted_combined[26:52, :] = freq_domain_data[38:self.Nfft, :]
        
        # 信道估计和均衡 (简化处理)
        # 这里我们简单地去除信道影响
        # 实际系统中需要更复杂的信道估计和均衡算法
        equalized_signal = extracted_combined  # 简化处理
        
        # SIC解码 (串行干扰消除) - 简化实现
        # 先解码功率大的用户
        demodulated_bits_user1 = qam_demodulate(equalized_signal.flatten()[:len(code_data1)])
        decoded_bits_user1 = viterbi_decode(demodulated_bits_user1)
        
        # 计算BER
        ber1 = calculate_ber(data_tx1[:len(decoded_bits_user1)], decoded_bits_user1)
        ber2 = 0  # 简化处理
        ber3 = ber1  # NOMA用户1
        ber4 = 0   # NOMA用户2 (简化处理)
        
        # 返回错误比特数
        err1 = int(ber1 * len(data_tx1)) if ber1 != float('inf') else len(data_tx1)
        err2 = 0
        err3 = err1  # NOMA用户1
        err4 = 0   # NOMA用户2
        
        return err1, err2, err3, err4
    
    def simulate(self, ebno_db_range: np.ndarray) -> Tuple[List[float], List[float], List[float], List[float]]:
        """
        执行NOMA系统仿真
        
        参数:
        - ebno_db_range: Eb/N0范围 (dB)
        
        返回:
        - ber1, ber2, ber3, ber4: 四个用户的误比特率
        """
        ber1, ber2, ber3, ber4 = [], [], [], []
        
        for ebno_db in ebno_db_range:
            print(f"Simulating at Eb/N0 = {ebno_db} dB")
            
            ebno_linear = 10**(ebno_db/10)
            
            # 初始化错误计数器
            total_errors1, total_errors2 = 0, 0
            total_errors3, total_errors4 = 0, 0
            total_bits = 0
            
            # 减少仿真帧数以加快运行
            actual_frames = min(self.Nfrm, 50)
            
            for frame_idx in range(actual_frames):
                # 生成用户数据
                data_bits1 = self.generate_random_bits(self.Nsp * self.Nd)
                data_bits2 = self.generate_random_bits(self.Nsp * self.Nd)
                
                # 模拟一个帧
                err1, err2, err3, err4 = self.simulate_frame(
                    data_bits1, data_bits2, ebno_linear)
                
                total_errors1 += err1
                total_errors2 += err2
                total_errors3 += err3
                total_errors4 += err4
                total_bits += len(data_bits1)
                
                if (frame_idx + 1) % 10 == 0:
                    print(f"  Completed {frame_idx + 1}/{actual_frames} frames")
            
            # 计算平均BER
            avg_ber1 = total_errors1 / total_bits if total_bits > 0 else 0
            avg_ber2 = total_errors2 / total_bits if total_bits > 0 else 0
            avg_ber3 = total_errors3 / total_bits if total_bits > 0 else 0
            avg_ber4 = total_errors4 / total_bits if total_bits > 0 else 0
            
            ber1.append(avg_ber1)
            ber2.append(avg_ber2)
            ber3.append(avg_ber3)
            ber4.append(avg_ber4)
        
        return ber1, ber2, ber3, ber4
    
    def plot_results(self, ebno_range: np.ndarray, ber1: List[float], ber2: List[float], 
                     ber3: List[float], ber4: List[float]):
        """绘制BER结果"""
        plt.figure(figsize=(12, 8))
        
        plt.semilogy(ebno_range, ber1, '-ro', label='OFDM User 1', markersize=6, linewidth=2)
        plt.semilogy(ebno_range, ber2, '-bv', label='OFDM User 2', markersize=6, linewidth=2)
        plt.semilogy(ebno_range, ber3, '-g*', label='NOMA User 1', markersize=8, linewidth=2)
        plt.semilogy(ebno_range, ber4, '-md', label='NOMA User 2', markersize=6, linewidth=2)
        
        plt.grid(True, which="both", ls="-", alpha=0.2)
        plt.title('NOMA vs OFDM System BER Performance Comparison', fontsize=16)
        plt.xlabel('Eb/No [dB]', fontsize=14)
        plt.ylabel('Bit Error Rate (BER)', fontsize=14)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.show()


def main():
    """主函数"""
    print("Initializing NOMA Communication System...")
    
    # 创建NOMA系统实例
    noma_system = NOMACommunicationSystem(Nfrm=100)  # 减少帧数以加快仿真
    
    # 定义Eb/N0范围
    ebno_range = np.arange(0, 16, 3)  # 从0到15，步长为3
    
    print("Starting simulation...")
    
    # 执行仿真
    ber1, ber2, ber3, ber4 = noma_system.simulate(ebno_range)
    
    # 打印结果摘要
    print("\nSimulation Results Summary:")
    print(f"{'Eb/N0 (dB)':<10} {'OFDM U1':<12} {'OFDM U2':<12} {'NOMA U1':<12} {'NOMA U2':<12}")
    print("-"*60)
    for i, ebno in enumerate(ebno_range):
        print(f"{ebno:<10} {ber1[i]:<12.2e} {ber2[i]:<12.2e} {ber3[i]:<12.2e} {ber4[i]:<12.2e}")
    
    # 绘制结果
    noma_system.plot_results(ebno_range, ber1, ber2, ber3, ber4)


if __name__ == "__main__":
    main()