# NOMA System Simulator

This project implements the design and simulation of NOMA (Non-Orthogonal Multiple Access) communication systems based on MATLAB and Python platforms, used for researching and validating the performance of NOMA technology in wireless communications.

## Project Overview

NOMA (Non-Orthogonal Multiple Access) is an advanced multiple access technology designed to improve spectral efficiency and the number of user connections in wireless communication systems. Unlike traditional orthogonal multiple access (OMA) technology, NOMA allows multiple users to transmit simultaneously on the same resource block (time and frequency) through power domain multiplexing and serial interference cancellation (SIC) technology at the receiver to achieve user separation.

## File Description

- `NOMA.m`: 原始MATLAB实现的NOMA通信系统
- `rayleigh.m`: MATLAB实现的瑞利衰落信道模型
- `noma_system.py`: Python版本的基础实现
- `noma_system_improved.py`: 改进版的Python实现，包含更完整的功能
- `README.md`: 项目说明文档

## Python Version Features

The Python version provides the following improvements:

1. **完整的模块化设计**: 将系统拆分为多个功能模块，便于理解和维护
2. **OFDM调制支持**: 实现了完整的OFDM调制和解调流程
3. **卷积编码和维特比译码**: 包含前向纠错编码功能
4. **QPSK调制**: 支持四相相移键控调制
5. **瑞利衰落信道**: 模拟真实无线信道环境
6. **SIC接收机**: 实现串行干扰消除算法
7. **BER性能分析**: 提供详细的误比特率分析

## Running Instructions

To run the Python version, please ensure the following dependencies are installed:

```bash
pip install numpy matplotlib scipy
```

然后运行：

```bash
python noma_system_improved.py
```

## Technical Features

- 支持多用户叠加编码与功率分配
- 接收端采用串行干扰消除（SIC）技术
- 基于瑞利衰落信道的误码率（BER）性能分析
- 自定义用户数量、调制方式、信噪比范围
- 完整的端到端仿真流程

## Application Scenarios

- 通信工程研究人员
- 高校学生学习NOMA技术
- 无线通信算法开发者
- 5G/5G+通信系统仿真

## System Requirements

- MATLAB R2010a 或更高版本（原版MATLAB代码）
- Python 3.6 或更高版本（Python版本）
- NumPy, Matplotlib, SciPy库（Python版本）

## Notes

- This project is mainly used for academic simulation and teaching demonstrations
- Simulation parameters can be adjusted according to specific requirements
- The Python version is a functional reproduction and improvement of the MATLAB version