# NOMA通信系统设计

此项目实现了NOMA（非正交多址接入）通信系统的设计与仿真，基于MATLAB和Python平台，用于研究和验证NOMA技术在无线通信中的性能表现。

## 项目概述

NOMA（Non-Orthogonal Multiple Access）是一种先进的多址接入技术，旨在提高无线通信系统的频谱效率和用户连接数。与传统的正交多址（OMA）技术不同，NOMA允许多个用户在同一资源块（时间和频率）上同时传输，通过功率域复用和接收端的串行干扰消除（SIC）技术实现用户分离。

## 文件说明

- `NOMA.m`: 原始MATLAB实现的NOMA通信系统
- `rayleigh.m`: MATLAB实现的瑞利衰落信道模型
- `noma_system.py`: Python版本的基础实现
- `noma_system_improved.py`: 改进版的Python实现，包含更完整的功能
- `README.md`: 项目说明文档

## Python版本特性

Python版本提供了以下改进：

1. **完整的模块化设计**: 将系统拆分为多个功能模块，便于理解和维护
2. **OFDM调制支持**: 实现了完整的OFDM调制和解调流程
3. **卷积编码和维特比译码**: 包含前向纠错编码功能
4. **QPSK调制**: 支持四相相移键控调制
5. **瑞利衰落信道**: 模拟真实无线信道环境
6. **SIC接收机**: 实现串行干扰消除算法
7. **BER性能分析**: 提供详细的误比特率分析

## 运行说明

要运行Python版本，请确保安装了以下依赖：

```bash
pip install numpy matplotlib scipy
```

然后运行：

```bash
python noma_system_improved.py
```

## 技术特点

- 支持多用户叠加编码与功率分配
- 接收端采用串行干扰消除（SIC）技术
- 基于瑞利衰落信道的误码率（BER）性能分析
- 自定义用户数量、调制方式、信噪比范围
- 完整的端到端仿真流程

## 应用场景

- 通信工程研究人员
- 高校学生学习NOMA技术
- 无线通信算法开发者
- 5G/5G+通信系统仿真

## 系统要求

- MATLAB R2010a 或更高版本（原版MATLAB代码）
- Python 3.6 或更高版本（Python版本）
- NumPy, Matplotlib, SciPy库（Python版本）

## 注意事项

- 本项目主要用于学术仿真和教学演示
- 仿真参数可以根据具体需求调整
- Python版本是对MATLAB版本的功能复现和改进