# Zoozve: A Strip-Mining-Free RISC-V Vector Extension with Arbitrary Register Grouping Compilation Support (WIP)
Zoozve is a compiler modified based on LLVM, aimed at reducing the additional overhead caused by software strip-mining.

The related work for this paper can be viewed on [here](https://arxiv.org/pdf/2504.15678) and [here](https://arxiv.org/pdf/2504.10832).

For detailed operational instructions of this repository, please refer to [README.pdf](https://github.com/ACELab-SHU/LCTES-25-Artifact/blob/main/README.pdf).

# Installation from Source
The repository is divided into two main parts:

1. Zoozve Compiler: This part includes our custom modifications to LLVM, forming the core of the artifact.

2. Official RISC-V GNU Toolchain: This part contains the official RISC-V GNU toolchain, which is linked to the Zoozve repository to provide necessary functionality. 

To begin, install the necessary dependencies and clone the repository to your local machine using the following commands:
```
sudo apt-get install autoconf automake autotools-dev curl python3 python3-pip libmpc-dev libmpfr-dev libgmp-dev gawk build-essential bison flex texinfo gperf libtool patchutils bc zlib1g-dev libexpat-dev ninja-build git cmake libglib2.0-dev ccache
git clone https://github.com/ACELab-SHU/LCTES-25-Artifact.git
cd LCTES-25-Artifact
git submodule update --init --remote riscv-gnu-toolchain
```

Once the repositories are cloned, navigate to the directory of the RISC-V GNU Toolchain to start the compilation process:
```
cd riscv-gnu-toolchain
./configure --with-arch=rv32im --prefix=/path/to/install --with-abi=ilp32
make -j16
```

After successfully compiling the RISC-V GNU Toolchain, the next step is to compile the Zoozve Compiler：
```
cd ../zoozve
mkdir -p build
cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/path/to/install -DLLVM_TARGETS_TO_BUILD="RISCV" -DLLVM_CCACHE_BUILD=ON -DLLVM_ENABLE_PROJECTS="clang;llvm" -DLLVM_USE_LINKER=gold -DLLVM_DEFAULT_TARGET_TRIPLE="riscv32-unknown-elf" ../llvm && ninja install
```

The process may take some time depending on your system. Once the compilation is complete, you will have successfully set up the Zoozve Compiler environment, ready for testing and evaluation.

# Citation
If you use this software, please cite it as:
```bibtex
@software{zoozve,
 title = {Zoozve},
 author = {Siyi Xu, Limin Jiang, Yintao Liu, Yihao Shen, Yi Shi, Shan Cao, Zhiyuan Jiang},
 note = {https://github.com/ACELab-SHU/LCTES-25-Artifact},
 year = {2024}
}
```
If you wish to cite the related work, please use the following citation format:
```bibtex
@article{xu2025zoozve,
  title={Zoozve: A Strip-Mining-Free {RISC-V} Vector Extension with Arbitrary Register Grouping Compilation Support {(WIP)}},
  author={Xu, Siyi and Jiang, Limin and Liu, Yintao and Shen, Yihao and Shi, Yi and Cao, Shan and Jiang, Zhiyuan},
  journal={arXiv preprint arXiv:2504.15678},
  year={2025}
}

@article{jiang2025unlimited,
  title={Unlimited Vector Processing for Wireless Baseband Based on {RISC-V} Extension},
  author={Jiang, Limin and Shi, Yi and Shen, Yihao and Cao, Shan and Jiang, Zhiyuan and Zhou, Sheng},
  journal={arXiv preprint arXiv:2504.10832},
  year={2025}
}
```
