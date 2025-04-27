INSTALL_DIR=/home/jim/bin

.PHONY: build_gcc
build_gcc:
	cd riscv-gnu-toolchain && \
		./configure --prefix=$(INSTALL_DIR)/riscv32im --with-arch=rv32im --with-abi=ilp32
	cd riscv-gnu-toolchain && \
		make -j16

.PHONY: build_llvm
build_llvm:
	cd zoozve && mkdir -p build
	cd zoozve/build && cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR)/zoozve_debug -DLLVM_TARGETS_TO_BUILD="RISCV" -DLLVM_CCACHE_BUILD=ON -DLLVM_ENABLE_PROJECTS="clang;llvm" -DLLVM_USE_LINKER=gold -DLLVM_DEFAULT_TARGET_TRIPLE="riscv32-unknown-elf" ../llvm && ninja install

.PHONY: test
test:
	cd zoozve/venus_test && make all LLVMPATH=$(INSTALL_DIR)/zoozve_debug/bin RVPATH=$(INSTALL_DIR)/riscv32im/bin
