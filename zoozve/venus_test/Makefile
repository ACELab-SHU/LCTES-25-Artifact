
# Change the directory of the compiler here
LLVMPATH=/home/jim/bin/zoozve_debug/bin
RVPATH=/home/jim/bin/riscv32im/bin

INC_DIR?=

CC=$(LLVMPATH)/clang
OPT=$(LLVMPATH)/opt
LLC=$(LLVMPATH)/llc
DUMP=$(LLVMPATH)/llvm-objdump
CPY=$(LLVMPATH)/llvm-objcopy
VENUSROW=2048
VENUSLANE=32
VENUS_VRFADDR=0x80100000
VenusInputStructAddr = 0x80021000
VENUSARCH=--target=riscv32-unknown-elf --gcc-toolchain=$(RVPATH)/../ -march=rv32imzvenus
VENUSCFLAGS=-mllvm --venus -mllvm --venus-nr-row=$(VENUSROW) -mllvm --venus-nr-lane=$(VENUSLANE) -mllvm --venus-task-name=mytask
VENUSIRFLAGS=-S -emit-llvm -Xclang -disable-O0-optnone
CFLAGS=$(VENUSCFLAGS) $(VENUSARCH) $(VENUSIRFLAGS)
VENUSOPTFLAGS= --mattr=+m,+zvenus --venus-nr-row=$(VENUSROW) --venus-nr-lane=$(VENUSLANE) --venus-task-name=mytask --venus-istruct-baseaddr=$(VenusInputStructAddr) -S -passes=mem2reg,venusplit 
VENUSLLCFLAGS= -O3 --mattr=+m,+zvenus --venus-nr-row=$(VENUSROW) --venus-nr-lane=$(VENUSLANE) --venus-vrf-baseaddr=$(VENUS_VRFADDR) -regalloc=basic --enable-venus-alu-anlysis
VENUSDUMPFLAGS=  -d --mattr=+m,+zvenus -M no-aliases

PYTHON=/usr/bin/python3
RM=/usr/bin/rm

TARGET=test
ASM=$(wildcard *.s)

.PHONY: all header ir clean

h2asm: header ir iropt assembly
target: h2asm bin tohex disassembly
all: header ir_target iropt assembly bin_target tohex disassembly

header:
	$(PYTHON) gen_venusbuiltin_h.py
	$(PYTHON) gen_venustype_h.py

ir:
	$(CC) $(CFLAGS) $(TARGET).c -I$(INC_DIR) -o $(TARGET).0.ll
ir_target:
	$(CC) $(CFLAGS) $(TARGET).c -o $(TARGET).0.ll

iropt:
	$(OPT) $(VENUSOPTFLAGS) $(TARGET).0.ll -o $(TARGET).split.ll

assembly:
	$(LLC) $(VENUSLLCFLAGS) --venus-disable-vm $(TARGET).split.ll -o $(TARGET)_before_merge.s
	$(LLC) $(VENUSLLCFLAGS) $(TARGET).split.ll -o $(TARGET).s

assemblynobind:
	$(LLC) $(VENUSLLCFLAGS_BINDREMOVE) $(TARGET).split.ll -o $(TARGET)_nobind.s

bin:
	$(CC) -v $(VENUSARCH) -o $(TARGET) $(ASM) -ffreestanding -nostdlib --std=gnu99 -Wl,-Bstatic,-T,firmware.lds,-Map,firmware.map,--strip-debug -lgcc
	# $(CPY) -O binary $(TARGET) $(TARGET).bin

bin_target:
	$(CC) -v $(VENUSARCH) -o $(TARGET) $(TARGET).s -ffreestanding -nostdlib --std=gnu99 -Wl,-Bstatic,-T,firmware.lds,-Map,firmware.map,--strip-debug -lgcc
	# $(CPY) -O binary $(TARGET) $(TARGET).bin

tohex:
	# $(PYTHON) makehex.py $(TARGET).bin 65536 > $(TARGET).hex

disassembly:
	$(DUMP) $(VENUSDUMPFLAGS) $(TARGET) > $(TARGET)_asm.txt

clean:
	$(RM) venusbuiltin.h venustype.h $(TARGET).0.ll $(TARGET).split.ll $(TARGET).s $(TARGET)
