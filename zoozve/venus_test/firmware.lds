OUTPUT_ARCH( "riscv" )
MEMORY
{
	ram   (wxa!ri) : ORIGIN = 0x80000000, LENGTH = 128M
}

SECTIONS {
	.memory : {
		. = 0x00005000;
		*(.init);
		KEEP (*(.init))
		*(.text);
		KEEP (*(.text))
		*(*);
		. = ALIGN(4);
		end = .;
	} >ram

	PROVIDE(_memory_start = ORIGIN(ram));
	PROVIDE(_memory_end = ORIGIN(ram) + LENGTH(ram));
}